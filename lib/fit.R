# lib/fit.R : model fitting and prediction with nested tuning.   Requires lib/prep.R (sourced below).
#
# API (stable; read by T3/T4)
#   MODELS : c("lr","lr_spline","ridge_lasso","cart","rf","xgb","mlp")
#   fit_predict(model, train, test, outcome, seed, tune = TRUE, predictors = NULL, val = NULL,
#               inner_k = 5, n_seeds_mlp = 5, ps = NULL, allow_sensitivity_only = FALSE)
#     model     : one of MODELS
#     train/test: data.frames (rows = patients). `outcome` is the 0/1 column NAME (e.g. "y_nonhome") and
#                 must be present and non-NA in train; test needs the predictor columns only.
#     seed      : integer. Controls inner folds (shared across models for the same seed -> paired tuning)
#                 and, offset per model, all model randomness. Same inputs + seed => identical output.
#     tune      : TRUE = small documented search with log-loss (inner stratified `inner_k`-fold CV inside
#                 `train`; or, if `val` is a data.frame, tuning on `val` instead = replication arm).
#                 FALSE = single default configuration (fast; used by tests).
#     predictors: NULL = ps$common (45). Any subset of ps$common is allowed (T4 sensitivities); anything
#                 else stops (leakage guard). allow_sensitivity_only=TRUE additionally admits ps$sensitivity_only.
#     Returns numeric vector of predicted probabilities (length nrow(test), clipped to [1e-6, 1-1e-6]) with attributes:
#       hparams (named list of chosen hyperparameters), tune_score (mean inner log-loss), seconds,
#       converged (logical), notes (character: warnings / fallbacks), n_train, n_events,
#       seed_preds (mlp only: n_test x n_seeds_mlp matrix of per-seed probabilities; the vector is their mean).
#   Preprocessing (median imputation, level encoding, scaling, spline knots) is fitted on `train` only.
#   For inner-CV tuning the unsupervised preprocessing fitted on the outer-train rows is reused (outcome-free,
#   never touches test rows); only the model fits are repeated inside inner folds.
#   No class weights, no resampling, no post hoc calibration. Single-threaded models (parallelise across tasks).
#   .model_defaults / .configs list the search spaces (see model_decisions.md).

if (!exists("prep_fit")) {
  .vv <- get0("V", envir = globalenv(), ifnotfound = Sys.getenv("ASD_V", unset = "."))
  source(file.path(.vv, "lib", "prep.R"))
}

suppressPackageStartupMessages({
  library(glmnet); library(rpart); library(ranger); library(xgboost); library(nnet)
})

MODELS <- c("lr", "lr_spline", "ridge_lasso", "cart", "rf", "xgb", "mlp")
.model_offset <- c(lr = 11L, lr_spline = 23L, ridge_lasso = 37L, cart = 41L, rf = 53L, xgb = 67L, mlp = 79L)

.clip <- function(p, eps = 1e-6) pmin(pmax(p, eps), 1 - eps)
.logloss <- function(y, p, eps = 1e-4) { p <- pmin(pmax(p, eps), 1 - eps); -mean(y * log(p) + (1 - y) * log(1 - p)) }

.collect <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(expr, warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
  list(value = val, warnings = unique(w))
}

# ------------------------------------------------------------------ search spaces
.rf_cfgs <- function(p, tune, seed) {
  if (!tune) return(list(list(mtry = max(1, floor(sqrt(p))), min_node = 20L)))
  g <- expand.grid(mtry = unique(pmax(1, round(c(sqrt(p), 0.2 * p, 0.35 * p, 0.5 * p)))),
                   min_node = c(5L, 15L, 30L, 60L, 120L))
  g <- g[with_seed(seed, sample.int(nrow(g), 6)), ]
  lapply(seq_len(nrow(g)), function(i) list(mtry = g$mtry[i], min_node = g$min_node[i]))
}
.xgb_cfgs <- function(tune, seed) {
  if (!tune) return(list(list(max_depth = 3L, min_child_weight = 5, subsample = 0.8, colsample = 0.8, lambda = 5, nrounds_fixed = 100L)))
  with_seed(seed, lapply(1:8, function(i) list(
    max_depth = sample(2:6, 1), min_child_weight = sample(c(1, 5, 10, 20), 1),
    subsample = sample(c(0.6, 0.8, 1), 1), colsample = sample(c(0.5, 0.75, 1), 1),
    lambda = sample(c(1, 5, 20), 1))))
}
.mlp_cfgs <- function(tune, seed) {
  if (!tune) return(list(list(size = 4L, decay = 0.1)))
  g <- expand.grid(size = c(2L, 4L, 8L, 12L), decay = c(0.01, 0.03, 0.1, 0.3, 1, 3))
  g <- g[with_seed(seed, sample.int(nrow(g), 6)), ]
  lapply(seq_len(nrow(g)), function(i) list(size = g$size[i], decay = g$decay[i]))
}
.cart_cfgs <- function(tune) {
  if (!tune) return(list(list(cp = 0.002, minbucket = 100L)))
  g <- expand.grid(cp = c(0.0003, 0.001, 0.003, 0.01), minbucket = c(30L, 100L, 250L))
  lapply(seq_len(nrow(g)), function(i) list(cp = g$cp[i], minbucket = g$minbucket[i]))
}

# ------------------------------------------------------------------ elementary fit / predict
.fit_cart <- function(X, y, cfg) {
  # Regression tree on the 0/1 outcome (method = "anova"): leaf means are probability estimates. rpart's
  # method = "class" returns a root-only tree for rare outcomes (misclassification-based stopping), so it is not used.
  d <- as.data.frame(X); names(d) <- make.names(names(d), unique = TRUE); d$.y <- as.numeric(y)
  rpart(.y ~ ., data = d, method = "anova",
        control = rpart.control(cp = cfg$cp, minbucket = cfg$minbucket, minsplit = 2 * cfg$minbucket, xval = 0, maxdepth = 12))
}
.pred_cart <- function(fit, X) { d <- as.data.frame(X); names(d) <- make.names(names(d), unique = TRUE)
  as.numeric(predict(fit, d)) }

.fit_rf <- function(X, y, cfg, seed, ntree) {
  colnames(X) <- make.names(colnames(X), unique = TRUE)
  ranger(x = X, y = factor(y, levels = c(0, 1)), num.trees = ntree, mtry = min(cfg$mtry, ncol(X)), min.node.size = cfg$min_node,
         probability = TRUE, num.threads = 1, seed = seed, verbose = FALSE)
}
.pred_rf <- function(fit, X) { colnames(X) <- make.names(colnames(X), unique = TRUE)
  as.numeric(predict(fit, X, num.threads = 1)$predictions[, "1"]) }

.xgb_params <- function(cfg) list(objective = "binary:logistic", eval_metric = "logloss", eta = 0.05,
  max_depth = cfg$max_depth, min_child_weight = cfg$min_child_weight, subsample = cfg$subsample,
  colsample_bytree = cfg$colsample, lambda = cfg$lambda, nthread = 1)
.XGB_MAXR <- 800L
.xgb_es <- function(Xa, ya, Xb, yb, cfg, seed) {   # early stopping on (Xb,yb): used ONLY inside inner folds / val set
  da <- xgb.DMatrix(Xa, label = ya); db <- xgb.DMatrix(Xb, label = yb)
  m <- with_seed(seed, xgb.train(.xgb_params(cfg), da, nrounds = .XGB_MAXR, evals = list(valid = db),
                                 early_stopping_rounds = 30, verbose = 0))
  el <- as.data.frame(attributes(m)$evaluation_log)   # columns: iter, valid_logloss
  i <- which.min(el[[2]])
  c(score = el[[2]][i], iter = el$iter[i])
}
.fit_xgb <- function(X, y, cfg, nrounds, seed) {
  with_seed(seed, xgb.train(.xgb_params(cfg), xgb.DMatrix(X, label = y), nrounds = nrounds, verbose = 0))
}
.pred_xgb <- function(fit, X) as.numeric(predict(fit, X))

.fit_mlp <- function(X, y, cfg, seed) {
  with_seed(seed, nnet(X, y, size = cfg$size, decay = cfg$decay, maxit = 300, MaxNWts = 100000,
                       entropy = TRUE, rang = 0.1, trace = FALSE))
}
.pred_mlp <- function(fit, X) as.numeric(predict(fit, X, type = "raw"))

.fit_glm <- function(X, y) {
  X1 <- cbind(`(Intercept)` = 1, X)
  r <- .collect(suppressMessages(stats::glm.fit(X1, y, family = stats::binomial(), control = stats::glm.control(maxit = 50))))
  f <- r$value
  b <- f$coefficients; b[is.na(b)] <- 0
  list(b = b, converged = isTRUE(f$converged), boundary = isTRUE(f$boundary), rank = f$rank, p = ncol(X1),
       warnings = r$warnings, fitted = f$fitted.values)
}
.pred_glm <- function(fit, X) as.numeric(stats::plogis(cbind(1, X) %*% fit$b))

.ALPHAS <- c(0, 0.5, 1)

# glmnet with penalty at the CROSS-VALIDATED MINIMUM (lambda.min, never 1-SE); alpha chosen among 0 (ridge),
# 0.5 (elastic net), 1 (lasso) by the same inner CV / validation log-loss.
.tune_glmnet <- function(X, y, foldid, tune, valX = NULL, valy = NULL) {
  alphas <- if (tune) .ALPHAS else 0.5
  best <- NULL
  for (a in alphas) {
    if (is.null(valX)) {
      cv <- cv.glmnet(X, y, family = "binomial", alpha = a, foldid = foldid, type.measure = "deviance",
                      standardize = FALSE)
      sc <- min(cv$cvm) / 2; lam <- cv$lambda.min; fit <- cv$glmnet.fit
    } else {
      fit <- glmnet(X, y, family = "binomial", alpha = a, standardize = FALSE)
      pr <- predict(fit, valX, type = "response")
      ll <- apply(pr, 2, function(p) .logloss(valy, p))
      i <- which.min(ll); sc <- ll[i]; lam <- fit$lambda[i]
    }
    if (is.null(best) || sc < best$score) best <- list(alpha = a, lambda = lam, fit = fit, score = sc)
  }
  best
}
.ridge_only <- function(X, y, foldid, valX = NULL, valy = NULL) {
  if (is.null(valX)) {
    cv <- cv.glmnet(X, y, family = "binomial", alpha = 0, foldid = foldid, type.measure = "deviance", standardize = FALSE)
    list(alpha = 0, lambda = cv$lambda.min, fit = cv$glmnet.fit, score = min(cv$cvm) / 2)
  } else {
    fit <- glmnet(X, y, family = "binomial", alpha = 0, standardize = FALSE)
    ll <- apply(predict(fit, valX, type = "response"), 2, function(pp) .logloss(valy, pp)); i <- which.min(ll)
    list(alpha = 0, lambda = fit$lambda[i], fit = fit, score = ll[i])
  }
}
.pred_glmnet <- function(b, X) as.numeric(predict(b$fit, X, s = b$lambda, type = "response"))

# ------------------------------------------------------------------ main entry
fit_predict <- function(model, train, test, outcome, seed, tune = TRUE, predictors = NULL, val = NULL,
                        inner_k = 5, n_seeds_mlp = 5, ps = NULL, allow_sensitivity_only = FALSE) {
  t0 <- proc.time()[["elapsed"]]
  stopifnot(model %in% MODELS)
  if (is.null(ps)) ps <- load_predictor_sets()
  if (is.null(predictors)) predictors <- ps$common
  check_predictors(predictors, ps, allow_sensitivity_only)
  train <- as.data.frame(train); test <- as.data.frame(test)
  y <- train[[outcome]]
  if (is.null(y) || anyNA(y) || !all(y %in% c(0, 1))) stop("train outcome missing, NA, or not 0/1")
  y <- as.numeric(y)
  if (!is.null(val)) { val <- as.data.frame(val); yv <- as.numeric(val[[outcome]]) } else yv <- NULL
  mseed <- as.integer(seed + .model_offset[[model]])
  notes <- character(0); converged <- TRUE; hp <- list(); score <- NA_real_

  cont <- intersect(ps$numeric, predictors)
  prep <- prep_fit(train, predictors, spline_vars = if (model == "lr_spline") cont else NULL)
  use_sp <- model == "lr_spline"
  scale <- model %in% c("lr", "lr_spline", "ridge_lasso", "mlp")
  Xtr <- prep_apply(prep, train, scale = scale, spline = use_sp)
  Xte <- prep_apply(prep, test, scale = scale, spline = use_sp)
  Xva <- if (!is.null(val)) prep_apply(prep, val, scale = scale, spline = use_sp) else NULL
  foldid <- make_folds(y, inner_k, seed)
  splits <- if (!tune) NULL else if (!is.null(val)) list(list(Xa = Xtr, ya = y, Xb = Xva, yb = yv)) else
    lapply(seq_len(inner_k), function(f) list(Xa = Xtr[foldid != f, , drop = FALSE], ya = y[foldid != f],
                                              Xb = Xtr[foldid == f, , drop = FALSE], yb = y[foldid == f]))
  tune_cfg <- function(cfgs, scorer) {
    if (length(cfgs) == 1L) return(list(cfg = cfgs[[1]], score = NA_real_, extra = NA_real_))
    res <- lapply(cfgs, function(cf) { r <- lapply(splits, function(s) scorer(cf, s)); do.call(rbind, r) })
    m <- vapply(res, function(r) mean(r[, 1]), numeric(1))
    i <- which.min(m)
    list(cfg = cfgs[[i]], score = m[i], extra = mean(res[[i]][, 2]))
  }

  if (model == "lr") {
    f <- .fit_glm(Xtr, y)
    p <- .pred_glm(f, Xte)
    converged <- f$converged && !f$boundary && f$rank == f$p
    if (length(f$warnings)) notes <- c(notes, f$warnings)
    if (f$rank < f$p) notes <- c(notes, sprintf("rank deficient (%d of %d)", f$rank, f$p))
    hp <- list(n_params = f$p)
  } else if (model == "lr_spline") {
    f <- .fit_glm(Xtr, y)
    unstable <- !f$converged || f$boundary || f$rank < f$p || max(abs(f$b)) > 15 ||
      any(f$fitted < 1e-8 | f$fitted > 1 - 1e-8)
    if (!unstable) {
      p <- .pred_glm(f, Xte); hp <- list(n_params = f$p, penalty = "none")
    } else {
      b <- .ridge_only(Xtr, y, foldid, Xva, yv)
      p <- .pred_glmnet(b, Xte); score <- b$score
      hp <- list(n_params = f$p, penalty = "ridge_fallback", lambda = b$lambda)
      notes <- c(notes, "unpenalised spline LR unstable (non-convergence, separation, rank or |coef|>15): ridge fallback", f$warnings)
    }
  } else if (model == "ridge_lasso") {
    b <- .tune_glmnet(Xtr, y, foldid, tune, if (!is.null(val)) Xva, yv)
    p <- .pred_glmnet(b, Xte); score <- b$score
    hp <- list(alpha = b$alpha, lambda = b$lambda,
               n_nonzero = sum(as.numeric(coef(b$fit, s = b$lambda))[-1] != 0))
  } else if (model == "cart") {
    cfgs <- .cart_cfgs(tune)
    tr <- if (tune) tune_cfg(cfgs, function(cf, s) c(.logloss(s$yb, .pred_cart(.fit_cart(s$Xa, s$ya, cf), s$Xb)), 0)) else
      list(cfg = cfgs[[1]], score = NA_real_)
    fit <- .fit_cart(Xtr, y, tr$cfg); p <- .pred_cart(fit, Xte)
    hp <- tr$cfg; score <- tr$score
  } else if (model == "rf") {
    cfgs <- .rf_cfgs(ncol(Xtr), tune, seed)
    tr <- if (tune) tune_cfg(cfgs, function(cf, s) c(.logloss(s$yb, .pred_rf(.fit_rf(s$Xa, s$ya, cf, mseed, 300), s$Xb)), 0)) else
      list(cfg = cfgs[[1]], score = NA_real_)
    fit <- .fit_rf(Xtr, y, tr$cfg, mseed, 500); p <- .pred_rf(fit, Xte)
    hp <- tr$cfg; score <- tr$score
  } else if (model == "xgb") {
    cfgs <- .xgb_cfgs(tune, seed)
    if (tune) {
      tr <- tune_cfg(cfgs, function(cf, s) .xgb_es(s$Xa, s$ya, s$Xb, s$yb, cf, mseed))
      nr <- max(10L, as.integer(round(tr$extra)))
    } else { tr <- list(cfg = cfgs[[1]], score = NA_real_); nr <- cfgs[[1]]$nrounds_fixed }
    fit <- .fit_xgb(Xtr, y, tr$cfg, nr, mseed); p <- .pred_xgb(fit, Xte)
    hp <- c(tr$cfg[setdiff(names(tr$cfg), "nrounds_fixed")], list(nrounds = nr)); score <- tr$score
    if (tune && nr >= .XGB_MAXR - 5) notes <- c(notes, "early stopping hit round cap")
  } else if (model == "mlp") {
    cfgs <- .mlp_cfgs(tune, seed)
    tr <- if (tune) tune_cfg(cfgs, function(cf, s) c(.logloss(s$yb, .pred_mlp(.fit_mlp(s$Xa, s$ya, cf, mseed), s$Xb)), 0)) else
      list(cfg = cfgs[[1]], score = NA_real_)
    sp <- matrix(NA_real_, nrow(Xte), n_seeds_mlp); nhit <- 0L
    for (s in seq_len(n_seeds_mlp)) {
      fit <- .fit_mlp(Xtr, y, tr$cfg, mseed + 1000L * s)
      sp[, s] <- .pred_mlp(fit, Xte)
      if (fit$convergence != 0) nhit <- nhit + 1L
    }
    p <- rowMeans(sp)
    hp <- c(tr$cfg, list(n_seeds = n_seeds_mlp, seeds_hit_maxit = nhit)); score <- tr$score
    if (nhit > 0) notes <- c(notes, sprintf("%d of %d final seeds reached maxit=300 without convergence", nhit, n_seeds_mlp))
  }

  out <- .clip(p)
  if (model == "mlp") attr(out, "seed_preds") <- .clip(sp)
  attr(out, "hparams") <- hp
  attr(out, "tune_score") <- score
  attr(out, "converged") <- converged
  attr(out, "notes") <- unique(notes)
  attr(out, "n_train") <- length(y)
  attr(out, "n_events") <- sum(y)
  attr(out, "seconds") <- proc.time()[["elapsed"]] - t0
  out
}
