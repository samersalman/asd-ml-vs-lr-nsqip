# lib/metrics.R : performance metrics and paired bootstrap helper (no model code; imported by T3 and T4).
#
# API (stable)   y = 0/1 numeric vector, p = predicted probability in [0,1]
#   auroc(y, p)                    rank-based (Mann-Whitney), ties = 0.5
#   auprc(y, p)                    average precision (step-wise, ties grouped); compare with prevalence mean(y)
#   brier(y, p); brier_scaled(y, p)   scaled = 1 - Brier / (prev * (1 - prev))
#   logloss(y, p, eps = 1e-6)
#   cal_slope(y, p); cal_intercept(y, p)
#         slope = coefficient of logit(p) in logistic regression of y on logit(p);
#         intercept = calibration-in-the-large: intercept of glm(y ~ offset(logit(p))), i.e. slope fixed at 1
#   oe_ratio(y, p)                 sum(y) / sum(p)
#   ici(y, p, span = 0.75)         Austin-Steyerberg integrated calibration index: mean |loess(y ~ p) - p|; NA if loess fails
#   net_benefit(y, p, t)           TP/n - FP/n * t/(1-t) at threshold t (vectorised over t);  nb_treat_all(y, t)
#   topk_metrics(y, p, frac)       flag the top `frac` of predicted risk (1 - tie-free, by rank): returns named numeric
#                                  (flagged, sens, spec, ppv, npv, lift = ppv/prevalence)
#   all_metrics(y, p, thresholds = NULL, topk = c(0.1, 0.2))  -> named numeric vector of all of the above
#   paired_bootstrap(y, P, stat_fun = all_metrics, B = 1000, seed, ref = "lr", cluster = NULL, cores = 1, strat = FALSE)
#         y : length-n 0/1;  P : n x M numeric matrix of predictions, columns named by model (SAME patients, so
#         every model is recomputed on the SAME resample); stat_fun(y, p) -> named numeric vector.
#         cluster : optional length-n vector (e.g. case_uid) to resample whole clusters; strat: resample within outcome class.
#         Returns list(point = M x S matrix, draws = B x M x S array, ci = data.table(model, stat, est, lo, hi),
#                      diff = data.table(model, stat, est_diff, lo, hi, p_boot) of model - ref (percentile 95% CIs;
#                      p_boot = two-sided bootstrap proportion))
#   oof_wide(oof, outcome_name, repeat_, arm_ = "cv")   long OOF table (columns outcome, repeat, arm, model, case_uid, y, p) -> list(y, P, case_uid) with one column per model

suppressPackageStartupMessages(library(data.table))

.cl <- function(p, eps = 1e-6) pmin(pmax(p, eps), 1 - eps)

auroc <- function(y, p) {
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

auprc <- function(y, p) {
  n1 <- sum(y == 1); if (n1 == 0) return(NA_real_)
  o <- order(p, decreasing = TRUE)
  ys <- y[o]; ps <- p[o]
  tp <- cumsum(ys); fp <- cumsum(1 - ys)
  last <- c(ps[-length(ps)] != ps[-1], TRUE)       # end of each tie group
  tp <- tp[last]; fp <- fp[last]
  prec <- tp / (tp + fp); rec <- tp / n1
  sum(diff(c(0, rec)) * prec)
}

brier <- function(y, p) mean((p - y)^2)
brier_scaled <- function(y, p) { pr <- mean(y); 1 - brier(y, p) / (pr * (1 - pr)) }
logloss <- function(y, p, eps = 1e-6) { p <- .cl(p, eps); -mean(y * log(p) + (1 - y) * log(1 - p)) }

cal_slope <- function(y, p) {
  lp <- stats::qlogis(.cl(p))
  f <- tryCatch(suppressWarnings(stats::glm.fit(cbind(1, lp), y, family = stats::binomial())), error = function(e) NULL)
  if (is.null(f)) NA_real_ else unname(f$coefficients[2])
}
cal_intercept <- function(y, p) {
  lp <- stats::qlogis(.cl(p))
  f <- tryCatch(suppressWarnings(stats::glm.fit(matrix(1, length(y), 1), y, offset = lp, family = stats::binomial())),
                error = function(e) NULL)
  if (is.null(f)) NA_real_ else unname(f$coefficients[1])
}
oe_ratio <- function(y, p) sum(y) / sum(p)

ici <- function(y, p, span = 0.75) {
  if (length(unique(p)) < 5) return(NA_real_)
  fit <- tryCatch(suppressWarnings(stats::loess(y ~ p, span = span, degree = 2)), error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  pc <- suppressWarnings(stats::predict(fit))
  mean(abs(pc - p))
}

net_benefit <- function(y, p, t) {
  n <- length(y)
  vapply(t, function(tt) { fl <- p >= tt; (sum(fl & y == 1) - sum(fl & y == 0) * tt / (1 - tt)) / n }, numeric(1))
}
nb_treat_all <- function(y, t) { pr <- mean(y); pr - (1 - pr) * t / (1 - t) }

topk_metrics <- function(y, p, frac) {
  n <- length(y); k <- max(1L, round(frac * n))
  fl <- logical(n); fl[order(p, decreasing = TRUE, method = "radix")[seq_len(k)]] <- TRUE
  tp <- sum(fl & y == 1); fp <- sum(fl & y == 0); fn <- sum(!fl & y == 1); tn <- sum(!fl & y == 0)
  c(flagged = k, sens = tp / (tp + fn), spec = tn / (tn + fp), ppv = tp / (tp + fp), npv = tn / (tn + fn),
    lift = (tp / (tp + fp)) / mean(y))
}

all_metrics <- function(y, p, thresholds = NULL, topk = c(0.1, 0.2)) {
  out <- c(auroc = auroc(y, p), auprc = auprc(y, p), prevalence = mean(y), brier = brier(y, p),
           brier_scaled = brier_scaled(y, p), logloss = logloss(y, p), cal_slope = cal_slope(y, p),
           cal_intercept = cal_intercept(y, p), oe_ratio = oe_ratio(y, p), ici = ici(y, p))
  if (length(thresholds)) out <- c(out, stats::setNames(net_benefit(y, p, thresholds), paste0("nb_", thresholds)))
  for (f in topk) {
    m <- topk_metrics(y, p, f)
    out <- c(out, stats::setNames(m[c("sens", "ppv", "lift")], paste0("top", round(100 * f), "_", c("sens", "ppv", "lift"))))
  }
  out
}

paired_bootstrap <- function(y, P, stat_fun = all_metrics, B = 1000, seed = 1L, ref = "lr", cluster = NULL,
                             cores = 1, strat = FALSE) {
  P <- as.matrix(P); n <- length(y); M <- ncol(P); mods <- colnames(P)
  stopifnot(nrow(P) == n, !is.null(mods))
  point <- t(vapply(seq_len(M), function(m) stat_fun(y, P[, m]), stat_fun(y, P[, 1])))
  rownames(point) <- mods; S <- ncol(point); snames <- colnames(point)
  draw_idx <- function(b) {
    set.seed(seed + b)
    if (!is.null(cluster)) {
      u <- unique(cluster); s <- sample(u, length(u), replace = TRUE)
      split_idx <- split(seq_len(n), cluster)
      unlist(split_idx[as.character(s)], use.names = FALSE)
    } else if (strat) {
      i1 <- which(y == 1); i0 <- which(y == 0)
      c(sample(i1, length(i1), replace = TRUE), sample(i0, length(i0), replace = TRUE))
    } else sample.int(n, n, replace = TRUE)
  }
  one <- function(b) {
    ix <- draw_idx(b); yy <- y[ix]
    t(vapply(seq_len(M), function(m) stat_fun(yy, P[ix, m]), numeric(S)))   # M x S
  }
  res <- if (cores > 1) parallel::mclapply(seq_len(B), one, mc.cores = cores, mc.preschedule = TRUE) else lapply(seq_len(B), one)
  draws <- array(NA_real_, dim = c(B, M, S), dimnames = list(NULL, mods, snames))
  for (b in seq_len(B)) draws[b, , ] <- res[[b]]
  q <- function(x) stats::quantile(x, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  ci <- rbindlist(lapply(seq_len(M), function(m) rbindlist(lapply(seq_len(S), function(s) {
    qq <- q(draws[, m, s]); data.table(model = mods[m], stat = snames[s], est = point[m, s], lo = qq[1], hi = qq[2]) }))))
  diff <- NULL
  if (!is.null(ref) && ref %in% mods) {
    diff <- rbindlist(lapply(setdiff(mods, ref), function(m) rbindlist(lapply(seq_len(S), function(s) {
      d <- draws[, m, s] - draws[, ref, s]; qq <- q(d)
      pb <- 2 * min(mean(d <= 0, na.rm = TRUE), mean(d >= 0, na.rm = TRUE))
      data.table(model = m, stat = snames[s], est_diff = point[m, s] - point[ref, s], lo = qq[1], hi = qq[2],
                 p_boot = min(1, pb)) }))))
  }
  list(point = point, draws = draws, ci = ci, diff = diff)
}

oof_wide <- function(oof, outcome_name, repeat_, arm_ = "cv") {
  o <- data.table::as.data.table(oof)
  o <- o[outcome == outcome_name & `repeat` == repeat_ & arm == arm_]
  w <- data.table::dcast(o, case_uid + y ~ model, value.var = "p")
  list(y = w$y, P = as.matrix(w[, setdiff(names(w), c("case_uid", "y")), with = FALSE]), case_uid = w$case_uid)
}
