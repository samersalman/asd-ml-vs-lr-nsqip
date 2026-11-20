# test_models.R : leakage guard, train-only preprocessing, determinism, OOF structure, metrics unit tests.
# Run: Rscript tests/test_models.R   (from anywhere). Uses a small subset for model runs (patient-level data stays in memory only).
Sys.setenv(VECLIB_MAXIMUM_THREADS = "1", OMP_NUM_THREADS = "1")
args_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
V <- if (length(args_file)) dirname(dirname(normalizePath(gsub("~\\+~", " ", args_file[1])))) else "."
suppressPackageStartupMessages({ library(arrow); library(data.table) })
source(file.path(V, "lib", "prep.R")); source(file.path(V, "lib", "fit.R"))
source(file.path(V, "lib", "cv.R")); source(file.path(V, "lib", "metrics.R"))

n_pass <- 0L
ok <- function(cond, msg) { if (!isTRUE(cond)) stop("FAILED: ", msg, call. = FALSE); n_pass <<- n_pass + 1L; cat("  pass:", msg, "\n") }
near <- function(a, b, tol = 1e-8) isTRUE(all(abs(a - b) <= tol))

dat <- load_modelling_data(V); d <- dat$d; ps <- dat$ps
cat("== leakage guard\n")
ok(check_predictors(ps$common, ps), "common set accepted")
for (bad in c("los_days", "y_nonhome", "death30", "puf_year", "era", "case_uid", "so_weight_loss", "dx_def", "in_primary"))
  ok(inherits(try(check_predictors(c("age", bad), ps), silent = TRUE), "try-error"), paste("rejected:", bad))
ok(!inherits(try(check_predictors(c("age", "so_weight_loss"), ps, allow_sensitivity_only = TRUE), silent = TRUE), "try-error"), "sensitivity-only allowed only by explicit flag")
ok(all(PARSIMONIOUS_PREDICTORS %in% ps$common), "parsimonious set inside common")
ok(all(RESTRICTED_PREDICTORS(ps) %in% ps$common) && length(RESTRICTED_PREDICTORS(ps)) == 29, "restricted set = common minus labs and indicators (29 variables)")
pp <- prep_fit(d, PARSIMONIOUS_PREDICTORS)
ok(ncol(prep_apply(pp, d)) == 10, "parsimonious set has exactly 10 model parameters")
ok(ncol(prep_apply(prep_fit(d, ps$common), d)) == ps$candidate_parameters[["full"]], "full set parameter count matches T1 (62)")
sm <- outcome_data(d, "y_nonhome"); sm <- sm[order(sm$case_uid), ]
ok(inherits(try(fit_predict("lr", sm[1:500, ], sm[501:600, ], "y_nonhome", 1L, predictors = c("age", "los_days")), silent = TRUE), "try-error"), "fit_predict refuses leaking predictor")
ok(inherits(try(fit_predict("lr", sm[1:500, ], sm[501:600, ], "y_nonhome", 1L, predictors = c("age", "y_ssi")), silent = TRUE), "try-error"), "fit_predict refuses another outcome as predictor")

cat("== preprocessing fitted on train only\n")
tr <- sm[1:1500, ]; te <- sm[1501:1800, ]
pr1 <- prep_fit(tr, ps$common, spline_vars = ps$numeric)
te2 <- te; te2$age <- 99; te2$lab_albumin <- NA; te2$bmi <- 60
pr2 <- prep_fit(tr, ps$common, spline_vars = ps$numeric)
ok(identical(pr1, pr2), "fitted object is a pure function of train rows")
ok(near(pr1$med$age, median(tr$age)) && near(pr1$med$lab_albumin, median(tr$lab_albumin, na.rm = TRUE)), "imputer medians equal TRAIN medians")
X1 <- prep_apply(pr1, te); X2 <- prep_apply(pr1, te2)
ok(identical(dim(X1), dim(X2)) && identical(colnames(X1), colnames(X2)), "same columns for any test rows")
ok(near(X1[, "lab_sodium"], X2[, "lab_sodium"]), "changing test rows does not change other test rows' encoding")
alb <- which(is.na(te$lab_albumin))[1]
ok(near(prep_apply(pr1, te, scale = FALSE)[alb, "lab_albumin"], pr1$med$lab_albumin), "missing test lab imputed with train median")
ok(near(colMeans(prep_apply(pr1, tr))[1:5], rep(0, 5), 1e-9), "scaling uses train mean (train columns centred)")
te3 <- te; te3$race_group <- factor("Martian", levels = c(levels(te$race_group), "Martian"))
X3 <- prep_apply(pr1, te3, scale = FALSE); rc <- grep("^race_group=", colnames(X3))
ok(all(X3[, rc] == 0), "unseen factor level -> all-zero dummies, no error")
p_full <- fit_predict("ridge_lasso", tr, te, "y_nonhome", 5L, tune = FALSE)
p_part <- fit_predict("ridge_lasso", tr, te[1:50, ], "y_nonhome", 5L, tune = FALSE)
ok(near(as.numeric(p_full)[1:50], as.numeric(p_part), 1e-10), "predictions for a test row do not depend on other test rows")
te_mod <- te; te_mod$age[1:50] <- 18
ok(near(as.numeric(fit_predict("xgb", tr, te_mod, "y_nonhome", 5L, tune = FALSE))[51:300],
        as.numeric(fit_predict("xgb", tr, te, "y_nonhome", 5L, tune = FALSE))[51:300], 1e-10), "editing some test rows leaves other predictions unchanged (xgb)")

cat("== determinism (small subset, tune = TRUE, all models)\n")
tr <- sm[1:900, ]; te <- sm[901:1200, ]
for (m in MODELS) {
  a <- fit_predict(m, tr, te, "y_nonhome", 11L, tune = TRUE); b <- fit_predict(m, tr, te, "y_nonhome", 11L, tune = TRUE)
  ok(identical(as.numeric(a), as.numeric(b)) && all(a > 0 & a < 1), paste("identical rerun:", m))
  if (m == "mlp") ok(identical(attr(a, "seed_preds"), attr(b, "seed_preds")) && ncol(attr(a, "seed_preds")) == 5, "mlp keeps 5 per-seed prediction columns, reproducible")
}
cat("== nested CV structure and core-independence\n")
ss <- sm[1:600, ]
r1 <- run_nested_cv(ss, "y_nonhome", c("lr", "ridge_lasso"), repeats = 2, k = 5, cores = 1, verbose = FALSE, ps = ps)
r2 <- run_nested_cv(ss, "y_nonhome", c("lr", "ridge_lasso"), repeats = 2, k = 5, cores = 2, verbose = FALSE, ps = ps)
ok(identical(r1$oof[order(model, `repeat`, case_uid)]$p, r2$oof[order(model, `repeat`, case_uid)]$p), "OOF identical for cores = 1 and cores = 2")
cnt <- r1$oof[, .N, by = .(model, `repeat`, case_uid)]
ok(all(cnt$N == 1) && nrow(cnt) == 2 * 2 * 600, "every patient has exactly one out-of-fold prediction per model and repeat")
ok(all(r1$oof[, uniqueN(fold), by = .(model, `repeat`)]$V1 == 5), "5 outer folds per repeat")
fo <- make_outer_folds(ss$y_nonhome, 5, 2, 20261116L)
ok(all(abs(tapply(ss$y_nonhome, fo[, 1], mean) - mean(ss$y_nonhome)) < 0.03), "outer folds are stratified (event rate balanced)")
ok(all(r1$runlog$n_failed == 0), "no failed tasks")

cat("== metrics unit tests\n")
y <- c(0, 0, 1, 1); p <- c(0.1, 0.4, 0.35, 0.8)
ok(near(auroc(y, p), 0.75), "AUROC textbook example = 0.75")
ok(near(auroc(c(0, 1, 0, 1), rep(0.5, 4)), 0.5), "AUROC all ties = 0.5")
set.seed(1); yy <- rbinom(2000, 1, 0.2); pp2 <- plogis(qlogis(0.2) + 1.2 * yy + rnorm(2000))
ok(near(auroc(yy, pp2), as.numeric(pROC::auc(yy, pp2, quiet = TRUE)), 1e-10), "AUROC equals pROC")
ok(near(auprc(c(1, 0), c(0.9, 0.1)), 1) && near(auprc(c(0, 1), c(0.9, 0.1)), 0.5), "AUPRC perfect = 1, worst of two = 0.5")
ok(near(auprc(y, p), 0.5 * 1 + 0.5 * (2 / 3)), "AUPRC small example (precision 1 then 2/3)")
ok(near(brier(y, p), mean((p - y)^2)) && near(brier(y, y), 0), "Brier")
ok(near(brier_scaled(yy, rep(mean(yy), length(yy))), 0, 1e-12), "scaled Brier of the prevalence model = 0")
ok(brier_scaled(yy, pmin(pmax(yy * 0.9 + 0.05, 0), 1)) > 0.8, "scaled Brier high for near-perfect model")
set.seed(2); N <- 200000; lp <- rnorm(N, -1.5, 1); pc <- plogis(lp); yc <- rbinom(N, 1, pc)
ok(abs(cal_slope(yc, pc) - 1) < 0.03 && abs(cal_intercept(yc, pc)) < 0.03, "slope ~ 1 and intercept ~ 0 for calibrated predictions")
ok(abs(cal_slope(yc, plogis(2 * lp)) - 0.5) < 0.02, "slope ~ 0.5 when predictions are overdispersed x2")
ok(abs(cal_intercept(yc, plogis(lp - 0.5)) - 0.5) < 0.03, "intercept ~ +0.5 when predictions are 0.5 too low on the logit scale")
ok(abs(oe_ratio(yc, pc) - 1) < 0.02 && near(oe_ratio(c(1, 1, 0, 0), rep(0.25, 4)), 2), "O:E")
set.seed(3); ok(ici(yc[1:5000], pc[1:5000]) < 0.02 && ici(yc[1:5000], plogis(lp[1:5000] + 1)) > 0.08, "ICI small if calibrated, larger if not")
ok(near(net_benefit(c(1, 1, 0, 0), c(0.9, 0.8, 0.7, 0.1), 0.5), (2 - 1 * 1) / 4), "net benefit by hand (TP 2, FP 1 at t = 0.5, 0.25)")
ok(near(nb_treat_all(c(1, 0, 0, 0), 0.2), 0.25 - 0.75 * 0.25), "treat-all net benefit")
tk <- topk_metrics(c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0), c(.9, .2, .8, .1, .1, .1, .1, .1, .1, .1), 0.2)
ok(tk[["flagged"]] == 2 && near(tk[["sens"]], 0.5) && near(tk[["ppv"]], 0.5) && near(tk[["lift"]], 2.5), "top 20% classification")
am <- all_metrics(yc[1:3000], pc[1:3000], thresholds = c(0.1, 0.2)); ok(all(is.finite(am)), "all_metrics returns finite values")
set.seed(4); yb <- rbinom(1500, 1, 0.3); good <- plogis(qlogis(0.3) + 1.5 * (yb - 0.3) + rnorm(1500, 0, .7)); bad <- runif(1500)
pb <- paired_bootstrap(yb, cbind(lr = good, same = good, worse = bad), stat_fun = function(y, p) c(auroc = auroc(y, p), brier = brier(y, p)), B = 200, seed = 1L)
ok(all(abs(pb$diff[model == "same", est_diff]) < 1e-12) && all(pb$diff[model == "same", lo] == 0), "paired bootstrap: identical model has zero difference")
ok(pb$diff[model == "worse" & stat == "auroc", hi] < 0, "paired bootstrap: clearly worse model has CI below zero")
ok(dim(pb$draws)[1] == 200 && all(c("est", "lo", "hi") %in% names(pb$ci)), "bootstrap draws and CI table shape")
cat(sprintf("\nALL %d TESTS PASSED\n", n_pass))
