# test_metrics.R : bootstrap CI covers point, paired differences antisymmetric, metrics vs pROC, determinism, real-OOF consistency.
args_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
V <- if (length(args_file)) dirname(dirname(normalizePath(gsub("~\\+~", " ", args_file[1])))) else "."
suppressPackageStartupMessages({ library(data.table); library(arrow) })
METRICS_NO_RUN <- TRUE
source(file.path(V, "analysis", "03_metrics.R"))
n_pass <- 0L
ok <- function(cond, msg) { if (!isTRUE(cond)) stop("FAILED: ", msg, call. = FALSE); n_pass <<- n_pass + 1L; cat("  pass:", msg, "\n") }

set.seed(11); n <- 600
x <- rnorm(n); y <- rbinom(n, 1, plogis(-1.5 + 0.9 * x))
mk <- function(b, s) plogis(b * x + s * rnorm(n) - 1.5)
Plist <- lapply(1:3, function(r) cbind(lr = mk(0.9, 0.2), xgb = mk(0.7, 0.4), other = mk(0.3, 0.5)))
sf <- stat_fun_factory(c(0.1, 0.2, 0.3))

cat("== metrics vs pROC\n")
p <- Plist[[1]][, "lr"]
ok(abs(auroc(y, p) - as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE, direction = "<")))) < 1e-10, "auroc matches pROC")
ok(abs(sf(y, p)["auroc"] - auroc(y, p)) < 1e-12, "stat function auroc = lib auroc")
ok(abs(sf(y, p)["brier"] - mean((p - y)^2)) < 1e-12, "brier matches definition")
tk <- topk_metrics(y, p, 0.1)
ok(abs(sf(y, p)["top10_spec"] - tk["spec"]) < 1e-12 && abs(sf(y, p)["top10_npv"] - tk["npv"]) < 1e-12, "top-k spec/npv wiring")

cat("== bootstrap (B=120)\n")
b1 <- boot_repeat_avg(y, Plist, sf, B = 120, seed = 5, ref = "lr", cores = 1)
ok(all(b1$ci[stat %in% c("auroc", "brier", "auprc", "cal_slope"), lo <= est & est <= hi]), "bootstrap CI contains point estimate (core metrics)")
ok(b1$ci[stat == "auroc" & model == "lr", est] > 0.6, "point sensible")
pt <- avg_stats(y, Plist, sf)
ok(isTRUE(all.equal(pt["lr", "auroc"], mean(sapply(Plist, function(P) auroc(y, P[, "lr"]))))), "repeat-averaged point = mean of per-repeat AUROC")

cat("== paired difference antisymmetry\n")
b2 <- boot_repeat_avg(y, Plist, sf, B = 120, seed = 5, ref = "xgb", cores = 1)
d1 <- b1$diff[model == "xgb" & stat == "auroc"]; d2 <- b2$diff[model == "lr" & stat == "auroc"]
ok(abs(d1$est_diff + d2$est_diff) < 1e-12, "diff(xgb - lr) = -diff(lr - xgb), point")
ok(abs(d1$lo + d2$hi) < 1e-10 && abs(d1$hi + d2$lo) < 1e-10, "CI bounds mirror")
ok(abs(d1$p_boot - d2$p_boot) < 1e-12, "bootstrap p identical")
ok(d1$lo <= d1$est_diff && d1$est_diff <= d1$hi, "diff CI contains diff point")

cat("== reproducibility\n")
b3 <- boot_repeat_avg(y, Plist, sf, B = 120, seed = 5, ref = "lr", cores = 1)
ok(isTRUE(all.equal(as.data.frame(b1$ci), as.data.frame(b3$ci), check.attributes = FALSE)), "same seed, same CIs")
b4 <- boot_repeat_avg(y, Plist, sf, B = 120, seed = 5, ref = "lr", cores = 4)
ok(isTRUE(all.equal(b1$ci$lo, b4$ci$lo)) && isTRUE(all.equal(b1$ci$hi, b4$ci$hi)), "independent of core count")
b5 <- boot_repeat_avg(y, Plist, sf, B = 120, seed = 6, ref = "lr", cores = 1)
ok(!isTRUE(all.equal(b1$ci$lo, b5$ci$lo)), "different seed changes CI")

cat("== margin class\n")
ok(identical(margin_class(c(0.03, -0.05, -0.01, -0.03), c(0.06, -0.03, 0.01, 0.01)),
             c("better_beyond_margin", "worse_beyond_margin", "within_margin", "inconclusive")), "margin classes")

cat("== real output consistency (skipped if results absent)\n")
f <- file.path(V, "results", "metrics_all.csv")
if (file.exists(f)) {
  all_dt <- fread(f)
  ok(all(all_dt[metric == "auroc", lo <= estimate + 1e-9 & estimate <= hi + 1e-9]), "all AUROC CIs contain estimates")
  ok(all(c("outcome", "tier", "model", "arm", "metric", "estimate", "lo", "hi", "n", "events") %in% names(all_dt)), "tidy columns")
  oof <- as.data.table(arrow::read_parquet(file.path(V, "derived/oof/oof_ssi.parquet")))
  w <- get_wide(oof, "ssi", "cv")
  ref_auc <- mean(sapply(w$Plist, function(P) as.numeric(pROC::auc(pROC::roc(w$y, P[, "xgb"], quiet = TRUE, direction = "<")))))
  ok(abs(all_dt[outcome == "ssi" & model == "xgb" & arm == "cv" & metric == "auroc", estimate] - ref_auc) < 1e-9, "ssi xgb AUROC in results = pROC repeat-average")
  dd <- fread(file.path(V, "results", "diff_vs_lr.csv"))
  a <- all_dt[outcome == "ssi" & arm == "cv" & metric == "auroc"]
  ok(abs(dd[outcome == "ssi" & arm == "cv" & metric == "auroc", diff] - (a[model == "xgb", estimate] - a[model == "ridge_lasso", estimate])) < 1e-9, "diff = difference of tidy estimates")
}
cat(sprintf("\nALL %d CHECKS PASSED\n", n_pass))
