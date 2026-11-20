# test_sens.R : checks for 04_sensitivity.R (era blocks, no train/test overlap, subgroup minimum-event rule, metric engine,
# output integrity). Run: Rscript tests/test_sens.R   (loads helpers only; quick settings, nothing is written to results/)
Sys.setenv(VECLIB_MAXIMUM_THREADS = "1", OMP_NUM_THREADS = "1", ASD_SENS_QUICK = "1", ASD_CORES = "2")
args_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
V <- if (length(args_file)) dirname(dirname(normalizePath(gsub("~\\+~", " ", args_file[1])))) else "."
SENS_NO_RUN <- TRUE
source(file.path(V, "analysis", "04_sensitivity.R"))

n_pass <- 0L
ok <- function(cond, msg) { if (!isTRUE(cond)) stop("FAILED: ", msg, call. = FALSE); n_pass <<- n_pass + 1L; cat("  pass:", msg, "\n") }

dat <- get_data(); d <- dat$d; ps <- dat$ps
cat("== era blocks\n")
ok(!anyNA(d$era) && all(as.character(d$era) %in% ERA_BLOCKS), "every primary-cohort row has one of the four era blocks")
tab <- table(factor(as.character(d$era), levels = ERA_BLOCKS))
ok(sum(tab) == nrow(d) && all(tab > 0), "era blocks partition the cohort (sizes sum to cohort size)")
ok(identical(unname(tapply(d$puf_year, as.character(d$era), min)[ERA_BLOCKS]), c(2012L, 2016L, 2019L, 2022L)) ||
   all(tapply(d$puf_year, as.character(d$era), min)[ERA_BLOCKS] == c(2012, 2016, 2019, 2022)), "era blocks start 2012, 2016, 2019, 2022")
for (ycol in c("y_nonhome", "y_los_p75", "y_readm_unplanned", "y_ssi", "y_cardiac")) {
  dd <- outcome_data(d, ycol); sp <- era_specs(dd)
  ok(length(sp) == 6, paste(ycol, ": four era blocks plus two forward-chaining splits"))
  te_all <- unlist(lapply(sp[1:4], `[[`, "test_idx"))
  ok(!anyDuplicated(te_all) && length(te_all) == nrow(dd), paste(ycol, ": held-out blocks are disjoint and cover every patient once"))
  for (s in sp) {
    ok(length(intersect(s$train_idx, s$test_idx)) == 0 && length(intersect(dd$case_uid[s$train_idx], dd$case_uid[s$test_idx])) == 0,
       paste(ycol, s$tag, ": no train/test overlap (row or case_uid)"))
    if (startsWith(s$tag, "era:")) ok(length(s$train_idx) + length(s$test_idx) == nrow(dd), paste(ycol, s$tag, ": train = all other blocks"))
    if (startsWith(s$tag, "forward")) ok(max(dd$puf_year[s$train_idx]) < min(dd$puf_year[s$test_idx]), paste(ycol, s$tag, ": training years strictly precede test years"))
  }
  ok(!any(dd$puf_year[sp[[6]]$test_idx] == 2021) && !any(dd$puf_year[sp[[6]]$train_idx] == 2021), paste(ycol, ": 2021 schema-break year excluded from forward chaining to 2022-2024"))
}

cat("== hold-out runner reproduces lib run_era_validation (same seeds)\n")
dd <- outcome_data(d, "y_nonhome"); dd <- dd[order(dd$case_uid), ]
a <- run_specs(dd, "y_nonhome", "lr", era_specs(dd)[1:4], ps$common, ps, "test", "test")
b <- run_era_validation(dd, "y_nonhome", "lr", era_col = "era", predictors = ps$common, cores = 1, ps = ps)$oof
a <- a[order(case_uid)]; b <- b[order(case_uid)]
ok(identical(a$case_uid, b$case_uid) && max(abs(a$p - b$p)) < 1e-12, "parallel hold-out runner equals lib run_era_validation for lr")

cat("== subgroup minimum-event rule\n")
ok(identical(subgroup_estimable(c(99, 100, 250, NA)), c(FALSE, TRUE, TRUE, FALSE)), "estimable only with >= 100 events (NA -> not estimable)")
ok(SUBGROUP_MIN_EVENTS == 100L, "threshold is 100 events")
sf <- file.path(P$results, "sens_subgroups.csv")
if (file.exists(sf)) {
  sg <- fread(sf)
  ok(all(sg[estimable == TRUE]$events >= 100), "every reported subgroup estimate has >= 100 events")
  ok(all(is.na(sg[estimable == FALSE]$estimate)) && all(sg[estimable == FALSE]$events < 100), "non-estimable rows carry no estimate and have < 100 events")
  ok(all(sg[subgroup_var != "overall" & estimable == TRUE & metric == "auroc"]$events >= 100), "AUROC rows only for subgroups with >= 100 events")
  ok(any(sg$subgroup == ">=75") && !any(sg$subgroup == ">"), "age >=75 label is kept intact (subset split on the first = only)")
} else cat("  (skipped: results/sens_subgroups.csv not built yet)\n")

cat("== metric engine\n")
set.seed(1); y <- rbinom(400, 1, 0.3); p <- plogis(qlogis(0.3) + 1.2 * (y - 0.3) + rnorm(400))
r1 <- sens_rows(y, list(m = matrix(p, ncol = 1), n = matrix(runif(400), ncol = 1)), "t", "o", "s", ref = "n", B = 40, seed = 3L, cores = 1)
r2 <- sens_rows(y, list(m = matrix(p, ncol = 1), n = matrix(runif(400), ncol = 1)), "t", "o", "s", ref = "n", B = 40, seed = 3L, cores = 1)
ok(abs(r1[model == "m" & metric == "auroc", estimate] - auroc(y, p)) < 1e-12, "AUROC point estimate equals lib auroc")
ok(all(r1[metric == "auroc" & model == "m", c(lo, hi)] > 0.5), "bootstrap interval sensible")
ok(all(names(r1) == TIDY_COLS), "tidy columns in the documented order")
rc <- sens_rows(y, list(m = matrix(p, ncol = 1)), "t", "o", "s", B = 20, seed = 3L, cores = 1)
rc2 <- sens_rows(y, list(m = matrix(p, ncol = 1)), "t", "o", "s", B = 20, seed = 3L, cores = 2)
ok(isTRUE(all.equal(rc$estimate, rc2$estimate)) && isTRUE(all.equal(rc$lo, rc2$lo)), "bootstrap output independent of core count")
yl <- c(rep(1, 60), rep(0, 340)); pl <- runif(400)
rl <- sens_rows(yl, list(m = matrix(pl, ncol = 1)), "t", "o", "s", B = 10, seed = 1L, cores = 1)
ok(all(is.na(rl[metric %in% c("cal_slope", "cal_intercept"), estimate])) && !is.na(rl[metric == "auroc", estimate]), "calibration slope/intercept withheld below 100 events, AUROC kept")
base <- list(m = matrix(p, ncol = 1))
rb <- sens_rows(y, list(m = matrix(p, ncol = 1)), "t", "o", "s", base = base, B = 10, seed = 1L, cores = 1)
ok(abs(rb[metric == "delta_auroc_vs_primary", estimate]) < 1e-12, "paired delta vs primary is zero for identical predictions")

cat("== predictor sets and leakage\n")
ok(length(CPT_VARS) == 9 && "n_cpt_codes" %in% CPT_VARS, "nine CPT-derived variables read from cohort_predictors.csv")
ok(all(CPT_VARS %in% ps$common), "CPT-derived variables are inside the common set (dropping is possible)")
src <- readLines(file.path(V, "analysis", "04_sensitivity.R"))
ok(!any(grepl("allow_sensitivity_only *= *TRUE", src)), "script never admits the so_* sensitivity-only fields")
ok(!any(grepl("y_los_p75|y_nonhome", gsub("OUT <- .*", "", grep("predictors *=|pred *=", src, value = TRUE)))), "no outcome column passed as a predictor")
ok(!any(c("optime", "dt_los_days") %in% ps$common), "primary predictor set contains no operative time or discharge-time variable")

cat("== output integrity (if built)\n")
sa <- file.path(P$results, "sens_all.csv")
if (file.exists(sa)) {
  s <- fread(sa)
  ok(identical(names(s), TIDY_COLS), "sens_all.csv has the documented columns")
  dtm <- s[analysis == "readmission_discharge_time_model" & grepl("^dt_", model)]
  ok(nrow(dtm) > 0 && all(grepl("dischargetime", dtm$outcome)), "discharge-time model rows carry dt_ model names and a separate outcome label")
  prim <- s[analysis == "readmission_discharge_time_model" & outcome == "readm_unplanned"]
  ok(!any(grepl("^dt_", prim$model)), "primary comparator rows never use dt_ names")
  ok(all(s[analysis == "day_of_surgery_addon_optime" & !grepl("primary", subset)]$subset == s[analysis == "day_of_surgery_addon_optime" & !grepl("primary", subset)]$subset) &&
     all(grepl("NOT a primary model", s[analysis == "day_of_surgery_addon_optime" & !grepl("primary", subset), unique(subset)])), "operative-time rows are labelled as a non-primary day-of-surgery add-on")
  ok(all(s[analysis == "era_validation" & outcome == "renal"]$model == "descriptive"), "renal is descriptive only in the era analysis")
  ok(!anyNA(s[metric == "auroc", estimate]), "no missing AUROC estimates")
  ok(all(s[metric == "auroc", estimate] > 0.4 & s[metric == "auroc", estimate] < 1), "AUROC estimates are in a plausible range")
  era <- s[analysis == "era_validation" & metric == "auroc" & grepl("^era:", subset)]
  ok(all(era[, .N, by = .(outcome, model)]$N == 4), "four held-out era AUROCs per outcome and model")
} else cat("  (skipped: results/sens_all.csv not built yet)\n")
cat(sprintf("\nALL %d CHECKS PASSED\n", n_pass))
