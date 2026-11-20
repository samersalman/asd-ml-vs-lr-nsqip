# 04_sensitivity.R : sensitivity, era, cohort, predictor-set, recalibration and subgroup analyses, instability (T4).
# Inputs : derived/analysis_ready.parquet, derived/predictor_sets.rds, derived/oof/oof_<outcome>.parquet (READ ONLY),
#          results/cohort_predictors.csv, results/events_by_era.csv, read-only NSQIP DuckDB (operative time, in-hospital events only)
# Outputs: results/sens_all.csv, sens_era.csv, sens_subgroups.csv, derived/sens/sens_instability.csv (point-level, not in results/), sens_instability_bins.csv, sens_runlog.csv, sens_headline.md,
#          derived/sens/* (patient-level or cached, never commit)
# Environment knobs: ASD_CORES, ASD_REPEATS (nested CV repeats for refits, default 5), ASD_SENS_RESUME=1 (reuse cached
#          sections and task files), ASD_SENS_SECTIONS (comma list of: era,cohort,predsets,mi,optime,recal,periods,subgroups,instability),
#          ASD_SENS_QUICK=1 (tiny settings, writes to a scratch folder, never touches results/), ASD_BOOT_INSTAB (default 200)
# Tidy row layout: analysis, outcome, model, subset, metric, estimate, lo, hi, n, events
# Sourced by run_all.R (V defined) or run directly: Rscript analysis/04_sensitivity.R. Helpers are testable:
# set SENS_NO_RUN <- TRUE before sourcing to load functions only.

Sys.setenv(VECLIB_MAXIMUM_THREADS = "1", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
if (!exists("V")) V <- Sys.getenv("ASD_V", unset = ".")
source(file.path(V, "analysis", "00_setup.R"))
suppressPackageStartupMessages(library(arrow))
source(file.path(V, "lib", "prep.R")); source(file.path(V, "lib", "fit.R"))
source(file.path(V, "lib", "cv.R")); source(file.path(V, "lib", "metrics.R"))

QUICK   <- Sys.getenv("ASD_SENS_QUICK", unset = "0") == "1"
CORES   <- as.integer(Sys.getenv("ASD_CORES", unset = max(1, parallel::detectCores() - 4)))
REPEATS <- as.integer(Sys.getenv("ASD_REPEATS", unset = if (QUICK) "1" else "3"))   # refits use repeats 1..REPEATS of the primary design (3 by default for runtime; primary used 5)
RESUME  <- Sys.getenv("ASD_SENS_RESUME", unset = "0") == "1"
B_CI    <- if (QUICK) 20L else 200L
B_ERA   <- if (QUICK) 20L else 500L
B_INST  <- as.integer(Sys.getenv("ASD_BOOT_INSTAB", unset = if (QUICK) "6" else "100"))
SECTIONS <- strsplit(Sys.getenv("ASD_SENS_SECTIONS", unset = "era,cohort,predsets,mi,optime,recal,periods,subgroups,instability"), ",")[[1]]
OOF_DIR  <- file.path(P$derived, "oof")
SCRATCH  <- Sys.getenv("ASD_SENS_SCRATCH", unset = tempdir())
SENS_DIR <- if (QUICK) file.path(SCRATCH, "sens_quick") else file.path(P$derived, "sens")
RES_DIR  <- if (QUICK) file.path(SCRATCH, "sens_quick_results") else P$results
TASK_DIR <- file.path(SENS_DIR, "tasks")
dir.create(SENS_DIR, showWarnings = FALSE, recursive = TRUE); dir.create(RES_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------- definitions
OUT <- list(  # label -> outcome column
  nonhome = "y_nonhome", los_p75 = "y_los_p75", readm_unplanned = "y_readm_unplanned", ssi = "y_ssi",
  cardiac = "y_cardiac", los_ge9 = "y_los_ge9")
LAB <- function(ycol) names(OUT)[match(ycol, unlist(OUT))]
# Prespecified decision-curve threshold used for subgroup net benefit: the nearest simple value to the outcome incidence
# in the primary cohort (nonhome 32.2%, LOS>=8 25.8%, readmission 7.3%, SSI 3.6%, cardiac 1.4%).
THR <- c(nonhome = 0.30, los_p75 = 0.25, readm_unplanned = 0.07, ssi = 0.035, cardiac = 0.015)
ERA_BLOCKS <- c("2012-2015", "2016-2018", "2019-2021", "2022-2024")
MIN_EVENTS <- 100L
SUBGROUP_MIN_EVENTS <- MIN_EVENTS
subgroup_estimable <- function(events, min_events = SUBGROUP_MIN_EVENTS) !is.na(events) & events >= min_events

TIDY_COLS <- c("analysis", "outcome", "model", "subset", "metric", "estimate", "lo", "hi", "n", "events")
tidy_empty <- function() data.table(analysis = character(), outcome = character(), model = character(), subset = character(),
                                    metric = character(), estimate = numeric(), lo = numeric(), hi = numeric(), n = integer(), events = integer())

# ---------------------------------------------------------------- metric engine
.METRICS <- c("auroc", "auprc", "brier", "brier_scaled", "cal_slope", "cal_intercept", "net_benefit")
.stat_one <- function(y, p, nb_t, do_cal) {
  c(auroc = auroc(y, p), auprc = auprc(y, p), brier = brier(y, p), brier_scaled = brier_scaled(y, p),
    cal_slope = if (do_cal) cal_slope(y, p) else NA_real_, cal_intercept = if (do_cal) cal_intercept(y, p) else NA_real_,
    net_benefit = if (is.null(nb_t)) NA_real_ else net_benefit(y, p, nb_t))
}
# P: n x R matrix (R = repeats, 1 for hold-out); statistic = mean over the R columns, bootstrap resamples patients.
.stat_mat <- function(y, P, nb_t, do_cal) {
  P <- as.matrix(P)
  rowMeans(vapply(seq_len(ncol(P)), function(r) .stat_one(y, P[, r], nb_t, do_cal), numeric(length(.METRICS))))
}

# sens_rows : tidy performance rows with paired percentile bootstrap CIs.
#   y        : 0/1 vector;  Plist : named list of n x R prediction matrices on the SAME patients (names = model labels)
#   ref      : name of the reference model in Plist for paired differences ("delta_<m>_vs_<ref>" for auroc, auprc, brier)
#   base     : optional named list of matrices (same patients) of the PRIMARY predictors; adds paired "delta_<m>_vs_primary"
#              for auroc, auprc, brier and cal_slope for every name present in both Plist and base
#   nb_t     : decision-curve threshold; adds net_benefit (and a "treat_all" row)
#   cal_min_events : slope/intercept reported only if events >= this (otherwise NA, 'not estimable')
sens_rows <- function(y, Plist, analysis, outcome, subset, ref = NULL, base = NULL, nb_t = NULL, B = B_CI, seed = 1L,
                      cal_min_events = MIN_EVENTS, cores = CORES, keep_draws = FALSE) {
  Plist <- lapply(Plist, as.matrix)
  n <- length(y); ev <- sum(y); do_cal <- ev >= cal_min_events
  stopifnot(all(vapply(Plist, nrow, 1L) == n), !anyNA(y), !any(vapply(Plist, anyNA, TRUE)))
  nm <- names(Plist)
  allP <- Plist
  if (!is.null(base)) { bn <- intersect(nm, names(base)); for (m in bn) allP[[paste0("..base..", m)]] <- as.matrix(base[[m]]) }
  an <- names(allP); M <- length(an); S <- length(.METRICS)
  one <- function(ix) {
    yy <- y[ix]
    out <- vapply(an, function(m) .stat_mat(yy, allP[[m]][ix, , drop = FALSE], nb_t, do_cal), numeric(S))   # S x M
    c(out, if (!is.null(nb_t)) nb_treat_all(yy, nb_t) else NA_real_)
  }
  point <- one(seq_len(n))
  draw <- function(b) { set.seed(seed + b); one(sample.int(n, n, replace = TRUE)) }
  res <- if (cores > 1) parallel::mclapply(seq_len(B), draw, mc.cores = cores, mc.preschedule = TRUE) else lapply(seq_len(B), draw)
  bad <- vapply(res, function(r) !is.numeric(r) || length(r) != S * M + 1L, TRUE)
  if (any(bad)) stop("bootstrap worker failed: ", paste(utils::head(which(bad)), collapse = ","))
  D <- simplify2array(res)                                    # (S*M+1) x B
  dim3 <- array(D[seq_len(S * M), ], dim = c(S, M, B), dimnames = list(.METRICS, an, NULL))
  pt <- matrix(point[seq_len(S * M)], S, M, dimnames = list(.METRICS, an))
  q <- function(x) stats::quantile(x, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  rows <- list()
  add <- function(model, metric, est, dr) {
    qq <- if (all(is.na(dr))) c(NA_real_, NA_real_) else q(dr)
    rows[[length(rows) + 1L]] <<- data.table(analysis = analysis, outcome = outcome, model = model, subset = subset, metric = metric,
                                             estimate = est, lo = qq[1], hi = qq[2], n = n, events = ev)
  }
  for (m in nm) for (s in .METRICS) {
    if (s == "net_benefit" && is.null(nb_t)) next
    add(m, s, pt[s, m], dim3[s, m, ])
  }
  if (!is.null(nb_t)) add("treat_all", "net_benefit", point[S * M + 1L], D[S * M + 1L, ])
  if (!is.null(ref) && ref %in% nm) for (m in setdiff(nm, ref)) for (s in c("auroc", "auprc", "brier"))
    add(m, paste0("delta_", s, "_vs_", ref), pt[s, m] - pt[s, ref], dim3[s, m, ] - dim3[s, ref, ])
  if (!is.null(base)) for (m in intersect(nm, names(base))) for (s in c("auroc", "auprc", "brier", "cal_slope")) {
    bm <- paste0("..base..", m)
    add(m, paste0("delta_", s, "_vs_primary"), pt[s, m] - pt[s, bm], dim3[s, m, ] - dim3[s, bm, ])
  }
  out <- rbindlist(rows)
  # slopes / intercepts: not estimable below the event floor
  if (!do_cal) out[metric %in% c("cal_slope", "cal_intercept") | grepl("cal_slope", metric), `:=`(estimate = NA_real_, lo = NA_real_, hi = NA_real_)]
  out[, `:=`(n = as.integer(n), events = as.integer(events))]
  if (keep_draws) attr(out, "draws") <- dim3
  out
}

# ---------------------------------------------------------------- OOF helpers
# long OOF (outcome, model, arm, repeat, fold, case_uid, y, p) -> list(case_uid, y, P = named list of n x R matrices), sorted by case_uid
oof_lists <- function(oof, arm_ = "cv", models = NULL) {
  o <- as.data.table(oof)[arm == arm_]
  if (!is.null(models)) o <- o[model %in% models]
  res <- list(); ids <- NULL; yy <- NULL
  for (m in unique(o$model)) {
    w <- data.table::dcast(o[model == m], case_uid + y ~ `repeat`, value.var = "p")
    setorder(w, case_uid)
    if (is.null(ids)) { ids <- w$case_uid; yy <- w$y } else stopifnot(identical(ids, w$case_uid))
    res[[m]] <- as.matrix(w[, setdiff(names(w), c("case_uid", "y")), with = FALSE])
  }
  list(case_uid = ids, y = yy, P = res)
}
subset_lists <- function(L, keep) list(case_uid = L$case_uid[keep], y = L$y[keep], P = lapply(L$P, function(m) m[keep, , drop = FALSE]))
# primary OOF (read only). Refit comparisons use primary repeats 1..REPEATS (same outer folds and seeds as the refits); `full = TRUE` keeps all 5 repeats
primary_oof <- function(lab, full = FALSE) { o <- read_oof(OOF_DIR, lab)$oof; if (full) o else o[arm != "cv" | `repeat` <= REPEATS] }
.pref <- function(models) if ("lr" %in% models) "lr" else if ("ridge_lasso" %in% models) "ridge_lasso" else models[1]

.RUNLOG <- list()
log_run <- function(section, label, res, wall) {
  rl <- res$runlog
  .RUNLOG[[length(.RUNLOG) + 1L]] <<- data.table(section = section, label = label, models = paste(unique(rl$model), collapse = "+"),
    n_tasks = sum(rl$n_tasks), n_failed = sum(rl$n_failed), n_not_converged = sum(rl$n_not_converged),
    n_with_notes = sum(rl$n_with_notes), cpu_minutes = round(sum(rl$cpu_seconds) / 60, 2), wall_minutes = round(wall / 60, 2))
}
cv_run <- function(section, label, d, ycol, models, predictors, ps) {
  t0 <- proc.time()[["elapsed"]]
  res <- run_nested_cv(d, ycol, models, predictors = predictors, repeats = REPEATS, k = 5, inner_k = 5, cores = CORES,
                       task_dir = TASK_DIR, resume = RESUME, ps = ps, label = label, verbose = FALSE)
  log_run(section, label, res, proc.time()[["elapsed"]] - t0)
  if (anyNA(res$oof$p)) stop("failed tasks in ", label, ": ", paste(unique(res$runlog[n_failed > 0]$model), collapse = ","))
  message(sprintf("    [%s] %s done (%.1f min wall)", section, label, (proc.time()[["elapsed"]] - t0) / 60))
  res
}

# parallel hold-out runner (identical seeds and calls as lib run_holdout / run_era_validation; flat over blocks x models)
run_specs <- function(d, ycol, models, specs, predictors, ps, section = "era", label = ycol) {
  d <- as.data.frame(d)
  tasks <- list()
  for (m in .heavy_first(models)) for (s in specs) tasks[[length(tasks) + 1L]] <- c(list(model = m), s)
  worker <- function(t) {
    tr <- d[t$train_idx, , drop = FALSE]; te <- d[t$test_idx, , drop = FALSE]
    .run_one(t$model, tr, te, ycol, .task_seed(ycol, 1L, t$fold_id, SEED), TRUE, predictors, NULL, 5L, ps)
  }
  t0 <- proc.time()[["elapsed"]]
  res <- .par_lapply(tasks, worker, CORES)
  oof <- rbindlist(lapply(seq_along(tasks), function(i) {
    t <- tasks[[i]]; r <- res[[i]]
    if (!is.list(r) || !isTRUE(r$ok)) stop("hold-out task failed: ", ycol, " ", t$model, " ", t$tag, " ", if (is.list(r)) r$notes else as.character(r))
    data.table(outcome = ycol, model = t$model, arm = t$tag, `repeat` = 1L, fold = t$fold_id, case_uid = d$case_uid[t$test_idx],
               y = d[[ycol]][t$test_idx], p = r$p, notes = r$notes, converged = r$converged)
  }))
  .RUNLOG[[length(.RUNLOG) + 1L]] <<- data.table(section = section, label = label, models = paste(models, collapse = "+"),
    n_tasks = length(tasks), n_failed = 0L, n_not_converged = sum(vapply(res, function(r) !isTRUE(r$converged), TRUE)),
    n_with_notes = sum(vapply(res, function(r) nzchar(r$notes), TRUE)), cpu_minutes = round(sum(vapply(res, function(r) r$seconds, 1), na.rm = TRUE) / 60, 2),
    wall_minutes = round((proc.time()[["elapsed"]] - t0) / 60, 2))
  oof
}

cache_section <- function(name, fun) {
  f <- file.path(SENS_DIR, paste0("sec_", name, ".rds"))
  if (RESUME && file.exists(f)) {
    message(">>> section ", name, " (cached)"); r <- readRDS(f)
    if (!is.null(r$runlog)) .RUNLOG[[length(.RUNLOG) + 1L]] <<- r$runlog
    return(r)
  }
  message(">>> section ", name); t0 <- proc.time()[["elapsed"]]; n0 <- length(.RUNLOG)
  r <- fun()
  r$minutes <- (proc.time()[["elapsed"]] - t0) / 60
  r$runlog <- if (length(.RUNLOG) > n0) rbindlist(.RUNLOG[(n0 + 1L):length(.RUNLOG)], fill = TRUE) else NULL
  saveRDS(r, f); message(sprintf("<<< section %s: %.1f min", name, r$minutes)); r
}

# rows for the primary OOF predictions of the same models (reference for cohort / subset comparisons)
primary_rows <- function(lab, models, analysis, subset = "primary", keep_case = NULL, ref = NULL, nb_t = NULL, seed = 7L) {
  L <- oof_lists(primary_oof(lab), "cv", models)
  if (!is.null(keep_case)) L <- subset_lists(L, L$case_uid %in% keep_case)
  sens_rows(L$y, L$P, analysis, lab, subset, ref = if (is.null(ref)) .pref(names(L$P)) else ref, nb_t = nb_t, seed = seed)
}

# ---------------------------------------------------------------- shared data
get_data <- function(cohort = "in_primary") load_modelling_data(V, cohort)
CPT_VARS <- local({ cp <- fread(file.path(P$results, "cohort_predictors.csv")); cp[cpt_derived == "yes", variable] })

# ---------------------------------------------------------------- 1. internal-external era validation + forward chaining
# era hold-out specifications: four internal-external blocks (train on the other three) + two forward-chaining splits
era_specs <- function(dd) {
  era <- as.character(dd$era); yr <- dd$puf_year
  specs <- list()
  for (i in seq_along(ERA_BLOCKS)) specs[[length(specs) + 1L]] <- list(tag = paste0("era:", ERA_BLOCKS[i]), fold_id = i,
    train_idx = which(era != ERA_BLOCKS[i]), test_idx = which(era == ERA_BLOCKS[i]))
  specs[[length(specs) + 1L]] <- list(tag = "forward:2012-2017->2018-2020", fold_id = 11L, train_idx = which(yr <= 2017), test_idx = which(yr >= 2018 & yr <= 2020))
  specs[[length(specs) + 1L]] <- list(tag = "forward:2012-2020->2022-2024", fold_id = 12L, train_idx = which(yr <= 2020), test_idx = which(yr >= 2022))
  stopifnot(all(vapply(specs, function(s) length(intersect(s$train_idx, s$test_idx)) == 0, TRUE)))
  specs
}

sec_era <- function() {
  dat <- get_data(); d_all <- dat$d; ps <- dat$ps
  full5 <- c("lr", "ridge_lasso", "xgb", "rf", "mlp")
  if (QUICK) full5 <- c("lr", "ridge_lasso", "xgb")
  cfg <- list(
    list(lab = "nonhome", models = full5, pred = ps$common, rename = NULL),
    list(lab = "los_p75", models = full5, pred = ps$common, rename = NULL),
    list(lab = "readm_unplanned", models = full5, pred = ps$common, rename = NULL),
    list(lab = "ssi", models = c("ridge_lasso", "xgb"), pred = RESTRICTED_PREDICTORS(ps), rename = NULL),
    list(lab = "cardiac", models = "ridge_lasso", pred = PARSIMONIOUS_PREDICTORS, rename = c(ridge_lasso = "pen_lr_parsimonious")))
  if (QUICK) cfg <- cfg[1:2]
  rows <- list(); het <- list(); desc <- list()
  for (cf in cfg) {
    ycol <- OUT[[cf$lab]]; dd <- outcome_data(d_all, ycol)
    specs <- era_specs(dd)
    message(sprintf("  era %s: n=%d events=%d", cf$lab, nrow(dd), sum(dd[[ycol]])))
    oof <- run_specs(dd, ycol, cf$models, specs, cf$pred, ps, "era", cf$lab)
    if (!is.null(cf$rename)) oof[model %in% names(cf$rename), model := cf$rename[model]]
    oof[, outcome := cf$lab]
    auc_draws <- list(); auc_pts <- list()
    for (s in specs) {
      L <- oof_lists(oof[arm == s$tag], s$tag)
      sr <- sens_rows(L$y, L$P, "era_validation", cf$lab, s$tag, ref = .pref(names(L$P)), B = B_ERA, seed = 1000L * s$fold_id + 17L,
                      keep_draws = startsWith(s$tag, "era:"))
      if (startsWith(s$tag, "era:")) { auc_draws[[s$tag]] <- attr(sr, "draws")["auroc", , ]; auc_pts[[s$tag]] <- sr[metric == "auroc", .(model, estimate, n, events)] }
      rows[[length(rows) + 1L]] <- sr
      desc[[length(desc) + 1L]] <- data.table(outcome = cf$lab, subset = s$tag, n_train = length(s$train_idx), n_test = length(s$test_idx),
                                              events_train = sum(dd[[ycol]][s$train_idx]), events_test = sum(L$y))
    }
    # heterogeneity across the four held-out blocks: SD (and range) of block AUROC; the interval resamples each block independently
    mods <- rownames(auc_draws[[1]])
    for (m in mods) {
      est <- vapply(auc_pts, function(z) z[model == m, estimate], 1)
      dr <- apply(vapply(auc_draws, function(a) a[m, ], numeric(B_ERA)), 1, sd)
      qq <- quantile(dr, c(0.025, 0.975), names = FALSE)
      nn <- sum(vapply(auc_pts, function(z) z[model == m, n], 1)); ee <- sum(vapply(auc_pts, function(z) z[model == m, events], 1))
      het[[length(het) + 1L]] <- data.table(analysis = "era_validation", outcome = cf$lab, model = m, subset = "four held-out era blocks",
        metric = c("auroc_sd_across_blocks", "auroc_range_across_blocks"), estimate = c(sd(est), diff(range(est))),
        lo = c(qq[1], min(est)), hi = c(qq[2], max(est)), n = as.integer(nn), events = as.integer(ee))
    }
    saveRDS(oof[, .(outcome, model, arm, fold, case_uid, y, p)], file.path(SENS_DIR, sprintf("era_oof_%s.rds", cf$lab)))
  }
  desc <- rbindlist(desc)
  # renal: descriptive only (2012-2020 legacy definition; 2021 removed by T1)
  rn <- d_all[!is.na(d_all$y_renal), ]
  renal <- rn[, c("era", "y_renal")] |> as.data.table()
  renal <- renal[, .(n = .N, events = sum(y_renal)), by = era][order(era)]
  renal_rows <- renal[, .(analysis = "era_validation", outcome = "renal", model = "descriptive", subset = as.character(era), metric = "incidence",
                          estimate = events / n, lo = NA_real_, hi = NA_real_, n = as.integer(n), events = as.integer(events))]
  list(rows = rbindlist(c(rows, het, list(renal_rows)), fill = TRUE), desc = desc)
}

# ---------------------------------------------------------------- 2. cohort sensitivity
sec_cohort <- function() {
  rows <- list()
  specs <- list(
    list(cohort = "in_dx_required", name = "dx_required"),
    list(cohort = "in_broad", name = "broad"),
    list(cohort = "in_strict", name = "strict"),
    list(cohort = "in_primary_nm", name = "primary_plus_neuromuscular"))
  if (QUICK) specs <- specs[1]
  for (sp in specs) {
    dat <- get_data(sp$cohort); d_all <- dat$d; ps <- dat$ps
    outs <- c("nonhome", "los_p75", "readm_unplanned")
    if (sp$name == "dx_required") outs <- c(outs, "ssi")
    if (QUICK) outs <- outs[1]
    for (lab in outs) {
      ycol <- OUT[[lab]]; dd <- outcome_data(d_all, ycol); ev <- sum(dd[[ycol]])
      pred <- if (lab == "ssi") RESTRICTED_PREDICTORS(ps) else ps$common
      models <- if (lab == "ssi") c("ridge_lasso", "xgb")
                else if (lab == "readm_unplanned" && ev < 500) c("ridge_lasso", "xgb") else c("lr", "ridge_lasso", "xgb")
      message(sprintf("  cohort %s / %s: n=%d events=%d models=%s", sp$name, lab, nrow(dd), ev, paste(models, collapse = ",")))
      res <- cv_run("cohort", paste("cohort", sp$name, lab), dd, ycol, models, pred, ps)
      L <- oof_lists(res$oof)
      rows[[length(rows) + 1L]] <- sens_rows(L$y, L$P, "cohort_sensitivity", lab, sp$name, ref = .pref(models), seed = 31L)
      # primary comparator: same models on the primary cohort (T2 OOF)
      pm <- intersect(models, unique(primary_oof(lab)$model))
      rows[[length(rows) + 1L]] <- primary_rows(lab, pm, "cohort_sensitivity", "primary")
    }
  }
  list(rows = rbindlist(rows))
}

# ---------------------------------------------------------------- 3. predictor-set sensitivity (a, b, c)
sec_predsets <- function() {
  dat <- get_data(); d_all <- dat$d; ps <- dat$ps
  sets <- list(
    drop_n_cpt_codes = setdiff(ps$common, "n_cpt_codes"),
    drop_all_cpt_derived = setdiff(ps$common, CPT_VARS),
    drop_labs = setdiff(ps$common, c(ps$labs, ps$missing_indicators)))
  rows <- list()
  runs <- list(list(lab = "nonhome", models = c("lr", "ridge_lasso", "xgb"), sets = names(sets)),
               list(lab = "los_p75", models = c("lr", "ridge_lasso", "xgb"), sets = names(sets)),
               list(lab = "readm_unplanned", models = c("lr", "ridge_lasso", "xgb"), sets = names(sets)),
               list(lab = "ssi", models = c("ridge_lasso", "xgb"), sets = c("drop_n_cpt_codes", "drop_all_cpt_derived")))
  if (QUICK) runs <- runs[1]
  for (r in runs) {
    ycol <- OUT[[r$lab]]; dd <- outcome_data(d_all, ycol)
    P0 <- oof_lists(primary_oof(r$lab), "cv", r$models)
    ord <- match(P0$case_uid, dd$case_uid)
    rows[[length(rows) + 1L]] <- primary_rows(r$lab, r$models, "predictor_set_sensitivity", "primary predictors")
    for (sn in r$sets) {
      pred <- sets[[sn]]
      if (r$lab == "ssi") pred <- setdiff(RESTRICTED_PREDICTORS(ps), if (sn == "drop_n_cpt_codes") "n_cpt_codes" else CPT_VARS)
      res <- cv_run("predsets", paste(sn, r$lab), dd, ycol, r$models, pred, ps)
      L <- oof_lists(res$oof); stopifnot(identical(L$case_uid, P0$case_uid))
      rows[[length(rows) + 1L]] <- sens_rows(L$y, L$P, "predictor_set_sensitivity", r$lab, sn, ref = .pref(r$models), base = P0$P, seed = 41L)
    }
  }
  list(rows = rbindlist(rows))
}

# ---------------------------------------------------------------- 3d. multiple imputation vs median + indicator
sec_mi <- function() {
  suppressPackageStartupMessages(library(mice))
  dat <- get_data(); d_all <- dat$d; ps <- dat$ps
  vars <- setdiff(ps$common, ps$missing_indicators)          # indicators are not needed after imputation
  X <- d_all[, vars, drop = FALSE]
  meth <- vapply(vars, function(v) if (is.numeric(X[[v]]) && anyNA(X[[v]])) "pmm" else "", "")
  M <- if (QUICK) 2L else 5L
  message(sprintf("  MI: m=%d, variables with missing values: %s", M, paste(vars[meth == "pmm"], collapse = ",")))
  imps <- parallel::mclapply(seq_len(M), function(i) {
    set.seed(SEED + 100L + i)
    imp <- suppressWarnings(mice(X, m = 1, maxit = 10, method = meth, seed = SEED + 100L + i, printFlag = FALSE))   # outcome is NOT in X
    mice::complete(imp, 1)
  }, mc.cores = min(M, CORES), mc.set.seed = FALSE)
  stopifnot(all(vapply(imps, function(z) !anyNA(z), TRUE)))
  rows <- list()
  runs <- list(list(lab = "nonhome", models = c("lr", "ridge_lasso", "xgb")), list(lab = "los_p75", models = c("lr", "ridge_lasso", "xgb")),
               list(lab = "readm_unplanned", models = c("lr", "ridge_lasso", "xgb")))
  if (QUICK) runs <- runs[1]
  for (r in runs) {
    ycol <- OUT[[r$lab]]; ok <- !is.na(d_all[[ycol]])
    P0 <- oof_lists(primary_oof(r$lab), "cv", r$models)
    rows[[length(rows) + 1L]] <- primary_rows(r$lab, r$models, "multiple_imputation", "median + missingness indicators (primary)")
    # (i) median imputation without the indicators (isolates the effect of the indicators)
    dd <- outcome_data(d_all, ycol)
    res0 <- cv_run("mi", paste("median_noind", r$lab), dd, ycol, r$models, vars, ps)
    L0 <- oof_lists(res0$oof)
    rows[[length(rows) + 1L]] <- sens_rows(L0$y, L0$P, "multiple_imputation", r$lab, "median, no indicators", ref = .pref(r$models), base = P0$P, seed = 51L)
    # (ii) MI: each imputed data set through the same nested CV; predictions averaged over the m imputations per patient
    allp <- list()
    for (i in seq_len(M)) {
      di <- d_all; di[vars] <- imps[[i]]; di <- outcome_data(di, ycol)
      res <- cv_run("mi", sprintf("mi%d %s", i, r$lab), di, ycol, r$models, vars, ps)
      allp[[i]] <- res$oof[, .(model, `repeat`, fold, case_uid, y, p)]
    }
    avg <- rbindlist(allp)[, .(p = mean(p), y = y[1], n_imp = .N), by = .(model, `repeat`, fold, case_uid)]
    stopifnot(all(avg$n_imp == M))
    avg[, `:=`(arm = "cv", outcome = r$lab)]
    L <- oof_lists(avg)
    rows[[length(rows) + 1L]] <- sens_rows(L$y, L$P, "multiple_imputation", r$lab, sprintf("multiple imputation m=%d (pmm, outcome excluded), predictions averaged", M),
                                           ref = .pref(r$models), base = P0$P, seed = 53L)
  }
  list(rows = rbindlist(rows))
}

# ---------------------------------------------------------------- db pull (operative time; in-hospital events for the discharge-time model)
DB_PULL_FILE <- function() file.path(SENS_DIR, "db_pull.rds")
pull_db <- function(ids) {
  if (file.exists(DB_PULL_FILE()) && RESUME) return(readRDS(DB_PULL_FILE()))
  comp <- c("SUPINFEC", "WNDINFD", "ORGSPCSSI", "DEHIS", "OUPNEUMO", "REINTUB", "PULEMBOL", "FAILWEAN", "RENAINSF", "OPRENAFL",
            "URNINFEC", "CNSCVA", "CDARREST", "CDMI", "OTHDVT", "OTHSYSEP", "OTHSESHOCK")
  dcols <- paste0("D", comp)
  q <- sprintf("SELECT case_uid, OPTIME, DISCHDEST, DOPTODIS, %s FROM main WHERE puf_year BETWEEN 2012 AND 2024 AND case_uid IN ('%s')",
               paste(c(comp, dcols), collapse = ", "), paste(ids, collapse = "','"))
  x <- nsqip_query(q)
  stopifnot(nrow(x) == length(ids), !anyDuplicated(x$case_uid))
  for (cc in names(x)[vapply(x, is.character, TRUE)]) set(x, j = cc, value = trimws(x[[cc]]))
  inhosp <- vapply(comp, function(cn) { v <- x[[cn]]; dy <- x[[paste0("D", cn)]]
    !is.na(v) & v != "No Complication" & v != "NULL" & !is.na(dy) & !is.na(x$DOPTODIS) & dy <= x$DOPTODIS }, logical(nrow(x)))
  HOME <- c("Home", "Facility Which was Home", "Home/Permanent Residence")
  NONHOME <- c("Skilled Care, Not Home", "Unskilled Facility Not Home", "Rehab", "Separate Acute Care", "Multi-level Senior Community",
               "Hospice", "Other Facility", "Acute Care Hospital")
  out <- data.table(case_uid = x$case_uid, optime = ifelse(!is.na(x$OPTIME) & x$OPTIME > 0 & x$OPTIME <= 1200, x$OPTIME, NA_real_),
                    dt_dest = factor(ifelse(x$DISCHDEST %in% HOME, "Home", ifelse(x$DISCHDEST %in% NONHOME, "Non-home", "Other or unknown")),
                                     levels = c("Home", "Non-home", "Other or unknown")),
                    dt_n_inhospital_events = rowSums(inhosp))
  out[, dt_any_inhospital_event := as.integer(dt_n_inhospital_events > 0)]
  saveRDS(out, DB_PULL_FILE()); out
}

# ---------------------------------------------------------------- 3e. operative time add-on (day-of-surgery; NEVER in primary models)
sec_optime <- function() {
  dat <- get_data(); d_all <- dat$d; ps <- dat$ps
  pull <- pull_db(d_all$case_uid)
  d_all$optime <- pull$optime[match(d_all$case_uid, pull$case_uid)]
  ps_ot <- ps; ps_ot$common <- c(ps$common, "optime"); ps_ot$numeric <- c(ps$numeric, "optime")   # local copy: primary ps untouched
  rows <- list()
  runs <- list(list(lab = "nonhome"), list(lab = "los_p75"), list(lab = "readm_unplanned"))
  if (QUICK) runs <- runs[1]
  for (r in runs) {
    ycol <- OUT[[r$lab]]; dd <- outcome_data(d_all, ycol); models <- c("lr", "ridge_lasso", "xgb")
    P0 <- oof_lists(primary_oof(r$lab), "cv", models)
    rows[[length(rows) + 1L]] <- primary_rows(r$lab, models, "day_of_surgery_addon_optime", "preoperative predictors only (primary)")
    res <- cv_run("optime", paste("optime", r$lab), dd, ycol, models, c(ps$common, "optime"), ps_ot)
    L <- oof_lists(res$oof); stopifnot(identical(L$case_uid, P0$case_uid))
    rows[[length(rows) + 1L]] <- sens_rows(L$y, L$P, "day_of_surgery_addon_optime", r$lab, "preoperative predictors + operative time (day-of-surgery add-on, NOT a primary model)",
                                           ref = "lr", base = P0$P, seed = 61L)
  }
  list(rows = rbindlist(rows), optime_missing = mean(is.na(d_all$optime)))
}

# ---------------------------------------------------------------- 4. fold-internal logistic recalibration of flexible models
sec_recal <- function() {
  dat <- get_data(); d_all <- dat$d; ps <- dat$ps
  runs <- list(list(lab = "nonhome", models = c("rf", "xgb"), pred = ps$common),
               list(lab = "los_p75", models = c("rf", "xgb"), pred = ps$common),
               list(lab = "readm_unplanned", models = c("rf", "xgb"), pred = ps$common),
               list(lab = "ssi", models = "xgb", pred = RESTRICTED_PREDICTORS(ps)))
  if (QUICK) runs <- list(list(lab = "nonhome", models = "xgb", pred = ps$common))
  rows <- list(); calpar <- list()
  for (r in runs) {
    ycol <- OUT[[r$lab]]; dd <- as.data.frame(outcome_data(d_all, ycol)); y <- dd[[ycol]]
    O <- primary_oof(r$lab, full = TRUE)[arm == "cv" & `repeat` == 1L]
    fold_of <- O[model == r$models[1]]$fold[match(dd$case_uid, O[model == r$models[1]]$case_uid)]
    stopifnot(!anyNA(fold_of))
    tasks <- list()
    for (m in .heavy_first(r$models)) for (f in 1:5) {
      tr_idx <- which(fold_of != f)
      cf <- make_folds(y[tr_idx], 5, .task_seed(ycol, 1L, 500L + f, SEED))
      for (j in 1:5) tasks[[length(tasks) + 1L]] <- list(model = m, f = f, j = j, tr_idx = tr_idx, hold = tr_idx[cf == j], fit = tr_idx[cf != j])
    }
    worker <- function(t) {
      pr <- .run_one(t$model, dd[t$fit, , drop = FALSE], dd[t$hold, , drop = FALSE], ycol, .task_seed(ycol, 1L, 1000L + 10L * t$f + t$j, SEED),
                     TRUE, r$pred, NULL, 5L, ps)
      list(ok = pr$ok, p = pr$p, hold = t$hold, model = t$model, f = t$f, secs = pr$seconds)
    }
    t0 <- proc.time()[["elapsed"]]
    res <- .par_lapply(tasks, worker, CORES)
    if (any(!vapply(res, function(z) is.list(z) && isTRUE(z$ok), TRUE))) stop("recalibration inner task failed: ", r$lab)
    .RUNLOG[[length(.RUNLOG) + 1L]] <<- data.table(section = "recal", label = r$lab, models = paste(r$models, collapse = "+"), n_tasks = length(tasks), n_failed = 0L,
      n_not_converged = 0L, n_with_notes = 0L, cpu_minutes = round(sum(vapply(res, function(z) z$secs, 1), na.rm = TRUE) / 60, 2),
      wall_minutes = round((proc.time()[["elapsed"]] - t0) / 60, 2))
    message(sprintf("    [recal] %s done (%.1f min wall)", r$lab, (proc.time()[["elapsed"]] - t0) / 60))
    recal_rows <- list()
    for (m in r$models) {
      pm <- O[model == m]; setorder(pm, case_uid)
      pcal <- rep(NA_real_, nrow(dd))
      for (f in 1:5) {
        pin <- rep(NA_real_, nrow(dd))
        for (z in Filter(function(z) z$model == m && z$f == f, res)) pin[z$hold] <- z$p
        tr_idx <- which(fold_of != f); stopifnot(!anyNA(pin[tr_idx]))
        lp <- qlogis(pmin(pmax(pin[tr_idx], 1e-6), 1 - 1e-6))
        fit <- suppressWarnings(glm.fit(cbind(1, lp), y[tr_idx], family = binomial()))
        ab <- unname(fit$coefficients); calpar[[length(calpar) + 1L]] <- data.table(outcome = r$lab, model = m, fold = f, intercept = ab[1], slope = ab[2])
        te <- which(fold_of == f)
        pte <- pm$p[match(dd$case_uid[te], pm$case_uid)]
        pcal[te] <- plogis(ab[1] + ab[2] * qlogis(pmin(pmax(pte, 1e-6), 1 - 1e-6)))
      }
      stopifnot(!anyNA(pcal))
      recal_rows[[m]] <- data.table(case_uid = dd$case_uid, y = y, p_uncal = pm$p[match(dd$case_uid, pm$case_uid)], p_cal = pcal)
    }
    ord <- order(dd$case_uid)
    ref <- .pref(unique(O$model))
    refm <- O[model == ref]; setorder(refm, case_uid)
    ids <- sort(dd$case_uid)
    Pl <- list(); Pl[[ref]] <- matrix(refm$p[match(ids, refm$case_uid)], ncol = 1)
    for (m in r$models) {
      rr <- recal_rows[[m]]
      Pl[[m]] <- matrix(rr$p_uncal[match(ids, rr$case_uid)], ncol = 1)
      Pl[[paste0(m, "_recal")]] <- matrix(rr$p_cal[match(ids, rr$case_uid)], ncol = 1)
    }
    yy <- y[match(ids, dd$case_uid)]
    rows[[length(rows) + 1L]] <- sens_rows(yy, Pl, "fold_internal_recalibration", r$lab, "outer-fold test predictions, repeat 1",
                                           ref = ref, nb_t = THR[[r$lab]], seed = 71L)
  }
  list(rows = rbindlist(rows), calpar = rbindlist(calpar))
}

# ---------------------------------------------------------------- 5. periods, posterior-only, LOS thresholds, discharge-time readmission
sec_periods <- function() {
  dat <- get_data(); d_all <- dat$d; ps <- dat$ps
  rows <- list(); notes <- list()
  base3 <- c("lr", "ridge_lasso", "xgb")
  outs <- c("nonhome", "los_p75", "readm_unplanned")
  if (QUICK) outs <- outs[1]
  # (a) pre-2020 vs 2020-2024: the same modelling pipeline run inside each period (nested CV), plus primary OOF evaluated in the period
  periods <- list(pre2020 = which(d_all$puf_year <= 2019), `2020_2024` = which(d_all$puf_year >= 2020))
  for (pn in names(periods)) {
    dp <- d_all[periods[[pn]], , drop = FALSE]
    for (lab in outs) {
      ycol <- OUT[[lab]]; dd <- outcome_data(dp, ycol)
      res <- cv_run("periods", paste("period", pn, lab), dd, ycol, base3, ps$common, ps)
      L <- oof_lists(res$oof)
      rows[[length(rows) + 1L]] <- sens_rows(L$y, L$P, "period_comparison", lab, paste0(pn, " (models refit within period)"), ref = "lr", seed = 81L)
      rows[[length(rows) + 1L]] <- primary_rows(lab, base3, "period_comparison", paste0(pn, " (primary OOF restricted to period)"), keep_case = dd$case_uid, ref = "lr")
    }
  }
  # (b) posterior-only (no anterior approach or anterior deformity arthrodesis CPT)
  post <- which(d_all$combined_ap == 0 & as.character(d_all$ant_deformity_arthrodesis) == "None")
  dpost <- d_all[post, , drop = FALSE]
  notes$posterior_n <- length(post)
  for (lab in outs) {
    ycol <- OUT[[lab]]; dd <- outcome_data(dpost, ycol)
    res <- cv_run("periods", paste("posterior_only", lab), dd, ycol, base3, ps$common, ps)
    L <- oof_lists(res$oof)
    rows[[length(rows) + 1L]] <- sens_rows(L$y, L$P, "posterior_only", lab, "posterior-only subset (models refit in subset)", ref = "lr", seed = 83L)
    rows[[length(rows) + 1L]] <- primary_rows(lab, base3, "posterior_only", "primary OOF restricted to posterior-only", keep_case = dd$case_uid, ref = "lr")
  }
  # (c) alternative LOS thresholds: >= 9 days and the cohort 90th percentile (same eligibility as the primary LOS outcome)
  elig <- !is.na(d_all$y_los_p75)
  q90 <- as.numeric(quantile(d_all$los_days[elig], 0.9, type = 1))
  d_all$y_los_p90 <- ifelse(elig, as.integer(d_all$los_days >= q90), NA_integer_)
  notes$los_q90 <- q90
  for (lab_y in c("y_los_ge9", "y_los_p90")) {
    dd <- outcome_data(d_all, lab_y); lab <- if (lab_y == "y_los_ge9") "los_ge9" else "los_p90"
    res <- cv_run("periods", paste("los_alt", lab), dd, lab_y, base3, ps$common, ps)
    L <- oof_lists(res$oof)
    rows[[length(rows) + 1L]] <- sens_rows(L$y, L$P, "alternative_los_threshold", lab, if (lab == "los_ge9") "LOS >= 9 days" else sprintf("LOS >= %d days (90th percentile)", q90),
                                           ref = "lr", seed = 85L)
  }
  # (d) readmission discharge-time model: SEPARATE secondary analysis; models are named dt_* and never share names with primary models
  pull <- pull_db(d_all$case_uid)
  d_dt <- d_all
  d_dt$dt_los_days <- d_dt$los_days
  d_dt$dt_dest <- pull$dt_dest[match(d_dt$case_uid, pull$case_uid)]
  d_dt$dt_n_inhospital_events <- pull$dt_n_inhospital_events[match(d_dt$case_uid, pull$case_uid)]
  d_dt$dt_any_inhospital_event <- pull$dt_any_inhospital_event[match(d_dt$case_uid, pull$case_uid)]
  add <- c("dt_los_days", "dt_dest", "dt_n_inhospital_events", "dt_any_inhospital_event")
  ps_dt <- ps; ps_dt$common <- c(ps$common, add); ps_dt$numeric <- c(ps$numeric, "dt_los_days", "dt_n_inhospital_events")
  ps_dt$categorical <- c(ps$categorical, "dt_dest"); ps_dt$binary <- c(ps$binary, "dt_any_inhospital_event")
  dd <- outcome_data(d_dt, "y_readm_unplanned")
  P0 <- oof_lists(primary_oof("readm_unplanned"), "cv", base3)
  res <- cv_run("periods", "discharge_time readm", dd, "y_readm_unplanned", base3, c(ps$common, add), ps_dt)
  res$oof[, model := paste0("dt_", model)]
  L <- oof_lists(res$oof)
  stopifnot(identical(L$case_uid, P0$case_uid))
  base <- setNames(P0$P, paste0("dt_", names(P0$P)))
  rows[[length(rows) + 1L]] <- sens_rows(L$y, L$P, "readmission_discharge_time_model", "readm_unplanned_dischargetime",
                                         "secondary: adds LOS, discharge destination, in-hospital events (NOT the primary preoperative model)", ref = "dt_lr", base = base, seed = 87L)
  rows[[length(rows) + 1L]] <- primary_rows("readm_unplanned", base3, "readmission_discharge_time_model", "primary preoperative-only models (comparator)", ref = "lr")
  notes$dt_inhospital_event_prevalence <- mean(pull$dt_any_inhospital_event[match(dd$case_uid, pull$case_uid)])
  list(rows = rbindlist(rows), notes = notes)
}

# ---------------------------------------------------------------- 6. subgroups (from primary OOF; read only)
subgroup_defs <- function(d) {
  ag <- cut(d$age, c(-Inf, 64.999, 74.999, Inf), labels = c("<65", "65-74", ">=75"))
  list(age_group = as.character(ag), sex = ifelse(d$sex_male == 1, "Male", "Female"),
       race = as.character(d$race_group), ethnicity = paste("Hispanic:", as.character(d$hispanic)))
}
sec_subgroups <- function() {
  dat <- get_data(); d_all <- dat$d
  rows <- list()
  runs <- list(list(lab = "nonhome", models = c("lr", "ridge_lasso", "xgb", "rf", "mlp")), list(lab = "los_p75", models = c("lr", "ridge_lasso", "xgb", "rf", "mlp")),
               list(lab = "readm_unplanned", models = c("lr", "ridge_lasso", "xgb", "rf", "mlp")), list(lab = "ssi", models = c("ridge_lasso", "xgb")),
               list(lab = "cardiac", models = "pen_lr_parsimonious"))
  if (QUICK) runs <- runs[1]
  for (r in runs) {
    L <- oof_lists(primary_oof(r$lab, full = TRUE), "cv", r$models)
    dd <- d_all[match(L$case_uid, d_all$case_uid), , drop = FALSE]
    sg <- subgroup_defs(dd)
    rows[[length(rows) + 1L]] <- sens_rows(L$y, L$P, "subgroup", r$lab, "overall", ref = .pref(r$models), nb_t = THR[[r$lab]], seed = 91L)
    for (v in names(sg)) for (lv in sort(unique(sg[[v]]))) {
      keep <- which(sg[[v]] == lv); ev <- sum(L$y[keep])
      sub <- paste0(v, "=", lv)
      if (!subgroup_estimable(ev)) {
        rows[[length(rows) + 1L]] <- data.table(analysis = "subgroup", outcome = r$lab, model = r$models, subset = sub, metric = "not estimable (events < 100)",
                                                estimate = NA_real_, lo = NA_real_, hi = NA_real_, n = length(keep), events = as.integer(ev))
        next
      }
      Ls <- subset_lists(L, keep)
      rows[[length(rows) + 1L]] <- sens_rows(Ls$y, Ls$P, "subgroup", r$lab, sub, ref = .pref(r$models), nb_t = THR[[r$lab]], seed = 93L)
    }
  }
  list(rows = rbindlist(rows, fill = TRUE))
}

# ---------------------------------------------------------------- 7. bootstrap instability (Riley and Collins)
sec_instability <- function() {
  dat <- get_data(); d_all <- dat$d; ps <- dat$ps
  models <- c("lr", "ridge_lasso", "xgb")
  outs <- c("nonhome", "los_p75", "readm_unplanned"); if (QUICK) outs <- outs[1]
  rows <- list(); pts <- list(); dec <- list()
  for (lab in outs) {
    ycol <- OUT[[lab]]; dd <- as.data.frame(outcome_data(d_all, ycol)); n <- nrow(dd)
    ev_idx <- with_seed(SEED + 5L + sum(utf8ToInt(lab)), sort(sample.int(n, min(n, if (QUICK) 100L else 1000L))))   # fixed evaluation subset
    # models: lr, ridge_lasso, xgb (full tuning repeated in every bootstrap fit) and xgb_fixed_config (xgb with the T2 default configuration,
    # tune = FALSE). The duplicated rows of a bootstrap sample leak across inner CV folds, which favours deeper xgb configurations and
    # inflates the instability of any tuned flexible model; the fixed-configuration arm shows instability without that artefact.
    inst_models <- c(models, "xgb_fixed_config")
    tasks <- list()
    for (m in inst_models) {
      tasks[[length(tasks) + 1L]] <- list(model = m, b = 0L)
      for (b in seq_len(B_INST)) tasks[[length(tasks) + 1L]] <- list(model = m, b = b)
    }
    tasks <- tasks[order(match(vapply(tasks, `[[`, "", "model"), c("xgb", "xgb_fixed_config", "ridge_lasso", "lr")))]
    worker <- function(t) {
      tr <- if (t$b == 0L) dd else dd[with_seed(SEED + 9001L * t$b + sum(utf8ToInt(lab)), sample.int(n, n, replace = TRUE)), , drop = FALSE]
      fixed <- t$model == "xgb_fixed_config"
      r <- .run_one(if (fixed) "xgb" else t$model, tr, dd[ev_idx, , drop = FALSE], ycol, .task_seed(ycol, 1L, 7000L + t$b, SEED), !fixed, ps$common, NULL, 5L, ps)
      list(ok = r$ok, p = r$p, model = t$model, b = t$b, secs = r$seconds)
    }
    t0 <- proc.time()[["elapsed"]]
    res <- .par_lapply(tasks, worker, CORES)
    if (any(!vapply(res, function(z) is.list(z) && isTRUE(z$ok), TRUE))) stop("instability task failed: ", lab)
    .RUNLOG[[length(.RUNLOG) + 1L]] <<- data.table(section = "instability", label = lab, models = paste(inst_models, collapse = "+"), n_tasks = length(tasks), n_failed = 0L,
      n_not_converged = 0L, n_with_notes = 0L, cpu_minutes = round(sum(vapply(res, function(z) z$secs, 1), na.rm = TRUE) / 60, 2),
      wall_minutes = round((proc.time()[["elapsed"]] - t0) / 60, 2))
    message(sprintf("    [instability] %s done (%.1f min wall)", lab, (proc.time()[["elapsed"]] - t0) / 60))
    for (m in inst_models) {
      rm_ <- Filter(function(z) z$model == m, res)
      p0 <- rm_[[which(vapply(rm_, function(z) z$b == 0L, TRUE))]]$p
      Pb <- do.call(cbind, lapply(rm_[vapply(rm_, function(z) z$b > 0L, TRUE)], function(z) z$p))
      mape_i <- rowMeans(abs(Pb - p0))
      qq <- quantile(mape_i, c(0.025, 0.975), names = FALSE)
      yv <- dd[[ycol]][ev_idx]
      # calibration instability: slope of each bootstrap model's predictions on the original-model logit
      rows[[length(rows) + 1L]] <- data.table(analysis = "bootstrap_instability", outcome = lab, model = m,
        subset = sprintf("%d bootstrap refits, fixed evaluation subset of %d patients", B_INST, length(ev_idx)), metric = "mape_mean_abs_prediction_error",
        estimate = mean(mape_i), lo = qq[1], hi = qq[2], n = as.integer(n), events = as.integer(sum(dd[[ycol]])))
      rows[[length(rows) + 1L]] <- data.table(analysis = "bootstrap_instability", outcome = lab, model = m,
        subset = sprintf("%d bootstrap refits, fixed evaluation subset of %d patients", B_INST, length(ev_idx)), metric = "mape_relative_to_mean_prediction",
        estimate = mean(mape_i) / mean(p0), lo = qq[1] / mean(p0), hi = qq[2] / mean(p0), n = as.integer(n), events = as.integer(sum(dd[[ycol]])))
      rng <- apply(Pb, 1, function(v) quantile(v, c(0.025, 0.5, 0.975), names = FALSE))
      pts[[length(pts) + 1L]] <- data.table(outcome = lab, model = m, point = seq_along(ev_idx), p_original = p0, p_boot_lo = rng[1, ], p_boot_median = rng[2, ],
                                            p_boot_hi = rng[3, ], mape = mape_i, observed = yv)
      brk <- unique(quantile(p0, seq(0, 1, 0.1), names = FALSE))
      bin <- cut(p0, brk, include.lowest = TRUE)
      dec[[length(dec) + 1L]] <- data.table(outcome = lab, model = m, bin = bin, p0 = p0, mape = mape_i, w = rng[3, ] - rng[1, ])[
        , .(mean_original_prediction = mean(p0), mean_mape = mean(mape), mean_boot_interval_width = mean(w), n = .N), by = bin]
    }
  }
  pts <- rbindlist(pts); dec <- rbindlist(dec)
  saveRDS(pts, file.path(SENS_DIR, "instability_points.rds"))
  list(rows = rbindlist(rows), points = pts, bins = dec)
}

# ---------------------------------------------------------------- driver
run_all_sections <- function() {
  if (!RESUME && dir.exists(TASK_DIR)) unlink(TASK_DIR, recursive = TRUE)
  dir.create(TASK_DIR, showWarnings = FALSE, recursive = TRUE)
  t_start <- proc.time()[["elapsed"]]
  fns <- list(era = sec_era, cohort = sec_cohort, predsets = sec_predsets, mi = sec_mi, optime = sec_optime, recal = sec_recal,
              periods = sec_periods, subgroups = sec_subgroups, instability = sec_instability)
  secs <- list()
  for (s in names(fns)) if (s %in% SECTIONS) {
    secs[[s]] <- cache_section(s, fns[[s]])
    # persist the run log after each section so an interrupted run keeps its timings
    saveRDS(rbindlist(.RUNLOG, fill = TRUE), file.path(SENS_DIR, "runlog_partial.rds"))
  }
  total <- sum(vapply(secs, function(z) z$minutes, 1))   # compute minutes of the sections (cached sections keep their original timing)
  list(secs = secs, total_minutes = total)
}

assemble_outputs <- function(secs, total_minutes = NA_real_) {
  rows <- rbindlist(lapply(secs, function(s) s$rows), fill = TRUE)
  rows <- rows[, ..TIDY_COLS]
  rows[, c("estimate", "lo", "hi") := lapply(.SD, function(x) round(x, 5)), .SDcols = c("estimate", "lo", "hi")]
  fwrite(rows, file.path(RES_DIR, "sens_all.csv"))
  fwrite(rows[analysis == "era_validation"], file.path(RES_DIR, "sens_era.csv"))
  sg <- rows[analysis == "subgroup"]
  sg[, `:=`(subgroup_var = sub("=.*$", "", subset), subgroup = sub("^[^=]*=", "", subset))]
  sg[subset == "overall", `:=`(subgroup_var = "overall", subgroup = "all")]
  sg[, estimable := !grepl("not estimable", metric)]
  sg[, note := fifelse(estimable, "", "fewer than 100 events in subgroup: no estimate reported")]
  fwrite(sg[, .(analysis, outcome, model, subgroup_var, subgroup, metric, estimate, lo, hi, n, events, estimable, note)], file.path(RES_DIR, "sens_subgroups.csv"))
  if (!is.null(secs$instability)) fwrite(secs$instability$points[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 5) else x)], file.path(SENS_DIR, "sens_instability.csv"))
  if (!is.null(secs$instability)) fwrite(secs$instability$bins[, mean_original_prediction := round(mean_original_prediction, 5)][, lapply(.SD, function(x) if (is.numeric(x)) round(x, 5) else x)],
                                        file.path(RES_DIR, "sens_instability_bins.csv"))
  rl <- rbindlist(.RUNLOG, fill = TRUE)
  if (nrow(rl)) { rl[, total_minutes_this_run := round(total_minutes, 1)]; fwrite(rl, file.path(RES_DIR, "sens_runlog.csv")) }
  invisible(rows)
}

# ---------------------------------------------------------------- headline summary (counts and plain-language reading)
.fmt <- function(x, d = 3) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = d))
.ci <- function(e, l, h, d = 3) sprintf("%s (%s to %s)", .fmt(e, d), .fmt(l, d), .fmt(h, d))
.class_delta <- function(est, lo, hi, margin = 0.02) {
  fifelse(is.na(est), "NA", fifelse(lo > margin, "ML better beyond 0.02 margin", fifelse(lo > 0, "ML better, within margin",
    fifelse(hi < 0, "LR better", "no clear difference"))))
}

write_headline <- function(rows, out) {
  s <- as.data.table(rows)
  L <- character()
  add <- function(...) L <<- c(L, paste0(...))
  add("# Sensitivity analyses (T4): headline summary")
  add("")
  add(sprintf("Generated by `analysis/04_sensitivity.R`; every number below is in `results/sens_all.csv` (tidy rows: analysis, outcome, model, subset, metric, estimate, lo, hi, n, events). Intervals are percentile bootstrap (200 draws; 500 for era blocks; patients resampled, paired across models). Instability uses %d bootstrap refits.", B_INST))
  add("")
  add(sprintf("Total run time of this script: %.1f minutes; nested-CV refits use %d repeats x 5 folds (inner 5-fold tuning) with the same fit_predict as T2.", out$total_minutes, REPEATS))
  add("")
  # ---- counts
  add("## Counts")
  add("")
  dd <- get_data()$d
  add("| Outcome | n | events |"); add("|---|---|---|")
  for (lab in c("nonhome", "los_p75", "readm_unplanned", "ssi", "cardiac")) {
    y <- dd[[OUT[[lab]]]]; add(sprintf("| %s | %d | %d |", lab, sum(!is.na(y)), sum(y, na.rm = TRUE)))
  }
  add("")
  rl <- rbindlist(.RUNLOG, fill = TRUE)
  if (nrow(rl)) add(sprintf("Model-fitting tasks: %d, failed: %d, not converged: %d.", sum(rl$n_tasks), sum(rl$n_failed), sum(rl$n_not_converged)))
  add("")
  # ---- ML minus LR classification across settings
  add("## Does the ML versus logistic regression conclusion change?")
  add("")
  add("Classification of the paired AUROC difference (model minus logistic regression; for Tier 2 readmission and SSI the reference is penalised LR): ML better beyond the 0.02 margin, ML better within margin, LR better, or no clear difference (95% interval spans 0).")
  add("")
  d <- s[metric %like% "^delta_auroc_vs_(lr|ridge_lasso|dt_lr)$" & model %in% c("xgb", "rf", "mlp", "ridge_lasso", "dt_xgb", "dt_ridge_lasso") & analysis != "subgroup" & analysis != "era_validation"]
  d[, cls := .class_delta(estimate, lo, hi)]
  prim <- d[subset %in% c("primary", "primary predictors", "median + missingness indicators (primary)") & analysis %in% c("cohort_sensitivity", "predictor_set_sensitivity")]
  prim <- unique(prim[, .(outcome, model, prim_est = estimate, prim_cls = cls)], by = c("outcome", "model"))
  d2 <- merge(d, prim, by = c("outcome", "model"), all.x = TRUE)
  d2 <- d2[!(subset %like% "primary" & analysis %in% c("cohort_sensitivity", "predictor_set_sensitivity"))]
  d2 <- d2[order(match(outcome, c("nonhome", "los_p75", "readm_unplanned")), analysis, subset, model)]
  add("| Analysis | Setting | Outcome | Model | AUROC difference vs reference (95% CI) | Class | Primary class | Same as primary |")
  add("|---|---|---|---|---|---|---|---|")
  for (i in seq_len(nrow(d2))) {
    r <- d2[i]
    add(sprintf("| %s | %s | %s | %s | %s | %s | %s | %s |", r$analysis, r$subset, r$outcome, r$model, .ci(r$estimate, r$lo, r$hi), r$cls,
                ifelse(is.na(r$prim_cls), "n/a", r$prim_cls), ifelse(is.na(r$prim_cls), "n/a", ifelse(r$cls == r$prim_cls, "yes", "NO"))))
  }
  add("")
  chk <- d2[!is.na(prim_cls)]
  if (nrow(chk)) {
    add(sprintf("Of %d model-by-setting comparisons that have a primary counterpart, %d (%.0f%%) fall in the same class as the primary analysis; %d change class.",
                nrow(chk), sum(chk$cls == chk$prim_cls), 100 * mean(chk$cls == chk$prim_cls), sum(chk$cls != chk$prim_cls)))
    chg <- chk[cls != prim_cls]
    if (nrow(chg)) add("Class changes: ", paste(sprintf("%s / %s / %s / %s (%s -> %s)", chg$analysis, chg$subset, chg$outcome, chg$model, chg$prim_cls, chg$cls), collapse = "; "), ".")
    add("")
  }
  # ---- era
  add("## Internal-external era validation")
  add("")
  e <- s[analysis == "era_validation" & metric == "auroc_sd_across_blocks"]
  add("| Outcome | Model | SD of block AUROCs (95% CI) |"); add("|---|---|---|")
  for (i in seq_len(nrow(e))) add(sprintf("| %s | %s | %s |", e$outcome[i], e$model[i], .ci(e$estimate[i], e$lo[i], e$hi[i])))
  add("")
  f <- s[analysis == "era_validation" & metric == "auroc" & grepl("^forward", subset)]
  add("Forward chaining AUROC:")
  add(""); add("| Outcome | Model | Split | AUROC (95% CI) | test n | test events |"); add("|---|---|---|---|---|---|")
  for (i in seq_len(nrow(f))) add(sprintf("| %s | %s | %s | %s | %d | %d |", f$outcome[i], f$model[i], f$subset[i], .ci(f$estimate[i], f$lo[i], f$hi[i]), f$n[i], f$events[i]))
  add("")
  add("Cardiac and SSI era blocks have fewer than 100 events per block (cardiac 23-32, SSI 56-80): AUROC is shown but calibration slope and intercept are withheld. Renal is descriptive only (incidence by era, 2012-2020 definition).")
  add("")
  # ---- subgroups
  add("## Subgroups")
  add("")
  sg <- s[analysis == "subgroup"]
  ne <- sg[grepl("not estimable", metric), .(outcome, subset, events)] |> unique()
  est <- sg[metric == "auroc" & subset != "overall", .(outcome, subset, events)] |> unique()
  add(sprintf("Estimable subgroup-by-outcome cells (>= 100 events): %d; not estimable: %d.", nrow(est), nrow(ne)))
  if (nrow(ne)) add("Not estimable: ", paste(sprintf("%s %s (%d events)", ne$outcome, ne$subset, ne$events), collapse = "; "), ".")
  add("")
  # ---- recalibration
  rc <- s[analysis == "fold_internal_recalibration" & metric %in% c("brier", "cal_slope") & grepl("^(rf|xgb|mlp)", model)]
  if (nrow(rc)) {
    add("## Fold-internal logistic recalibration (flexible models)")
    add(""); add("| Outcome | Model | Metric | Estimate (95% CI) |"); add("|---|---|---|---|")
    for (i in seq_len(nrow(rc))) add(sprintf("| %s | %s | %s | %s |", rc$outcome[i], rc$model[i], rc$metric[i], .ci(rc$estimate[i], rc$lo[i], rc$hi[i])))
    add("")
  }
  # ---- instability
  ins <- s[analysis == "bootstrap_instability" & metric == "mape_mean_abs_prediction_error"]
  if (nrow(ins)) {
    add("## Bootstrap instability")
    add(""); add("Mean absolute prediction error across bootstrap refits (interval = 2.5th to 97.5th percentile of per-patient values):"); add("")
    add("| Outcome | Model | MAPE (spread) |"); add("|---|---|---|")
    for (i in seq_len(nrow(ins))) add(sprintf("| %s | %s | %s |", ins$outcome[i], ins$model[i], .ci(ins$estimate[i], ins$lo[i], ins$hi[i], 4)))
    add("")
  }
  add("## Reading")
  add("")
  add("Day-of-surgery operative time and the discharge-time readmission model are labelled secondary and are never part of the primary models. The multiple-imputation analysis imputes predictors only (outcome excluded) on the whole analytic cohort before nested CV, so test-fold covariates (not outcomes) inform the imputation; this is a documented departure from fully fold-internal preprocessing.")
  writeLines(L, file.path(RES_DIR, "sens_headline.md"))
  invisible(L)
}

if (!exists("SENS_NO_RUN")) {
  out <- run_all_sections()
  rows <- assemble_outputs(out$secs, out$total_minutes)
  if (!QUICK) write_headline(rows, out)
  message(sprintf("04_sensitivity.R complete: %.1f min total, %d cores, %d tidy rows", out$total_minutes, CORES, nrow(rows)))
}
