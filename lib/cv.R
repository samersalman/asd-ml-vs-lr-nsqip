# lib/cv.R : repeated stratified nested CV, replication split, hold-out (era) runs.  Requires lib/prep.R, lib/fit.R.
#
# API (stable; read by T3/T4)
#   load_modelling_data(V, cohort = "in_primary")   -> list(d = data.frame (cohort rows only), ps = predictor sets)
#   outcome_data(d, outcome)                         -> rows of d with non-NA outcome
#   make_outer_folds(y, k = 5, repeats = 5, seed)    -> integer matrix n x repeats of stratified fold ids
#   run_nested_cv(d, outcome, models, predictors = NULL, repeats = 5, k = 5, inner_k = 5, tune = TRUE, seed = SEED,
#                 cores = 1, task_dir = NULL, resume = FALSE, ps = NULL, label = NULL, verbose = TRUE)
#       d must already be outcome-complete (use outcome_data). Every (model, repeat, fold) is one task, run with
#       dynamic scheduling over `cores` forked workers; each task is deterministic from (outcome, repeat, fold, seed),
#       independent of cores/ordering. task_dir (optional) stores one .rds per task so an interrupted run resumes.
#       Returns list(oof = data.table(outcome, model, arm="cv", repeat, fold, case_uid, y, p),
#                    seed_preds = data.table (mlp only; per-seed p1..p5, same keys),
#                    hparams = data.table(outcome, model, repeat, fold, hparams (string), tune_score, converged, notes, seconds, n_train, n_events),
#                    runlog = data.table(outcome, model, n_tasks, n_failed, cpu_seconds, n_not_converged, notes_summary))
#       Each repeat gives every patient exactly one out-of-fold prediction per model (pooled OOF per repeat).
#   run_split_arm(d, outcome, models, predictors = NULL, seed, tune = TRUE, cores = 1, ps = NULL, label = NULL)
#       Replication arm: ONE stratified 60/20/20 split (train/val/test); tuning on the validation set, evaluation on
#       the test set. Same return layout, arm = "split", repeat = 1, fold = 0; only the 20% test rows are returned.
#   run_holdout(d, outcome, models, train_idx, test_idx, predictors = NULL, seed, tune = TRUE, cores = 1, tag = NA,
#               ps = NULL, inner_k = 5)
#       Generic train-on-blocks / test-on-held-out run for T4 (era validation, forward chaining): `train_idx` and
#       `test_idx` are row indices into d (outcome-complete). Inner 5-fold CV tuning inside train. Same layout, arm = tag.
#   run_era_validation(d, outcome, models, era_col = "era", ...)  internal-external CV: for each era level, train on the others
#       and test on it; passes through to run_holdout. `fold` column holds the era index, `arm` = "era:<level>".
#   write_oof(res, dir, outcome) / read_oof(dir, outcome)   parquet I/O of res$oof, res$seed_preds
#   Unit of resampling for paired comparisons is the patient (one surgery = one row; NSQIP has no patient id).

suppressPackageStartupMessages({ library(data.table) })
if (!exists("fit_predict")) {
  .vv <- get0("V", envir = globalenv(), ifnotfound = Sys.getenv("ASD_V", unset = "."))
  source(file.path(.vv, "lib", "fit.R"))
}

load_modelling_data <- function(V, cohort = "in_primary") {
  suppressPackageStartupMessages(library(arrow))
  d <- as.data.frame(arrow::read_parquet(file.path(V, "derived", "analysis_ready.parquet")))
  ps <- load_predictor_sets(V)
  stopifnot(cohort %in% names(d))
  list(d = d[which(d[[cohort]] %in% TRUE), , drop = FALSE], ps = ps)
}
outcome_data <- function(d, outcome) { d <- d[!is.na(d[[outcome]]), , drop = FALSE]; rownames(d) <- NULL; d }

.outcome_id <- function(outcome) sum(utf8ToInt(outcome) * seq_along(utf8ToInt(outcome))) %% 9973L

make_outer_folds <- function(y, k = 5, repeats = 5, seed) {
  vapply(seq_len(repeats), function(r) make_folds(y, k, seed + 7919L * r), integer(length(y)))
}
.task_seed <- function(outcome, rep, fold, seed) as.integer((seed + 100003L * .outcome_id(outcome) + 1009L * rep + 17L * fold) %% 2147483000L)

.hp_string <- function(h) paste(names(h), vapply(h, function(v) format(v, digits = 6), ""), sep = "=", collapse = ";")

.par_lapply <- function(X, FUN, cores) {
  if (cores > 1) parallel::mclapply(X, FUN, mc.cores = cores, mc.preschedule = FALSE, mc.set.seed = FALSE) else lapply(X, FUN)
}

.run_one <- function(model, train, test, outcome, seed, tune, predictors, val, inner_k, ps) {
  tryCatch({
    p <- fit_predict(model, train, test, outcome, seed, tune = tune, predictors = predictors, val = val, inner_k = inner_k, ps = ps)
    list(ok = TRUE, p = as.numeric(p), seed_preds = attr(p, "seed_preds"),
         hp = .hp_string(attr(p, "hparams")), score = attr(p, "tune_score"), converged = attr(p, "converged"),
         notes = paste(attr(p, "notes"), collapse = " | "), seconds = attr(p, "seconds"),
         n_train = attr(p, "n_train"), n_events = attr(p, "n_events"))
  }, error = function(e) list(ok = FALSE, p = rep(NA_real_, nrow(test)), seed_preds = NULL, hp = NA_character_,
                              score = NA_real_, converged = FALSE, notes = paste("ERROR:", conditionMessage(e)),
                              seconds = NA_real_, n_train = nrow(train), n_events = sum(train[[outcome]])))
}

.assemble <- function(tasks, results, d, outcome, arm_fun) {
  oof <- list(); sp <- list(); hp <- list()
  for (i in seq_along(results)) {
    r <- results[[i]]; t <- tasks[[i]]
    if (inherits(r, "try-error") || is.null(r)) r <- list(ok = FALSE, p = rep(NA_real_, length(t$test_idx)), notes = "task crashed", hp = NA_character_,
                                                         score = NA_real_, converged = FALSE, seconds = NA_real_, n_train = NA, n_events = NA, seed_preds = NULL)
    ti <- t$test_idx
    key <- data.table(outcome = outcome, model = t$model, arm = arm_fun(t), `repeat` = t$rep, fold = t$fold,
                      case_uid = d$case_uid[ti], y = d[[outcome]][ti])
    oof[[i]] <- cbind(key, p = r$p)
    if (!is.null(r$seed_preds)) { m <- as.data.table(r$seed_preds); setnames(m, paste0("p", seq_len(ncol(m)))); sp[[i]] <- cbind(key, m) }
    hp[[i]] <- data.table(outcome = outcome, model = t$model, `repeat` = t$rep, fold = t$fold, arm = arm_fun(t),
                          hparams = r$hp, tune_score = r$score, converged = r$converged, ok = r$ok, notes = r$notes,
                          seconds = r$seconds, n_train = r$n_train, n_events = r$n_events)
  }
  oof <- rbindlist(oof); hp <- rbindlist(hp); sp <- if (length(sp)) rbindlist(sp) else data.table()
  runlog <- hp[, .(n_tasks = .N, n_failed = sum(!ok), cpu_seconds = sum(seconds, na.rm = TRUE),
                   n_not_converged = sum(!converged), n_with_notes = sum(nzchar(notes) & !is.na(notes)),
                   notes_summary = paste(head(unique(unlist(strsplit(notes[nzchar(notes) & !is.na(notes)], " \\| "))), 3), collapse = " || ")),
               by = .(outcome, model)]
  list(oof = oof, seed_preds = sp, hparams = hp, runlog = runlog)
}

.heavy_first <- function(models) models[order(match(models, c("mlp", "rf", "xgb", "cart", "ridge_lasso", "lr_spline", "lr")))]

run_nested_cv <- function(d, outcome, models, predictors = NULL, repeats = 5, k = 5, inner_k = 5, tune = TRUE,
                          seed = get0("SEED", ifnotfound = 20261116L), cores = 1, task_dir = NULL, resume = FALSE,
                          ps = NULL, label = NULL, verbose = TRUE) {
  d <- as.data.frame(d); y <- d[[outcome]]
  stopifnot(!anyNA(y), all(y %in% 0:1))
  if (is.null(ps)) ps <- load_predictor_sets()
  if (is.null(predictors)) predictors <- ps$common
  check_predictors(predictors, ps)
  fo <- make_outer_folds(y, k, repeats, seed + .outcome_id(outcome))
  tasks <- list()
  for (m in .heavy_first(models)) for (r in seq_len(repeats)) for (f in seq_len(k))
    tasks[[length(tasks) + 1]] <- list(model = m, rep = r, fold = f, test_idx = which(fo[, r] == f), train_idx = which(fo[, r] != f))
  if (!is.null(task_dir)) dir.create(task_dir, showWarnings = FALSE, recursive = TRUE)
  lab <- if (is.null(label)) outcome else label
  tfile <- function(t) if (is.null(task_dir)) NULL else file.path(task_dir, sprintf("%s__%s__r%d_f%d.rds", lab, t$model, t$rep, t$fold))
  worker <- function(t) {
    tf <- tfile(t)
    if (resume && !is.null(tf) && file.exists(tf)) return(readRDS(tf))
    tr <- d[t$train_idx, , drop = FALSE]; te <- d[t$test_idx, , drop = FALSE]
    r <- .run_one(t$model, tr, te, outcome, .task_seed(outcome, t$rep, t$fold, seed), tune, predictors, NULL, inner_k, ps)
    if (!is.null(tf)) saveRDS(r, tf)
    r
  }
  t0 <- proc.time()[["elapsed"]]
  results <- .par_lapply(tasks, worker, cores)
  wall <- proc.time()[["elapsed"]] - t0
  out <- .assemble(tasks, results, d, outcome, function(t) "cv")
  out$runlog[, wall_seconds_all_models := wall]
  if (verbose) message(sprintf("[%s] %d tasks, %d models, wall %.1f min", lab, length(tasks), length(models), wall / 60))
  out
}

run_split_arm <- function(d, outcome, models, predictors = NULL, seed = get0("SEED", ifnotfound = 20261116L),
                          tune = TRUE, cores = 1, ps = NULL, label = NULL) {
  d <- as.data.frame(d); y <- d[[outcome]]
  if (is.null(ps)) ps <- load_predictor_sets()
  if (is.null(predictors)) predictors <- ps$common
  check_predictors(predictors, ps)
  sid <- make_folds(y, 5, seed + 31L + .outcome_id(outcome))        # five stratified parts: 3 train, 1 val, 1 test
  part <- ifelse(sid <= 3, "train", ifelse(sid == 4, "val", "test"))
  tr <- d[part == "train", , drop = FALSE]; va <- d[part == "val", , drop = FALSE]; te_idx <- which(part == "test")
  tasks <- lapply(.heavy_first(models), function(m) list(model = m, rep = 1L, fold = 0L, test_idx = te_idx))
  worker <- function(t) .run_one(t$model, tr, d[t$test_idx, , drop = FALSE], outcome,
                                 .task_seed(outcome, 1L, 0L, seed), tune, predictors, va, 5, ps)
  t0 <- proc.time()[["elapsed"]]
  res <- .par_lapply(tasks, worker, cores)
  out <- .assemble(tasks, res, d, outcome, function(t) "split")
  out$split_sizes <- c(train = nrow(tr), val = nrow(va), test = length(te_idx))
  out
}

run_holdout <- function(d, outcome, models, train_idx, test_idx, predictors = NULL, seed = get0("SEED", ifnotfound = 20261116L),
                        tune = TRUE, cores = 1, tag = NA_character_, ps = NULL, inner_k = 5, fold_id = 0L) {
  d <- as.data.frame(d)
  if (is.null(ps)) ps <- load_predictor_sets()
  if (is.null(predictors)) predictors <- ps$common
  check_predictors(predictors, ps)
  tr <- d[train_idx, , drop = FALSE]
  tasks <- lapply(.heavy_first(models), function(m) list(model = m, rep = 1L, fold = fold_id, test_idx = test_idx))
  worker <- function(t) .run_one(t$model, tr, d[t$test_idx, , drop = FALSE], outcome,
                                 .task_seed(outcome, 1L, fold_id, seed), tune, predictors, NULL, inner_k, ps)
  res <- .par_lapply(tasks, worker, cores)
  .assemble(tasks, res, d, outcome, function(t) tag)
}

run_era_validation <- function(d, outcome, models, era_col = "era", ...) {
  d <- as.data.frame(d); lv <- levels(factor(d[[era_col]]))
  out <- lapply(seq_along(lv), function(i) {
    te <- which(d[[era_col]] == lv[i]); tr <- which(d[[era_col]] != lv[i])
    run_holdout(d, outcome, models, tr, te, tag = paste0("era:", lv[i]), fold_id = i, ...)
  })
  list(oof = rbindlist(lapply(out, `[[`, "oof")), seed_preds = rbindlist(lapply(out, `[[`, "seed_preds"), fill = TRUE),
       hparams = rbindlist(lapply(out, `[[`, "hparams")), runlog = rbindlist(lapply(out, `[[`, "runlog")))
}

write_oof <- function(res, dir, outcome) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  arrow::write_parquet(res$oof, file.path(dir, sprintf("oof_%s.parquet", outcome)))
  if (nrow(res$seed_preds)) arrow::write_parquet(res$seed_preds, file.path(dir, sprintf("oof_%s_mlp_seeds.parquet", outcome)))
  invisible(TRUE)
}
read_oof <- function(dir, outcome) {
  list(oof = as.data.table(arrow::read_parquet(file.path(dir, sprintf("oof_%s.parquet", outcome)))),
       seed_preds = {f <- file.path(dir, sprintf("oof_%s_mlp_seeds.parquet", outcome)); if (file.exists(f)) as.data.table(arrow::read_parquet(f)) else data.table()})
}
