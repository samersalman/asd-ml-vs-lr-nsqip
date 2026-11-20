# 02_models.R : nested-CV modelling for every outcome, by tier (primary cohort).
# Inputs : derived/analysis_ready.parquet, derived/predictor_sets.rds, results/tiers.csv
# Outputs: derived/oof/oof_<outcome>.parquet (long OOF predictions, patient level, never commit),
#          derived/oof/oof_<outcome>_mlp_seeds.parquet, results/models_hparams.csv, results/models_runlog.csv,
#          results/models_mlp_seed_spread.csv, results/models_arm_sizes.csv
# Environment knobs: ASD_CORES (default detectCores()-2), ASD_REPEATS (default 5), ASD_RESUME=1, ASD_OUTCOMES (comma list)
# Sourced by run_all.R (V defined) or run directly: Rscript analysis/02_models.R

Sys.setenv(VECLIB_MAXIMUM_THREADS = "1", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
if (!exists("V")) V <- Sys.getenv("ASD_V", unset = ".")
source(file.path(V, "analysis", "00_setup.R"))
suppressPackageStartupMessages(library(arrow))
source(file.path(V, "lib", "prep.R")); source(file.path(V, "lib", "fit.R"))
source(file.path(V, "lib", "cv.R")); source(file.path(V, "lib", "metrics.R"))

CORES   <- as.integer(Sys.getenv("ASD_CORES", unset = max(1, parallel::detectCores() - 2)))
REPEATS <- as.integer(Sys.getenv("ASD_REPEATS", unset = "5"))
RESUME  <- Sys.getenv("ASD_RESUME", unset = "0") == "1"
ONLY    <- Sys.getenv("ASD_OUTCOMES", unset = "")
OOF_DIR <- file.path(P$derived, "oof")
dir.create(OOF_DIR, showWarnings = FALSE, recursive = TRUE)
TASK_DIR <- file.path(OOF_DIR, "tasks")
if (!RESUME && dir.exists(TASK_DIR)) unlink(TASK_DIR, recursive = TRUE)

dat <- load_modelling_data(V, "in_primary")
d_all <- dat$d; ps <- dat$ps
tiers <- fread(file.path(P$results, "tiers.csv"))[cohort == "primary"]

ALL7 <- c("lr", "lr_spline", "ridge_lasso", "cart", "rf", "xgb", "mlp")
# plan: outcome column, label, models, predictors, tier, split arm?, exploratory?, notes
plan <- list(
  list(y = "y_nonhome",         lab = "nonhome",         tier = "Tier 1", models = ALL7, pred = ps$common, split = TRUE,  expl = FALSE, note = ""),
  list(y = "y_los_p75",         lab = "los_p75",         tier = "Tier 1", models = ALL7, pred = ps$common, split = TRUE,  expl = FALSE,
       note = "LOS >= 8 days; Riley met at 15% but not at 10% of max R2 (T1)"),
  list(y = "y_readm_unplanned", lab = "readm_unplanned", tier = "Tier 1", models = ALL7, pred = ps$common, split = TRUE,  expl = TRUE,
       note = "Riley criteria NOT met for 62 parameters (8.98 events/parameter): flexible-model results exploratory"),
  list(y = "y_los_ge9",         lab = "los_ge9",         tier = "Tier 1 comparator", models = c("lr", "xgb"), pred = ps$common, split = FALSE, expl = FALSE,
       note = "LOS >= 9 days comparator, lr and xgb only"),
  list(y = "y_ssi",             lab = "ssi",             tier = "Tier 2", models = c("ridge_lasso", "xgb"), pred = RESTRICTED_PREDICTORS(ps), split = FALSE, expl = TRUE,
       note = "Restricted predictors (no labs, 46 parameters); Riley not met; exploratory"),
  list(y = "y_cardiac",         lab = "cardiac",         tier = "Tier 3", models = c("ridge_lasso"), pred = PARSIMONIOUS_PREDICTORS, split = FALSE, expl = FALSE,
       note = "Prespecified parsimonious penalised LR, 10 parameters (Riley met for 10 parameters); reported as pen_lr_parsimonious")
)
if (nzchar(ONLY)) plan <- Filter(function(p) p$lab %in% strsplit(ONLY, ",")[[1]], plan)
# renal (Descriptive only) deliberately skipped

hp_all <- list(); rl_all <- list(); spread <- list(); sizes <- list()
T_START <- proc.time()[["elapsed"]]
for (pl in plan) {
  trow <- tiers[outcome == sub("^y_", "", pl$y)]
  if (nrow(trow) == 1 && !grepl(sub(" comparator", "", pl$tier), trow$tier) && pl$lab != "los_ge9") stop("tier mismatch for ", pl$lab)
  dd <- outcome_data(d_all, pl$y)
  stopifnot(all(pl$pred %in% ps$common))
  message(sprintf(">>> %s (%s): n=%d events=%d, models=%s", pl$lab, pl$tier, nrow(dd), sum(dd[[pl$y]]), paste(pl$models, collapse = ",")))
  t0 <- proc.time()[["elapsed"]]
  cv <- run_nested_cv(dd, pl$y, pl$models, predictors = pl$pred, repeats = REPEATS, k = 5, inner_k = 5, cores = CORES,
                      task_dir = TASK_DIR, resume = RESUME, ps = ps, label = pl$lab)
  res <- cv
  if (pl$split) {
    sp <- run_split_arm(dd, pl$y, pl$models, predictors = pl$pred, cores = CORES, ps = ps, label = pl$lab)
    res$oof <- rbind(cv$oof, sp$oof)
    res$seed_preds <- rbind(cv$seed_preds, sp$seed_preds, fill = TRUE)
    res$hparams <- rbind(cv$hparams, sp$hparams)
    res$runlog <- rbind(cv$runlog, sp$runlog, fill = TRUE)
    sizes[[length(sizes) + 1]] <- data.table(outcome = pl$lab, as.data.table(as.list(sp$split_sizes)))
  }
  wall <- proc.time()[["elapsed"]] - t0
  if (pl$lab == "cardiac") for (nm in c("oof", "seed_preds", "hparams", "runlog")) if (nrow(res[[nm]])) res[[nm]][model == "ridge_lasso", model := "pen_lr_parsimonious"]
  write_oof(res, OOF_DIR, pl$lab)
  res$hparams[, outcome := pl$lab]; res$runlog[, outcome := pl$lab]; res$oof[, outcome := pl$lab]
  # exploratory flag and metadata in the run log
  rl <- res$runlog
  rl[, arm := ifelse(is.na(wall_seconds_all_models), "split", "cv")]
  rl[, `:=`(tier = pl$tier, n = nrow(dd), events = sum(dd[[pl$y]]), n_predictors = length(pl$pred), repeats = REPEATS,
            exploratory = pl$expl, outcome_note = pl$note, cores = CORES, wall_seconds_outcome = wall)]
  rl[pl$lab == "readm_unplanned" & model %in% c("lr_spline", "cart", "rf", "xgb", "mlp"), exploratory := TRUE]
  rl[pl$lab == "readm_unplanned" & model %in% c("lr", "ridge_lasso"), exploratory := FALSE]
  rl_all[[length(rl_all) + 1]] <- rl
  hp_all[[length(hp_all) + 1]] <- res$hparams
  # mlp per-seed AUROC spread (CV arm)
  if (nrow(res$seed_preds)) {
    sps <- res$seed_preds[arm == "cv"]
    pc <- grep("^p[0-9]+$", names(sps), value = TRUE)
    spread[[length(spread) + 1]] <- rbindlist(lapply(unique(sps$`repeat`), function(r) {
      s <- sps[`repeat` == r]
      a <- vapply(pc, function(cn) auroc(s$y, s[[cn]]), numeric(1))
      m <- res$oof[arm == "cv" & model == "mlp" & `repeat` == r]
      data.table(outcome = pl$lab, `repeat` = r, auroc_mean_of_seeds_prediction = auroc(m$y, m$p),
                 auroc_seed_min = min(a), auroc_seed_max = max(a), auroc_seed_sd = sd(a))
    }))
  }
  fails <- res$runlog[n_failed > 0]
  if (nrow(fails)) message("!!! FAILED TASKS: ", paste(fails$outcome, fails$model, fails$n_failed, collapse = "; "))
  message(sprintf("    done %s in %.1f min", pl$lab, wall / 60))
  rm(res, cv); invisible(gc())
}
TOTAL <- proc.time()[["elapsed"]] - T_START
hp <- rbindlist(hp_all, fill = TRUE)
rl <- rbindlist(rl_all, fill = TRUE)
rl[, total_pipeline_minutes := round(TOTAL / 60, 1)]
if (!nzchar(ONLY)) {
  fwrite(hp[, .(outcome, model, arm, `repeat`, fold, hparams, tune_score = round(tune_score, 5), converged, ok, seconds = round(seconds, 2), n_train, n_events, notes)], file.path(P$results, "models_hparams.csv"))
  fwrite(rl[, .(outcome, model, arm, tier, exploratory, n, events, n_predictors, repeats, n_tasks, n_failed, n_not_converged, n_with_notes,
                cpu_seconds = round(cpu_seconds, 1), wall_seconds_all_models = round(wall_seconds_all_models, 1), wall_seconds_outcome = round(wall_seconds_outcome, 1),
                cores, total_pipeline_minutes, outcome_note, notes_summary)], file.path(P$results, "models_runlog.csv"))
  if (length(spread)) fwrite(rbindlist(spread), file.path(P$results, "models_mlp_seed_spread.csv"))
  if (length(sizes)) fwrite(rbindlist(sizes), file.path(P$results, "models_arm_sizes.csv"))
}
message(sprintf("02_models.R complete: %.1f min total, %d cores", TOTAL / 60, CORES))
