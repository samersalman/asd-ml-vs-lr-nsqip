# 03_metrics.R : T3 evaluation. Pooled out-of-fold metrics, patient-level bootstrap CIs, paired differences vs reference,
# decision curves, calibration/ROC/PR coordinates, split-arm replication (DeLong + Holm).
# Aggregation (see qa/metrics_decisions.md): every metric is computed WITHIN each CV repeat on that repeat's pooled OOF predictions
# (each patient once per repeat) and then averaged over the 5 repeats. Bootstrap: resample patients (case_uid) once per replicate,
# recompute every model x repeat on the SAME resample, average over repeats, percentile 95% CI.
suppressPackageStartupMessages({ library(data.table); library(arrow) })
if (!exists("V")) {
  args_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  V <- if (length(args_file)) dirname(dirname(normalizePath(gsub("~\\+~", " ", args_file[1])))) else
    "."
}
source(file.path(V, "lib", "metrics.R"))

B_BOOT <- 1000L
SEED_M <- 20261116L
MARGIN <- 0.02
RANGES <- list(nonhome = c(0.05, 0.50), los_p75 = c(0.05, 0.50), los_ge9 = c(0.05, 0.50),
               readm_unplanned = c(0.02, 0.20), ssi = c(0.01, 0.10), cardiac = c(0.005, 0.05))
REFS <- list(nonhome = "lr", los_p75 = "lr", readm_unplanned = "lr", los_ge9 = "lr", ssi = "ridge_lasso", cardiac = NA_character_)
OUTCOMES <- names(RANGES)

nb_grid <- function(outcome, n) round(seq(RANGES[[outcome]][1], RANGES[[outcome]][2], length.out = n), 4)

stat_fun_factory <- function(thr) function(y, p) {
  m <- all_metrics(y, p, thresholds = thr, topk = c(0.1, 0.2))
  ex <- unlist(lapply(c(0.1, 0.2), function(f) {
    t <- topk_metrics(y, p, f); setNames(t[c("spec", "npv")], paste0("top", round(100 * f), "_", c("spec", "npv")))
  }))
  c(m, ex)
}

# Wide OOF for one outcome/arm: list(y, Plist (one n x M matrix per repeat), models)
get_wide <- function(oof, outcome, arm_) {
  reps <- sort(unique(oof[arm == arm_, `repeat`]))
  w <- lapply(reps, function(r) oof_wide(oof, paste0("y_", outcome), r, arm_))   # file stores outcome as y_<name>
  for (k in seq_along(w)) stopifnot(identical(w[[k]]$case_uid, w[[1]]$case_uid), identical(w[[k]]$y, w[[1]]$y),
                                    identical(colnames(w[[k]]$P), colnames(w[[1]]$P)))
  list(y = w[[1]]$y, Plist = lapply(w, `[[`, "P"), models = colnames(w[[1]]$P), reps = reps, case_uid = w[[1]]$case_uid)
}

# Repeat-averaged point estimates (M x S) and bootstrap draws (B x M x S)
avg_stats <- function(y, Plist, stat_fun, ix = NULL) {
  if (is.null(ix)) ix <- seq_along(y)
  yy <- y[ix]; M <- ncol(Plist[[1]])
  arr <- vapply(Plist, function(P) t(vapply(seq_len(M), function(m) stat_fun(yy, P[ix, m]), stat_fun(yy, P[ix, 1]))),
                matrix(0, M, length(stat_fun(yy, Plist[[1]][ix, 1]))))
  out <- apply(arr, c(1, 2), function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE))
  dimnames(out) <- list(colnames(Plist[[1]]), names(stat_fun(yy, Plist[[1]][ix, 1])))
  out
}

boot_repeat_avg <- function(y, Plist, stat_fun, B = B_BOOT, seed = SEED_M, ref = NULL, cores = 1L) {
  n <- length(y); mods <- colnames(Plist[[1]])
  point <- avg_stats(y, Plist, stat_fun); S <- ncol(point); M <- nrow(point); sn <- colnames(point)
  one <- function(b) { set.seed(seed + b); avg_stats(y, Plist, stat_fun, ix = sample.int(n, n, replace = TRUE)) }
  res <- if (cores > 1) parallel::mclapply(seq_len(B), one, mc.cores = cores, mc.preschedule = TRUE) else lapply(seq_len(B), one)
  bad <- vapply(res, function(r) !is.matrix(r), logical(1)); if (any(bad)) stop("bootstrap worker failure")
  draws <- array(NA_real_, c(B, M, S), dimnames = list(NULL, mods, sn))
  for (b in seq_len(B)) draws[b, , ] <- res[[b]]
  q <- function(x) stats::quantile(x, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  ci <- rbindlist(lapply(seq_len(M), function(m) rbindlist(lapply(seq_len(S), function(s) {
    qq <- if (all(is.na(draws[, m, s]))) c(NA_real_, NA_real_) else q(draws[, m, s])
    data.table(model = mods[m], stat = sn[s], est = point[m, s], lo = qq[1], hi = qq[2]) }))))
  diff <- NULL
  if (!is.null(ref) && !is.na(ref) && ref %in% mods)
    diff <- rbindlist(lapply(setdiff(mods, ref), function(m) rbindlist(lapply(seq_len(S), function(s) {
      d <- draws[, m, s] - draws[, ref, s]
      if (all(is.na(d))) return(data.table(model = m, stat = sn[s], est_diff = NA_real_, lo = NA_real_, hi = NA_real_, p_boot = NA_real_))
      qq <- q(d); pb <- 2 * min(mean(d <= 0, na.rm = TRUE), mean(d >= 0, na.rm = TRUE))
      data.table(model = m, stat = sn[s], est_diff = point[m, s] - point[ref, s], lo = qq[1], hi = qq[2], p_boot = min(1, pb)) }))))
  list(point = point, draws = draws, ci = ci, diff = diff)
}

margin_class <- function(lo, hi) ifelse(is.na(lo), NA_character_,
  ifelse(lo > MARGIN, "better_beyond_margin", ifelse(hi < -MARGIN, "worse_beyond_margin",
  ifelse(lo > -MARGIN & hi < MARGIN, "within_margin", "inconclusive"))))

# Display-only curves from predictions pooled across repeats
roc_pr_coords <- function(y, p, npts = 200L) {
  o <- order(p, decreasing = TRUE); ys <- y[o]; ps <- p[o]
  tp <- cumsum(ys); fp <- cumsum(1 - ys); last <- c(ps[-length(ps)] != ps[-1], TRUE)
  tp <- tp[last]; fp <- fp[last]; n1 <- sum(y == 1); n0 <- sum(y == 0)
  idx <- unique(round(seq(1, length(tp), length.out = min(npts, length(tp)))))
  list(roc = data.table(x = c(0, fp[idx] / n0), y = c(0, tp[idx] / n1)),
       pr = data.table(x = c(0, tp[idx] / n1), y = c(1, tp[idx] / (tp[idx] + fp[idx]))))
}
calib_coords <- function(y, p, nbins = 10L, ngrid = 40L) {
  brk <- unique(stats::quantile(p, seq(0, 1, length.out = nbins + 1), names = FALSE))
  out <- list()
  if (length(brk) >= 3) {
    g <- cut(p, brk, include.lowest = TRUE)
    bt <- data.table(g = g, y = y, p = p)[, .(x = mean(p), y = mean(y), n = .N, k = sum(y)), by = g][order(x)]
    z <- 1.96; ph <- bt$y; nn <- bt$n
    bt[, `:=`(lo = (ph + z^2 / (2 * nn) - z * sqrt(ph * (1 - ph) / nn + z^2 / (4 * nn^2))) / (1 + z^2 / nn),
              hi = (ph + z^2 / (2 * nn) + z * sqrt(ph * (1 - ph) / nn + z^2 / (4 * nn^2))) / (1 + z^2 / nn))]
    out$bin <- bt[, .(type = "decile_bin", x, y, n, lo, hi)]
  }
  if (length(unique(p)) >= 5) {
    fit <- tryCatch(suppressWarnings(stats::loess(y ~ p, span = 0.75, degree = 2)), error = function(e) NULL)
    if (!is.null(fit)) {
      r <- stats::quantile(p, c(0.01, 0.99), names = FALSE); gx <- seq(r[1], r[2], length.out = ngrid)
      out$lo <- data.table(type = "loess", x = gx, y = pmin(1, pmax(0, suppressWarnings(stats::predict(fit, data.frame(p = gx))))),
                           n = NA_integer_, lo = NA_real_, hi = NA_real_)
    }
  }
  rbindlist(out, fill = TRUE)
}

run_metrics_main <- function(cores = 12L) {
  Sys.setenv(VECLIB_MAXIMUM_THREADS = "1", OMP_NUM_THREADS = "1")
  derived <- file.path(V, "derived", "oof"); resd <- file.path(V, "results")
  tiers <- fread(file.path(resd, "tiers.csv"))[cohort == "primary", .(outcome, tier)]
  tier_of <- function(o) tiers[outcome == o, tier][1]
  all_rows <- list(); diff_rows <- list(); dca_rows <- list(); dca_cnt <- list(); cal_rows <- list(); rp_rows <- list()
  topk_rows <- list(); delong_rows <- list(); seed_rows <- list(); note_rows <- list()
  for (o in OUTCOMES) {
    message("== ", o)
    oof <- as.data.table(arrow::read_parquet(file.path(derived, sprintf("oof_%s.parquet", o))))
    ref <- REFS[[o]]; tier <- tier_of(o)
    thr_b <- nb_grid(o, 10L); thr_f <- nb_grid(o, 46L)
    for (arm_ in intersect(c("cv", "split"), unique(oof$arm))) {
      w <- get_wide(oof, o, arm_); y <- w$y; n <- length(y); ev <- sum(y)
      message("   arm ", arm_, ": n=", n, " events=", ev, " repeats=", length(w$reps), " models=", paste(w$models, collapse = ","))
      t0 <- proc.time()[3]
      bb <- boot_repeat_avg(y, w$Plist, stat_fun_factory(thr_b), B = B_BOOT, seed = SEED_M, ref = ref, cores = cores)
      message("   bootstrap ", round(proc.time()[3] - t0), " s")
      ci <- copy(bb$ci); ci[, `:=`(outcome = o, tier = tier, arm = arm_, n = n, events = ev)]
      setnames(ci, c("stat", "est"), c("metric", "estimate"))
      all_rows[[length(all_rows) + 1]] <- ci[, .(outcome, tier, model, arm, metric, estimate, lo, hi, n, events)]
      if (!is.null(bb$diff)) {
        dd <- copy(bb$diff); dd[, `:=`(outcome = o, tier = tier, arm = arm_, ref = ref)]
        setnames(dd, c("stat", "est_diff"), c("metric", "diff"))
        dd[, margin_class := ifelse(metric == "auroc", margin_class(lo, hi), NA_character_)]
        dd[, abs_diff_ge_margin := ifelse(metric == "auroc", abs(diff) >= MARGIN, NA)]
        diff_rows[[length(diff_rows) + 1]] <- dd[, .(outcome, tier, arm, ref, model, metric, diff, lo, hi, p_boot, margin_class, abs_diff_ge_margin)]
      }
      # decision curves on the fine grid (repeat-averaged NB, point only) + counts vs reference / treat-all
      nbm <- lapply(w$Plist, function(P) sapply(seq_len(ncol(P)), function(m) net_benefit(y, P[, m], thr_f)))
      nbm <- Reduce(`+`, nbm) / length(nbm); colnames(nbm) <- w$models
      ta <- nb_treat_all(y, thr_f)
      d1 <- rbindlist(lapply(w$models, function(m) data.table(outcome = o, arm = arm_, model = m, threshold = thr_f, nb = nbm[, m])))
      d1 <- rbind(d1, data.table(outcome = o, arm = arm_, model = "treat_all", threshold = thr_f, nb = ta),
                  data.table(outcome = o, arm = arm_, model = "treat_none", threshold = thr_f, nb = 0))
      dca_rows[[length(dca_rows) + 1]] <- d1
      dca_cnt[[length(dca_cnt) + 1]] <- rbindlist(lapply(w$models, function(m) data.table(
        outcome = o, arm = arm_, model = m, ref = ref, range_lo = min(thr_f), range_hi = max(thr_f), n_thresholds = length(thr_f),
        n_beats_ref = if (is.na(ref) || m == ref) NA_integer_ else sum(nbm[, m] > nbm[, ref]),
        n_beats_treat_all = sum(nbm[, m] > ta), n_beats_treat_none = sum(nbm[, m] > 0),
        n_beats_all_and_none = sum(nbm[, m] > pmax(ta, 0)))))
      # calibration + ROC/PR coordinates (pooled over repeats, display only)
      ys <- rep(y, length(w$Plist))
      for (m in w$models) {
        ps <- unlist(lapply(w$Plist, function(P) P[, m]), use.names = FALSE)
        cc <- calib_coords(ys, ps); if (nrow(cc)) cal_rows[[length(cal_rows) + 1]] <- cbind(data.table(outcome = o, arm = arm_, model = m), cc)
        rc <- roc_pr_coords(ys, ps)
        rp_rows[[length(rp_rows) + 1]] <- rbind(cbind(data.table(outcome = o, arm = arm_, model = m, curve = "roc"), rc$roc),
                                                cbind(data.table(outcome = o, arm = arm_, model = m, curve = "pr"), rc$pr))
      }
      # split arm: DeLong (paired) + Holm for AUROC vs reference
      if (arm_ == "split" && !is.na(ref)) {
        P <- w$Plist[[1]]; pv <- numeric(0); rr <- list()
        for (m in setdiff(w$models, ref)) {
          rr[[m]] <- tryCatch({
            r1 <- pROC::roc(y, P[, ref], quiet = TRUE, direction = "<"); r2 <- pROC::roc(y, P[, m], quiet = TRUE, direction = "<")
            tt <- pROC::roc.test(r2, r1, method = "delong", paired = TRUE)
            data.table(model = m, auc_model = as.numeric(pROC::auc(r2)), auc_ref = as.numeric(pROC::auc(r1)),
                       diff = as.numeric(pROC::auc(r2)) - as.numeric(pROC::auc(r1)), z = unname(tt$statistic), p_raw = tt$p.value)
          }, error = function(e) data.table(model = m, auc_model = NA_real_, auc_ref = NA_real_, diff = NA_real_, z = NA_real_, p_raw = NA_real_))
        }
        dl <- rbindlist(rr); dl[, p_holm := p.adjust(p_raw, method = "holm")]
        dl[, `:=`(outcome = o, ref = ref, n = n, events = ev)]
        delong_rows[[length(delong_rows) + 1]] <- dl[, .(outcome, ref, model, auc_model, auc_ref, diff, z, p_raw, p_holm, n, events)]
      }
    }
    # MLP seed spread (cv arm)
    sf <- file.path(derived, sprintf("oof_%s_mlp_seeds.parquet", o))
    if (file.exists(sf)) {
      sp <- as.data.table(arrow::read_parquet(sf)); sp <- sp[arm == "cv"]
      pc <- grep("^p[0-9]+$", names(sp), value = TRUE)
      per <- rbindlist(lapply(split(sp, sp$`repeat`), function(g) {
        v <- vapply(pc, function(cn) auroc(g$y, g[[cn]]), numeric(1))
        as.data.table(c(list(`repeat` = g$`repeat`[1]), as.list(v), list(ensemble = auroc(g$y, rowMeans(as.matrix(g[, pc, with = FALSE]))))))
      }))
      mm <- melt(per, id.vars = "repeat", variable.name = "which", value.name = "auroc")
      sm <- mm[, .(auroc_mean_over_repeats = mean(auroc), auroc_min = min(auroc), auroc_max = max(auroc)), by = which]
      rng <- range(sm[which != "ensemble", auroc_mean_over_repeats]); sm[, `:=`(outcome = o, spread_single_seed_range = rng[2] - rng[1])]
      seed_rows[[length(seed_rows) + 1]] <- sm
    }
    # top-k pulled out of the long tidy file for convenience
  }
  all_dt <- rbindlist(all_rows); diff_dt <- rbindlist(diff_rows)
  fwrite(all_dt, file.path(resd, "metrics_all.csv"))
  fwrite(diff_dt, file.path(resd, "diff_vs_lr.csv"))
  fwrite(all_dt[metric %in% c("auroc", "auprc", "prevalence", "brier", "brier_scaled", "cal_slope", "cal_intercept", "oe_ratio", "ici")],
         file.path(resd, "metrics_core.csv"))
  fwrite(all_dt[grepl("^top(10|20)_", metric)], file.path(resd, "metrics_topk.csv"))
  fwrite(rbindlist(dca_rows), file.path(resd, "dca_net_benefit.csv"))
  cnt <- rbindlist(dca_cnt); fwrite(cnt, file.path(resd, "dca_counts.csv"))
  fwrite(rbindlist(cal_rows), file.path(resd, "calibration_curve_coords.csv"))
  rp <- rbindlist(rp_rows)
  fwrite(rp[arm == "cv"], file.path(resd, "roc_pr_curves_cv.csv")); fwrite(rp[arm == "split"], file.path(resd, "roc_pr_curves_split.csv"))
  if (length(delong_rows)) fwrite(rbindlist(delong_rows), file.path(resd, "metrics_delong_split.csv"))
  if (length(seed_rows)) fwrite(rbindlist(seed_rows, fill = TRUE), file.path(resd, "metrics_mlp_seed_spread.csv"))
  # notes / flags
  rl <- tryCatch(fread(file.path(resd, "models_runlog.csv")), error = function(e) NULL)
  notes <- data.table(outcome = c("readm_unplanned", "readm_unplanned", "readm_unplanned", "nonhome;los_p75;readm_unplanned", "all", "ssi;cardiac", "los_ge9"),
    model = c("lr_spline;cart;rf;xgb;mlp", "cart", "lr;ridge_lasso", "mlp", "all", "all", "lr;xgb"),
    note = c("exploratory: Riley criterion not met for 62 parameters; flexible-model results are unstable and not confirmatory",
             "near root-only tree (few distinct probabilities); ICI is not comparable and is NA when under 5 unique predictions; top-k ties broken by row order",
             "Riley not met at 62 parameters (8.98 events per parameter); interpret with caution",
             "MLP reached nnet maxit=300 in some final seeds (see models_runlog.csv); predictions used as is",
             "los_p75 cutoff (>=8 d) was computed on the whole cohort (75th percentile), not inside training folds",
             "exploratory low-event outcomes; wide CIs",
             "sensitivity comparator definition of prolonged LOS (>=9 d)"))
  fwrite(notes, file.path(resd, "metrics_notes.csv"))
  write_headline(all_dt, diff_dt, cnt, resd)
  invisible(TRUE)
}

write_headline <- function(all_dt, diff_dt, cnt, resd) {
  L <- c("# T3 headline metrics (generated by analysis/03_metrics.R)", "",
         "Primary = pooled out-of-fold predictions from repeated stratified nested CV (5 x 5), metrics computed per repeat then averaged; 95% CIs are percentile patient-level bootstrap (B = 1000). Margin for AUROC differences = 0.02. Counts only; no claim of superiority is made from point estimates.", "")
  g <- function(o, mod, a, met) { r <- all_dt[outcome == o & model == mod & arm == a & metric == met]; if (nrow(r)) sprintf("%.3f (%.3f-%.3f)", r$estimate, r$lo, r$hi) else NA }
  for (o in OUTCOMES) {
    r <- all_dt[outcome == o & arm == "cv" & metric == "auroc"]
    L <- c(L, sprintf("## %s (tier: %s; n=%d, events=%d, prevalence %.4f)", o, r$tier[1], r$n[1], r$events[1], r$events[1] / r$n[1]), "",
           "| model | AUROC (95% CI) | AUPRC | Brier scaled | slope | O:E | ICI |", "|---|---|---|---|---|---|---|")
    for (m in r$model) L <- c(L, sprintf("| %s | %s | %s | %s | %s | %s | %s |", m, g(o, m, "cv", "auroc"), g(o, m, "cv", "auprc"),
                                          g(o, m, "cv", "brier_scaled"), g(o, m, "cv", "cal_slope"), g(o, m, "cv", "oe_ratio"), g(o, m, "cv", "ici")))
    d <- diff_dt[outcome == o & arm == "cv" & metric == "auroc"]
    if (nrow(d)) {
      L <- c(L, "", sprintf("AUROC difference vs %s (CV arm):", d$ref[1]), "")
      for (i in seq_len(nrow(d))) L <- c(L, sprintf("- %s: %+.3f (%.3f to %.3f), class = %s", d$model[i], d$diff[i], d$lo[i], d$hi[i], d$margin_class[i]))
    }
    L <- c(L, "")
  }
  d <- diff_dt[arm == "cv" & metric == "auroc"]
  L <- c(L, "## Counts against the 0.02 AUROC margin (CV arm, all non-reference comparisons)", "",
         sprintf("- comparisons: %d", nrow(d)),
         sprintf("- point-estimate |diff| >= 0.02: %d (of which positive: %d)", sum(d$abs_diff_ge_margin), sum(d$abs_diff_ge_margin & d$diff > 0)),
         sprintf("- whole 95%% CI above +0.02 (better beyond margin): %d", sum(d$margin_class == "better_beyond_margin")),
         sprintf("- whole 95%% CI below -0.02 (worse beyond margin): %d", sum(d$margin_class == "worse_beyond_margin")),
         sprintf("- whole 95%% CI inside +/-0.02 (equivalent within margin): %d", sum(d$margin_class == "within_margin")),
         sprintf("- inconclusive: %d", sum(d$margin_class == "inconclusive")), "")
  if (nrow(d)) L <- c(L, "Comparisons whose point estimate exceeded the margin:", "",
                      sprintf("- %s: %s vs %s %+.3f (%.3f to %.3f) [%s]", d[abs_diff_ge_margin == TRUE, outcome], d[abs_diff_ge_margin == TRUE, model],
                              d[abs_diff_ge_margin == TRUE, ref], d[abs_diff_ge_margin == TRUE, diff], d[abs_diff_ge_margin == TRUE, lo],
                              d[abs_diff_ge_margin == TRUE, hi], d[abs_diff_ge_margin == TRUE, margin_class]), "")
  c1 <- cnt[arm == "cv" & !is.na(n_beats_ref)]
  L <- c(L, "## Decision curves (CV arm): thresholds where model net benefit exceeds the reference (counts, grid of 46)", "")
  L <- c(L, sprintf("- %s / %s: %d of %d (beats treat-all and none: %d)", c1$outcome, c1$model, c1$n_beats_ref, c1$n_thresholds, c1$n_beats_all_and_none), "")
  dl <- file.path(resd, "metrics_delong_split.csv")
  if (file.exists(dl)) { x <- fread(dl); L <- c(L, "## Split replication arm (DeLong, Holm-adjusted within outcome)", "",
    sprintf("- %s %s vs %s: dAUROC %+.3f, Holm p = %.3f", x$outcome, x$model, x$ref, x$diff, x$p_holm), "") }
  L <- c(L, "## Caveats", "",
         "- Readmission flexible models (lr_spline, cart, rf, xgb, mlp) are exploratory: Riley criterion not met; unpenalised LR/MLP/RF slopes well below 1.",
         "- MLP reached nnet maxit = 300 in some final seeds; predictions used as is (seed-average).",
         "- CART readmission is near root-only; its ICI is not comparable (NA when fewer than 5 distinct predictions) and top-k ties are broken by row order.",
         "- The LOS cutoff (>= 8 d) was computed on the whole cohort, not inside training folds.",
         "- SSI and cardiac are exploratory, low-event outcomes; cardiac is a single parsimonious model (performance only).",
         "- Split arm has about 1,520 test patients (about 110 readmission events); imprecise.")
  writeLines(L, file.path(resd, "metrics_headline.md"))
}

if (!exists("METRICS_NO_RUN")) run_metrics_main(cores = min(12L, parallel::detectCores()))
