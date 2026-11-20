# lib/prep.R : fold-internal preprocessing (fitted on TRAINING rows only).
#
# API (stable; read by T3/T4)
#   prep_fit(train, predictors, spline_vars = NULL)
#       train      : data.frame/tibble/data.table with the predictor columns (outcome not needed)
#       predictors : character vector; must already have passed check_predictors()
#       spline_vars: optional numeric predictors that get restricted-cubic-spline (rcs) nonlinear terms
#       Fits on train only: numeric/binary median (NA -> median), factor levels (NA -> level "Missing"),
#       column means/sds (for scaling), spline knots (quantiles of train), constant-column drop.
#       Returns an object of class "asd_prep".
#   prep_apply(obj, df, scale = TRUE, spline = FALSE)
#       Applies the fitted object to ANY rows (train or test). Returns a numeric matrix with the SAME
#       columns in the same order for every call. Factor levels unseen in train become all-zero dummies.
#       scale=TRUE standardises with the TRAIN mean/sd (glmnet, nnet, lr); spline=TRUE appends rcs terms.
#   check_predictors(predictors, ps, allow_sensitivity_only = FALSE)
#       Leakage guard: stops unless every predictor is in ps$common (or, if allowed, ps$sensitivity_only).
#   make_folds(y, k, seed)      : stratified fold ids (1..k), local seed (global RNG state is restored)
#   with_seed(seed, expr)       : evaluate expr under set.seed(seed), restore previous RNG state
#   load_predictor_sets(V)      : readRDS of derived/predictor_sets.rds (cached)
#   PARSIMONIOUS_PREDICTORS     : the prespecified <=10-parameter predictor set used for Tier 3 (cardiac)
#   RESTRICTED_PREDICTORS(ps)   : the "no labs" restricted set used for Tier 2 (SSI)

.asd_cache <- new.env(parent = emptyenv())

load_predictor_sets <- function(V = get0("V", envir = globalenv())) {
  if (is.null(.asd_cache$ps)) .asd_cache$ps <- readRDS(file.path(V, "derived", "predictor_sets.rds"))
  .asd_cache$ps
}

# Tier 3 prespecified parsimonious set (clinical rationale, chosen without looking at outcome data):
# classic cardiac-risk factors. Dummy-coded parameters: age 1, sex 1, bmi 1, heart failure 1,
# hypertension 1, COPD 1, diabetes 2 (non-insulin, insulin), dialysis 1, smoker 1 = 10.
PARSIMONIOUS_PREDICTORS <- c("age", "sex_male", "bmi", "heart_failure", "hypertension", "copd",
                             "diabetes", "dialysis_preop", "smoker")
# Tier 2 restricted set (SSI): common set without the 8 laboratory values and their indicators (46 parameters).
RESTRICTED_PREDICTORS <- function(ps) setdiff(ps$common, c(ps$labs, ps$missing_indicators))

check_predictors <- function(predictors, ps = load_predictor_sets(), allow_sensitivity_only = FALSE) {
  ok <- ps$common
  if (allow_sensitivity_only) ok <- c(ok, ps$sensitivity_only)
  bad <- setdiff(predictors, ok)
  if (length(bad)) stop("Leakage guard: predictors outside predictor_sets$common: ", paste(bad, collapse = ", "))
  if (anyDuplicated(predictors)) stop("Duplicated predictors")
  invisible(TRUE)
}

with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv())
  on.exit({
    if (had) assign(".Random.seed", old, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  }, add = TRUE)
  set.seed(seed)
  force(expr)
}

make_folds <- function(y, k, seed) {
  with_seed(seed, {
    f <- integer(length(y))
    for (cl in unique(y)) {
      idx <- which(y == cl)
      idx <- idx[sample.int(length(idx))]
      f[idx] <- rep_len(sample.int(k), length(idx))
    }
    f
  })
}

# restricted cubic spline basis (Harrell), k knots -> k-2 nonlinear columns
.rcs_basis <- function(x, knots) {
  k <- length(knots)
  pos3 <- function(z) pmax(z, 0)^3
  tk <- knots[k]; tk1 <- knots[k - 1]
  out <- matrix(0, length(x), k - 2)
  for (j in seq_len(k - 2)) {
    out[, j] <- (pos3(x - knots[j]) - pos3(x - tk1) * (tk - knots[j]) / (tk - tk1) +
                   pos3(x - tk) * (tk1 - knots[j]) / (tk - tk1)) / (tk - knots[1])^2
  }
  out
}

.is_cat <- function(x) is.factor(x) || is.character(x)

.raw_matrix <- function(obj, df) {
  df <- as.data.frame(df)
  n <- nrow(df)
  cols <- list()
  for (v in obj$num_vars) {
    x <- as.numeric(df[[v]])
    x[is.na(x)] <- obj$med[[v]]
    cols[[v]] <- x
  }
  for (v in obj$cat_vars) {
    x <- as.character(df[[v]])
    x[is.na(x)] <- "Missing"
    for (l in obj$levels[[v]][-1]) cols[[paste0(v, "=", l)]] <- as.numeric(x == l)
  }
  m <- matrix(unlist(cols, use.names = FALSE), nrow = n, ncol = length(cols))
  colnames(m) <- names(cols)
  m
}

prep_fit <- function(train, predictors, spline_vars = NULL) {
  train <- as.data.frame(train)
  miss_cols <- setdiff(predictors, names(train))
  if (length(miss_cols)) stop("Predictors missing from train: ", paste(miss_cols, collapse = ", "))
  tr <- train[, predictors, drop = FALSE]
  is_cat <- vapply(tr, .is_cat, logical(1))
  obj <- list(predictors = predictors, cat_vars = predictors[is_cat], num_vars = predictors[!is_cat])
  obj$med <- lapply(obj$num_vars, function(v) { m <- stats::median(tr[[v]], na.rm = TRUE); if (is.na(m)) 0 else m })
  names(obj$med) <- obj$num_vars
  obj$levels <- lapply(obj$cat_vars, function(v) {
    x <- tr[[v]]
    l <- if (is.factor(x)) levels(x) else sort(unique(as.character(x[!is.na(x)])))
    l <- l[l %in% as.character(x)]
    if (anyNA(x)) l <- c(l, "Missing")
    l
  })
  names(obj$levels) <- obj$cat_vars
  raw <- .raw_matrix(obj, tr)
  sds <- apply(raw, 2, stats::sd)
  keep <- which(is.finite(sds) & sds > 0)
  obj$dropped_constant <- colnames(raw)[-keep]
  obj$cols <- colnames(raw)[keep]
  raw <- raw[, keep, drop = FALSE]
  obj$mean <- colMeans(raw)
  obj$sd <- apply(raw, 2, stats::sd)
  # spline knots (train quantiles); fall back to 3 knots, then to linear if too few unique values
  obj$spline_vars <- intersect(spline_vars, obj$num_vars)
  obj$knots <- list()
  for (v in obj$spline_vars) {
    x <- raw[, v]
    kn <- unique(as.numeric(stats::quantile(x, c(0.05, 0.35, 0.65, 0.95), type = 7)))
    if (length(kn) < 4) kn <- unique(as.numeric(stats::quantile(x, c(0.10, 0.50, 0.90), type = 7)))
    if (length(kn) >= 3) obj$knots[[v]] <- kn
  }
  obj$spline_vars <- names(obj$knots)
  if (length(obj$spline_vars)) {
    sp <- .spline_matrix(obj, raw)
    obj$sp_mean <- colMeans(sp)
    obj$sp_sd <- apply(sp, 2, stats::sd)
    bad <- !(is.finite(obj$sp_sd) & obj$sp_sd > 0)
    obj$sp_keep <- colnames(sp)[!bad]
  } else obj$sp_keep <- character(0)
  class(obj) <- "asd_prep"
  obj
}

.spline_matrix <- function(obj, raw) {
  out <- list()
  for (v in obj$spline_vars) {
    b <- .rcs_basis(raw[, v], obj$knots[[v]])
    colnames(b) <- paste0(v, "_rcs", seq_len(ncol(b)))
    out[[v]] <- b
  }
  do.call(cbind, out)
}

prep_apply <- function(obj, df, scale = TRUE, spline = FALSE) {
  raw <- .raw_matrix(obj, df)[, obj$cols, drop = FALSE]
  if (spline && length(obj$spline_vars)) {
    sp <- .spline_matrix(obj, raw)[, obj$sp_keep, drop = FALSE]
    if (scale) sp <- sweep(sweep(sp, 2, obj$sp_mean[obj$sp_keep]), 2, obj$sp_sd[obj$sp_keep], "/")
  } else sp <- NULL
  if (scale) raw <- sweep(sweep(raw, 2, obj$mean), 2, obj$sd, "/")
  if (!is.null(sp)) raw <- cbind(raw, sp)
  raw
}
