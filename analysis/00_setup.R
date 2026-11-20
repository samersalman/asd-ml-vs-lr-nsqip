# 00_setup.R : paths, seed, packages, read-only NSQIP DuckDB helper.
# Source this first in every analysis script: source(file.path(V, "analysis", "00_setup.R"))
# or, when run from run_all.R, V is already defined.

if (!exists("V")) {
  V <- Sys.getenv("ASD_V", unset = ".")
}
stopifnot(dir.exists(V))

P <- list(
  root = V,
  analysis = file.path(V, "analysis"),
  lib = file.path(V, "lib"),
  derived = file.path(V, "derived"),
  results = file.path(V, "results"),
  figures = file.path(V, "figures"),
  tables = file.path(V, "tables"),
  manuscript = file.path(V, "manuscript"),
  qa = file.path(V, "qa"),
  refs = file.path(V, "refs"),
  tests = file.path(V, "tests"),
  db = Sys.getenv("NSQIP_DUCKDB", unset = "nsqip.duckdb")
)
for (d in c("derived", "results", "figures", "tables", "manuscript", "qa")) {
  dir.create(P[[d]], showWarnings = FALSE, recursive = TRUE)
}

SEED <- 20261116L
set.seed(SEED)
STUDY_YEARS <- 2012:2024

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
  library(data.table)
})

# Optional packages are attached by the scripts that need them, e.g.
# xgboost, ranger, nnet, rpart, rms, dcurves, pmsampsize, pROC, ggplot2, ragg, officer, flextable.

# Read-only connection to the curated NSQIP view `main` (Adult PUF 2010-2024).
# Identity is (CASEID, puf_year) or case_uid. Always filter on puf_year.
nsqip_con <- function() {
  DBI::dbConnect(duckdb::duckdb(), dbdir = P$db, read_only = TRUE)
}

# Run a query on a short-lived read-only connection and return a data.table.
nsqip_query <- function(sql) {
  con <- nsqip_con()
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  data.table::as.data.table(DBI::dbGetQuery(con, sql))
}

# Pre-operative CPT slots to scan for cohort definition: principal CPT plus
# OTHERCPT1-10 and CONCPT1-10 (all VARCHAR).
CPT_SLOTS <- c("CPT", paste0("OTHERCPT", 1:10), paste0("CONCPT", 1:10))
