# run_all.R : sequential driver for the analysis steps 00-04.
# Usage (from the repository root):
#   Rscript run_all.R
# Step 01 needs the ACS NSQIP Adult Participant Use Files loaded into a DuckDB file
# (path in the NSQIP_DUCKDB environment variable). Steps 02-04 use the patient-level
# derived/ files written by step 01; these are never included in this repository.
# Environment knobs of the heavy steps: ASD_CORES, ASD_REPEATS, ASD_RESUME, ASD_OUTCOMES (02);
# ASD_SENS_RESUME, ASD_SENS_QUICK, ASD_BOOT_INSTAB (04).

args_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
V <- if (length(args_file)) dirname(normalizePath(gsub("~\\+~", " ", args_file[1]))) else getwd()
Sys.setenv(ASD_V = V)
setwd(V)
rs <- file.path(R.home("bin"), "Rscript")
steps <- c("00_setup.R", "01_build_cohort.R", "02_models.R", "03_metrics.R", "04_sensitivity.R")
for (s in steps) {
  message(">>> ", s, "  (", format(Sys.time(), "%H:%M:%S"), ")")
  st <- system2(rs, shQuote(file.path(V, "analysis", s)))
  if (st != 0) stop("step failed: ", s, " (exit ", st, ")", call. = FALSE)
}
message("run_all.R finished.")
