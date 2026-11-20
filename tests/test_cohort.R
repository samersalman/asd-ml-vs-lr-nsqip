# test_cohort.R : T1 checks, including the leakage unit test. Run: Rscript tests/test_cohort.R
if (!exists("V")) {
  a <- gsub("~\\+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
  V <- if (length(a)) dirname(dirname(normalizePath(a[1]))) else
    "."
}
suppressPackageStartupMessages({ library(arrow); library(data.table) })
ps <- readRDS(file.path(V, "derived", "predictor_sets.rds"))
ar <- as.data.table(read_parquet(file.path(V, "derived", "analysis_ready.parquet")))
res <- fread(file.path(V, "results", "cohort_predictors.csv"))
fails <- character()
chk <- function(ok, msg) { if (!isTRUE(ok)) fails <<- c(fails, msg); cat(if (isTRUE(ok)) "PASS " else "FAIL ", msg, "\n") }

# ---- leakage: no postoperative or outcome-derived column in the predictor set
leak_tokens <- c("optime", "operative", "wrvu", "workrvu", "rvu", "doptodis", "los", "length", "tothlos", "hdisdt",
  "dischdest", "discharge", "readm", "readmission", "reop", "reoperation", "returnor", "return", "death", "dopertod",
  "expired", "still", "complication", "complications", "supinfec", "wndinfd", "orgspcssi", "ssi", "dehis",
  "dehiscence", "cdmi", "cdarrest", "mi", "arrest", "renainsf", "oprenafl", "intraop", "anesth", "anesthesia",
  "closure", "puf", "year", "era", "operyr", "uid", "caseid", "y", "elig", "days", "postop")
allp <- c(ps$common, ps$sensitivity_only)
toks <- strsplit(tolower(allp), "_")
hits <- allp[vapply(toks, function(t) any(t %in% leak_tokens), logical(1))]
chk(length(hits) == 0, paste("no leakage tokens in predictor names", paste(hits, collapse = ",")))
# source-field map for every predictor: none may come from a postoperative field
src <- c(age = "AGE", bmi = "HEIGHT,WEIGHT", sex_male = "SEX", race_group = "RACE_NEW", hispanic = "ETHNICITY_HISPANIC",
  functional_status = "FNSTATUS2", asa_class = "ASACLAS", inpatient = "INOUT", admission_origin = "TRANST",
  diabetes = "DIABETES", smoker = "SMOKE", copd = "HXCOPD", heart_failure = "HXCHF", hypertension = "HYPERMED",
  dialysis_preop = "DIALYSIS", disseminated_cancer = "DISCANCR", steroid_immunosuppression = "STEROID",
  bleeding_disorder = "BLEEDDIS", preop_transfusion = "TRANSFUS", preop_sirs = "PRSEPIS",
  lab_sodium = "PRSODM", lab_bun = "PRBUN", lab_creatinine = "PRCREAT", lab_albumin = "PRALBUM",
  lab_hematocrit = "PRHCT", lab_platelets = "PRPLATE", lab_wbc = "PRWBC", lab_inr = "PRINR",
  so_weight_loss = "WTLOSS", so_open_wound = "WNDINF", so_dyspnea = "DYSPNEA", so_wound_class = "WNDCLAS")
cpt_derived <- res[cpt_derived == "yes", variable]
forbidden_src <- c("OPTIME", "WORKRVU", "OTHERWRVU", "CONWRVU", "ANESTHES", "WOUND_CLOSURE", "DOPTODIS", "TOTHLOS",
  "HDISDT", "HTOODAY", "DISCHDEST", "STILLINHOSP", "DOPERTOD", "READMISSION1", "RETURNOR", "REOPERATION1",
  "SUPINFEC", "WNDINFD", "ORGSPCSSI", "DEHIS", "CDMI", "CDARREST", "RENAINSF", "OPRENAFL", "UNPLANNEDREADMISSION1")
mapped <- unlist(strsplit(unname(src), ","))
chk(!any(mapped %in% forbidden_src), "no predictor sourced from a postoperative field")
covered <- c(names(src), ps$missing_indicators, cpt_derived)
chk(all(allp %in% covered), paste("every predictor has a documented source:", paste(setdiff(allp, covered), collapse = ",")))
# outcome columns are never predictors
chk(!any(ps$outcomes %in% allp), "outcome columns are not predictors")
chk(!any(c("los_days", "death30", "puf_year", "era", "case_uid") %in% allp), "los_days, death30, year, era, id not predictors")
# CPT-derived plan variables only use CPT codes, never wRVU or time fields (check the build script text)
code <- readLines(file.path(V, "analysis", "01_build_cohort.R"))
pred_lines <- code[grepl("hasc\\(|n_cpt", code)]
chk(!any(grepl("WORKRVU|OPTIME|WRVU", pred_lines)), "CPT-derived predictors reference no wRVU/time field")

# ---- structure
chk(all(allp %in% names(ar)), "all predictors present in analysis-ready file")
chk(!anyDuplicated(ar$case_uid), "one row per case_uid (one surgery per row)")
chk(all(ar$in_primary | ar$in_broad | ar$in_strict | ar$in_primary_nm), "every row belongs to at least one cohort")
chk(all(ar[in_primary == TRUE, in_broad]), "primary is a subset of broad")
chk(all(ar[in_strict == TRUE, in_broad]), "strict is a subset of broad")
chk(all(ar[in_primary == TRUE, in_primary_nm]), "primary is a subset of primary + neuromuscular")
for (y in ps$outcomes) chk(all(ar[[y]] %in% c(0L, 1L, NA)), paste(y, "is 0/1/NA"))
chk(all(is.na(ar[puf_year >= 2021, y_renal])), "renal outcome NA from 2021")
chk(all(ar[in_primary == TRUE & !is.na(y_renal), dialysis_preop] == 0L), "renal denominator has no preop dialysis")
chk(min(ar$age) >= 18 && !anyNA(ar$age), "age >= 18 and non-missing")
chk(all(ar$sex_male %in% c(0L, 1L)), "sex male/female only")
chk(isTRUE(ps$los_cutoff_days >= 1), "LOS cutoff recorded")
fl <- fread(file.path(V, "results", "cohort_flow.csv"))
chk(fl[section == "cohort" & cohort == "primary", tail(n_remaining, 1)] == sum(ar$in_primary), "flow end count equals primary n")
chk(all(fl[section == "cohort", n_removed][-1] >= 0, na.rm = TRUE), "no negative flow steps")
ti <- fread(file.path(V, "results", "tiers.csv"))
chk(all(ti[events < 100 & outcome != "los_ge9", tier] == "Descriptive only"), "tiers: <100 events are descriptive only")
chk(file.exists(file.path(V, "results", "sample_size.csv")), "sample_size.csv exists")
# ---- T1 fix additions (D1, D2)
chk("in_dx_required" %in% names(ar), "in_dx_required column present")
chk(all(ar[in_dx_required == TRUE, in_primary]), "dx_required is a subset of primary")
chk(all(ar[in_dx_required == TRUE, dx_def]), "dx_required cases all have a deformity diagnosis")
chk(sum(ar$in_dx_required) == sum(ar$in_primary & ar$dx_def), "dx_required = primary with deformity diagnosis")
chk(fl[section == "cohort" & cohort == "dx_required", tail(n_remaining, 1)] == sum(ar$in_dx_required), "flow end count equals dx_required n")
chk(all(c("primary", "dx_required") %in% ti$cohort), "tiers.csv has primary and dx_required rows")
cm <- fread(file.path(V, "results", "cohort_case_mix.csv"))
chk(sum(cm$n) == sum(ar$in_primary), "case mix table sums to primary n")
req_cpt <- c("n_cpt_codes", "interbody_fusion", "three_column_osteotomy", "posterior_column_osteotomy", "pelvic_fixation",
             "instrumentation", "post_deformity_arthrodesis", "ant_deformity_arthrodesis", "combined_ap")
chk(all(res[variable %in% req_cpt, cpt_derived] == "yes") && all(req_cpt %in% res$variable), "cpt_derived = yes for every CPT-derived predictor")
chk(res[cpt_derived == "yes", .N] == length(req_cpt), "no extra variables flagged cpt_derived")
ev <- fread(file.path(V, "results", "events_by_outcome.csv"))
chk(ev[cohort == "primary" & outcome == "nonhome", events] == 2448, "primary non-home events unchanged (2448)")
chk(sum(ar$in_primary) == 7699, "primary n unchanged (7699)")
if (length(fails)) { cat("\nFAILED:", length(fails), "\n"); quit(status = 1) } else cat("\nALL TESTS PASSED\n")
