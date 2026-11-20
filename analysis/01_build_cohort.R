# 01_build_cohort.R : T1 cohort build + Phase 0 gate (ASD ML vs LR, ACS-NSQIP 2012-2024)
# Run via source() from run_all.R (V defined) or: Rscript analysis/01_build_cohort.R
# Writes: derived/analysis_ready.parquet (+ .rds), derived/predictor_sets.rds,
#         results/cohort_*.csv, events_*.csv, tiers.csv, sample_size.csv, cohort_summary.md
# No model is fit here. Tiers are assigned from counts only.

if (!exists("V")) {
  args_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  args_file <- gsub("~\\+~", " ", args_file)
  V <- if (length(args_file)) dirname(dirname(normalizePath(args_file[1]))) else
    "."
}
source(file.path(V, "analysis", "00_setup.R"))
suppressPackageStartupMessages({ library(arrow); library(pmsampsize) })
set.seed(SEED)

# ---------------------------------------------------------------- constants
DEF_ARTH <- c("22800", "22802", "22804", "22808", "22810", "22812")
LONG_INSTR <- c("22843", "22844")
CAND_CODES <- c(DEF_ARTH, LONG_INSTR)
GE7_CODES <- c("22802", "22804", "22812", "22843", "22844", "22206", "22207", "22848")
CERV_CODES <- c("22548", "22551", "22552", "22554", "22585", "22590", "22595", "22600",
                "22856", "22858", "22318", "22319", "63001", "63015", "63020", "63040",
                "63045", "63050", "63051", "63075", "63076", "63081", "63082", "22210")
TL_CODES <- c("22533", "22534", "22556", "22558", "22610", "22612", "22630", "22632",
              "22633", "22634", "22206", "22207", "22208", "22212", "22214", "22216",
              "22800", "22802", "22804", "22808", "22810", "22812", "22830", "22848",
              "63005", "63011", "63012", "63017", "63030", "63035", "63046", "63047",
              "63048", "63055", "63056", "63057", "63064", "63085", "63087", "63090")
# note: 22800-22812 and 22843/22844 are region-agnostic in CPT; deformity arthrodesis codes
# are counted as thoracolumbar-compatible so that a cervical-only rule cannot remove them.
ANT_CODES <- c("22808", "22810", "22812", "22558", "22556", "22533", "22534")
POST_CODES <- c("22800", "22802", "22804", "22610", "22612", "22614", "22630", "22632",
                "22633", "22634", "22840", "22842", "22843", "22844", "22848")
ERA_LEVELS <- c("2012-2015", "2016-2018", "2019-2021", "2022-2024")
era_of <- function(y) factor(ifelse(y <= 2015, ERA_LEVELS[1], ifelse(y <= 2018, ERA_LEVELS[2],
                      ifelse(y <= 2021, ERA_LEVELS[3], ERA_LEVELS[4]))), levels = ERA_LEVELS)
outcome_csv <- function(x, name) fwrite(x, file.path(P$results, name))

# ---------------------------------------------------------------- pull candidates
slots_cols <- CPT_SLOTS
slot_pred <- paste(sprintf("TRIM(%s) IN ('%s')", slots_cols, paste(CAND_CODES, collapse = "','")),
                   collapse = " OR ")
sel <- c("case_uid", "CASEID", "puf_year", "AGE", "SEX", "RACE_NEW", "ETHNICITY_HISPANIC",
         "HEIGHT", "WEIGHT", "FNSTATUS2", "ASACLAS", "TRANST", "INOUT", "DIABETES", "SMOKE",
         "HXCOPD", "HXCHF", "HYPERMED", "DIALYSIS", "DISCANCR", "STEROID", "BLEEDDIS",
         "TRANSFUS", "PRSEPIS", "VENTILAT", "PRSODM", "PRBUN", "PRCREAT", "PRALBUM", "PRHCT",
         "PRPLATE", "PRWBC", "PRINR", "WTLOSS", "WNDINF", "DYSPNEA", "WNDCLAS",
         "ELECTSURG", "EMERGNCY", "CASETYPE", "PODIAG", "PODIAG10",
         "DISCHDEST", "DOPTODIS", "STILLINHOSP", "DOPERTOD", "READMISSION1",
         "UNPLANNEDREADMISSION1", "READMPODAYS1", "SUPINFEC", "WNDINFD", "ORGSPCSSI",
         "SSSIPATOS", "DSSIPATOS", "OSSIPATOS", "DEHIS", "CDMI", "CDARREST", "RENAINSF",
         "OPRENAFL", slots_cols)
sql <- sprintf("SELECT %s FROM main WHERE puf_year BETWEEN 2012 AND 2024 AND (%s)",
               paste(unique(sel), collapse = ", "), slot_pred)
D <- nsqip_query(sql)
setorder(D, puf_year, CASEID)
stopifnot(!anyDuplicated(D$case_uid))
chr_cols <- names(D)[vapply(D, is.character, logical(1))]
for (cc in chr_cols) set(D, j = cc, value = trimws(D[[cc]]))

# ---------------------------------------------------------------- CPT scan (all 21 slots)
slot_mat <- as.matrix(D[, ..slots_cols])
has_any <- function(codes) Reduce(`|`, lapply(slots_cols, function(s) D[[s]] %in% codes))
has_code <- function(code) has_any(code)
D[, `:=`(
  c_def_arth = has_any(DEF_ARTH),
  c_long_instr = has_any(LONG_INSTR),
  c_ge7 = has_any(GE7_CODES),
  c_cerv = has_any(CERV_CODES),
  c_tl = has_any(TL_CODES)
)]
D[, n_cpt := apply(slot_mat, 1, function(r) length(unique(r[!is.na(r) & r != ""])))]

# ---------------------------------------------------------------- diagnosis classes
D[, `:=`(dx = toupper(ifelse(!is.na(PODIAG10), PODIAG10, PODIAG)),
         icd = ifelse(!is.na(PODIAG10), 10L, ifelse(!is.na(PODIAG), 9L, NA_integer_)))]
D[, dx_nd := gsub(".", "", dx, fixed = TRUE)]
rx <- function(p, x) !is.na(x) & grepl(p, x, perl = TRUE)
D[, dx_def := (icd == 10L & rx("^(M40|M41(?!4)|M438|Q763|Q7641)", dx_nd)) |
              (icd == 9L & rx("^737", dx))]
D[is.na(dx_def), dx_def := FALSE]
D[, dx_neuromusc := icd == 10L & rx("^M414", dx_nd)]
D[is.na(dx_neuromusc), dx_neuromusc := FALSE]
D[, dx_tumor := (icd == 10L & rx("^(C|D16|D32|D33|D36|D42|D43|D48)", dx_nd)) |
                (icd == 9L & rx("^(1[4-9][0-9]|2[0-3][0-9])(\\.|$)", dx))]
D[, dx_infect := (icd == 10L & rx("^(M46[2-5]|G061|G062|M86|T814|T846|A180)", dx_nd)) |
                 (icd == 9L & rx("^(730|324)(\\.|$)|^996\\.6|^998\\.5", dx))]
D[, dx_fracture := (icd == 10L & rx("^(S12|S22|S32|M80|M484|M485)", dx_nd)) |
                   (icd == 9L & rx("^(805|806)(\\.|$)|^733\\.1", dx))]
for (cc in c("dx_tumor", "dx_infect", "dx_fracture")) set(D, which(is.na(D[[cc]])), cc, FALSE)

# ---------------------------------------------------------------- base-exclusion flags
D[, elective := ifelse(!is.na(CASETYPE), CASETYPE == "Elective",
                       ELECTSURG %in% "Yes" & fcoalesce(EMERGNCY, "No") == "No")]
D[is.na(elective), elective := FALSE]
D[, non_emergent := ifelse(!is.na(CASETYPE), CASETYPE != "Emergent", fcoalesce(EMERGNCY, "No") == "No")]
D[is.na(non_emergent), non_emergent := FALSE]
D[, sex_ok := SEX %in% c("male", "female")]
D[, era := era_of(puf_year)]

# ---------------------------------------------------------------- cohort rules + flow
base_steps <- list(
  list("Age missing", function(d) !is.na(d$AGE)),
  list("Age < 18 years", function(d) d$AGE >= 18),
  list("Sex missing, non-binary or intersex", function(d) d$sex_ok),
  list("Not elective (ELECTSURG/EMERGNCY to 2020, CASETYPE from 2021)", function(d) d$elective),
  list("Preoperative sepsis or septic shock", function(d) !(d$PRSEPIS %in% c("Sepsis", "Septic Shock"))),
  list("Ventilator dependence", function(d) !(d$VENTILAT %in% "Yes")),
  list("ASA class 5", function(d) !(d$ASACLAS %in% "5-Moribund")),
  list("Principal diagnosis: tumour", function(d) !d$dx_tumor),
  list("Principal diagnosis: spinal infection", function(d) !d$dx_infect),
  list("Principal diagnosis: acute vertebral fracture", function(d) !d$dx_fracture),
  list("Cervical-only procedure", function(d) !(d$c_cerv & !d$c_tl))
)
elig <- list(
  primary = list("No deformity arthrodesis code and no (22843/22844 plus deformity diagnosis)",
                 function(d) d$c_def_arth | (d$c_long_instr & d$dx_def)),
  broad = list("No deformity arthrodesis or long-instrumentation code", function(d) d$c_def_arth | d$c_long_instr),
  strict = list("Not (deformity diagnosis and (>=7 segments, 22206/22207 or 22848))",
                function(d) d$dx_def & d$c_ge7),
  primary_nm = list("No deformity arthrodesis code and no (22843/22844 plus deformity diagnosis incl. M41.4x)",
                    function(d) d$c_def_arth | (d$c_long_instr & (d$dx_def | d$dx_neuromusc)))
)
neuro_step <- list("Neuromuscular scoliosis (M41.4x)", function(d) !d$dx_neuromusc)
cohort_steps <- list(
  primary = c(base_steps, list(elig$primary, neuro_step)),
  broad = c(base_steps, list(elig$broad)),
  strict = c(base_steps, list(elig$strict, neuro_step)),
  primary_nm = c(base_steps, list(elig$primary_nm)),
  dx_required = c(base_steps, list(elig$primary, neuro_step,
    list("No deformity diagnosis (arthrodesis-code-only cases removed; D1 prespecified sensitivity)", function(d) d$dx_def)))
)
flow_rows <- list()
masks <- list()
for (nm in names(cohort_steps)) {
  keep <- rep(TRUE, nrow(D))
  flow_rows[[length(flow_rows) + 1]] <- data.table(section = "cohort", cohort = nm, step = 0L,
    description = "Candidate cases: any of 22800/22802/22804/22808/22810/22812/22843/22844 in any of 21 CPT slots, 2012-2024",
    n_removed = NA_integer_, n_remaining = sum(keep))
  i <- 0L
  for (s in cohort_steps[[nm]]) {
    i <- i + 1L
    pass <- s[[2]](D); pass[is.na(pass)] <- FALSE
    removed <- sum(keep & !pass)
    keep <- keep & pass
    flow_rows[[length(flow_rows) + 1]] <- data.table(section = "cohort", cohort = nm, step = i,
      description = s[[1]], n_removed = removed, n_remaining = sum(keep))
  }
  masks[[nm]] <- keep
}
D[, `:=`(in_primary = masks$primary, in_broad = masks$broad, in_strict = masks$strict,
         in_primary_nm = masks$primary_nm, in_dx_required = masks$dx_required)]
A <- D[in_primary | in_broad | in_strict | in_primary_nm]
stopifnot(all(A[in_dx_required == TRUE, in_primary]))
cat("Union analysis rows:", nrow(A), " primary:", sum(A$in_primary), "\n")

# ---------------------------------------------------------------- predictors
num_or_na <- function(x, lo, hi) { x <- as.numeric(x); x[!is.na(x) & (x < lo | x > hi)] <- NA; x }
A[, height_in := num_or_na(HEIGHT, 48, 84)]
A[, weight_lb := num_or_na(WEIGHT, 70, 700)]
A[, bmi := num_or_na(703 * weight_lb / height_in^2, 12, 80)]
A[, age := as.numeric(AGE)]
A[, sex_male := as.integer(SEX == "male")]
A[, race_group := factor(fcase(
  RACE_NEW == "White", "White",
  RACE_NEW == "Black or African American", "Black",
  RACE_NEW == "Asian", "Asian",
  is.na(RACE_NEW) | RACE_NEW %in% c("Unknown/Not Reported", "Hispanic or Latino"), "Unknown",
  default = "Other"), levels = c("White", "Black", "Asian", "Other", "Unknown"))]
A[, hispanic := factor(fcase(ETHNICITY_HISPANIC == "Yes", "Yes", ETHNICITY_HISPANIC == "No", "No",
                             default = "Unknown"), levels = c("No", "Yes", "Unknown"))]
A[, functional_status := factor(fcase(FNSTATUS2 == "Independent", "Independent",
  FNSTATUS2 == "Partially Dependent", "Partially dependent",
  FNSTATUS2 == "Totally Dependent", "Totally dependent", default = "Unknown"),
  levels = c("Independent", "Partially dependent", "Totally dependent", "Unknown"))]
A[, asa_class := factor(fcase(ASACLAS == "1-No Disturb", "ASA 1", ASACLAS == "2-Mild Disturb", "ASA 2",
  ASACLAS == "3-Severe Disturb", "ASA 3", ASACLAS == "4-Life Threat", "ASA 4",
  default = "Not assigned"), levels = c("ASA 1", "ASA 2", "ASA 3", "ASA 4", "Not assigned"))]
A[, inpatient := as.integer(INOUT == "Inpatient")]
A[, admission_origin := factor(fcase(
  TRANST %in% c("Not transferred (admitted from home)", "Home/Permanent residence"), "Home",
  is.na(TRANST) | TRANST == "Unknown", "Unknown", default = "Not home"),
  levels = c("Home", "Not home", "Unknown"))]
A[, diabetes := factor(fcase(DIABETES == "INSULIN", "Insulin", DIABETES == "NON-INSULIN", "Non-insulin",
  default = "None"), levels = c("None", "Non-insulin", "Insulin"))]
yn <- function(x) as.integer(x %in% "Yes")
A[, `:=`(smoker = yn(SMOKE), copd = yn(HXCOPD), heart_failure = yn(HXCHF), hypertension = yn(HYPERMED),
         dialysis_preop = yn(DIALYSIS), disseminated_cancer = yn(DISCANCR),
         steroid_immunosuppression = yn(STEROID), bleeding_disorder = yn(BLEEDDIS),
         preop_transfusion = yn(TRANSFUS), preop_sirs = as.integer(PRSEPIS %in% "SIRS"))]
A[, `:=`(lab_sodium = num_or_na(PRSODM, 100, 170), lab_bun = num_or_na(PRBUN, 1, 150),
         lab_creatinine = num_or_na(PRCREAT, 0.2, 15), lab_albumin = num_or_na(PRALBUM, 1, 6),
         lab_hematocrit = num_or_na(PRHCT, 10, 65), lab_platelets = num_or_na(PRPLATE, 5, 1500),
         lab_wbc = num_or_na(PRWBC, 0.5, 100), lab_inr = num_or_na(PRINR, 0.5, 10))]
LABS <- c("lab_sodium", "lab_bun", "lab_creatinine", "lab_albumin", "lab_hematocrit",
          "lab_platelets", "lab_wbc", "lab_inr")
for (l in LABS) set(A, j = sub("^lab_", "miss_", l), value = as.integer(is.na(A[[l]])))
MISS <- sub("^lab_", "miss_", LABS)
# sensitivity-only legacy fields (collected to 2020 only; NA from 2021)
A[, `:=`(so_weight_loss = ifelse(is.na(WTLOSS), NA_integer_, yn(WTLOSS)),
         so_open_wound = ifelse(is.na(WNDINF), NA_integer_, yn(WNDINF)),
         so_dyspnea = factor(DYSPNEA), so_wound_class = factor(WNDCLAS))]
# surgical plan from CPT (assumed planned, as in the template study)
hasc <- function(codes) {
  m <- as.matrix(A[, ..slots_cols]); Reduce(`|`, lapply(seq_along(slots_cols), function(j) m[, j] %in% codes))
}
A[, instrumentation := factor(fcase(hasc("22844"), "22844 (13 or more)", hasc("22843"), "22843 (7-12)",
  hasc("22842"), "22842 (3-6)", default = "None recorded"),
  levels = c("None recorded", "22842 (3-6)", "22843 (7-12)", "22844 (13 or more)"))]
A[, post_deformity_arthrodesis := factor(fcase(hasc("22804"), "22804 (13 or more)", hasc("22802"), "22802 (7-12)",
  hasc("22800"), "22800 (up to 6)", default = "None"),
  levels = c("None", "22800 (up to 6)", "22802 (7-12)", "22804 (13 or more)"))]
A[, ant_deformity_arthrodesis := factor(fcase(hasc("22812"), "22812 (8 or more)", hasc("22810"), "22810 (4-7)",
  hasc("22808"), "22808 (2-3)", default = "None"),
  levels = c("None", "22808 (2-3)", "22810 (4-7)", "22812 (8 or more)"))]
A[, three_column_osteotomy := as.integer(hasc(c("22206", "22207")))]
A[, posterior_column_osteotomy := as.integer(hasc(c("22212", "22214")))]
A[, pelvic_fixation := as.integer(hasc("22848"))]
A[, interbody_fusion := as.integer(hasc(c("22558", "22630", "22633")))]
A[, combined_ap := as.integer(hasc(ANT_CODES) & hasc(POST_CODES))]
A[, n_cpt_codes := as.integer(n_cpt)]
PRED_CAT <- c("race_group", "hispanic", "functional_status", "asa_class", "admission_origin", "diabetes",
              "instrumentation", "post_deformity_arthrodesis", "ant_deformity_arthrodesis")
PRED_BIN <- c("sex_male", "inpatient", "smoker", "copd", "heart_failure", "hypertension", "dialysis_preop",
              "disseminated_cancer", "steroid_immunosuppression", "bleeding_disorder",
              "preop_transfusion", "preop_sirs", "three_column_osteotomy", "posterior_column_osteotomy",
              "pelvic_fixation", "interbody_fusion", "combined_ap")
PRED_NUM <- c("age", "bmi", LABS, "n_cpt_codes")
PRED_COMMON <- c(PRED_NUM, PRED_BIN, PRED_CAT, MISS)
PRED_SENS_ONLY <- c("so_weight_loss", "so_open_wound", "so_dyspnea", "so_wound_class")

# ---------------------------------------------------------------- outcomes
A[, still := (STILLINHOSP %in% "Yes") | (!is.na(DOPTODIS) & DOPTODIS > 30)]
A[, expired := DISCHDEST %in% "Expired"]
NONHOME <- c("Skilled Care, Not Home", "Unskilled Facility Not Home", "Rehab", "Separate Acute Care",
             "Multi-level Senior Community", "Hospice", "Other Facility", "Acute Care Hospital")
HOME <- c("Home", "Facility Which was Home", "Home/Permanent Residence")
A[, los_days := fifelse(still, 31, as.numeric(DOPTODIS))]
# non-home discharge among patients alive at discharge and discharged within 30 days
A[, elig_nonhome := !expired & !still & DISCHDEST %in% c(NONHOME, HOME)]
A[, y_nonhome := fifelse(elig_nonhome, as.integer(DISCHDEST %in% NONHOME), NA_integer_)]
# LOS: alive at discharge, a known stay length
A[, elig_los := !expired & !is.na(los_days)]
# readmission denominator
A[, elig_readm := !expired & !still & !is.na(DISCHDEST) & !is.na(READMISSION1)]
A[, readm_in30 := READMISSION1 %in% "Yes" & (is.na(READMPODAYS1) | READMPODAYS1 <= 30)]
A[, y_readm_any := fifelse(elig_readm, as.integer(readm_in30), NA_integer_)]
A[, y_readm_unplanned := fifelse(elig_readm, as.integer(readm_in30 & UNPLANNEDREADMISSION1 %in% "Yes"), NA_integer_)]
# SSI
A[, ssi_known := !is.na(SUPINFEC) & !is.na(WNDINFD) & !is.na(ORGSPCSSI)]
A[, ssi_patos := SSSIPATOS %in% "Yes" | DSSIPATOS %in% "Yes" | OSSIPATOS %in% "Yes"]
A[, open_wound := WNDINF %in% "Yes"]
A[, elig_ssi := ssi_known & !ssi_patos & !open_wound]
A[, y_ssi := fifelse(elig_ssi, as.integer(SUPINFEC != "No Complication" | WNDINFD != "No Complication" |
                                           ORGSPCSSI != "No Complication"), NA_integer_)]
A[, y_dehiscence := fifelse(elig_ssi & !is.na(DEHIS), as.integer(DEHIS == "Wound Disruption"), NA_integer_)]
# cardiac
A[, elig_cardiac := !is.na(CDMI) & !is.na(CDARREST)]
A[, y_mi := fifelse(elig_cardiac, as.integer(CDMI == "Myocardial Infarction"), NA_integer_)]
A[, y_cardiac_arrest := fifelse(elig_cardiac, as.integer(CDARREST == "Cardiac Arrest Requiring CPR"), NA_integer_)]
A[, y_cardiac := fifelse(elig_cardiac, as.integer(y_mi == 1L | y_cardiac_arrest == 1L), NA_integer_)]
# renal (legacy era only)
# a legacy-era case is eligible when it is not on dialysis and at least the ARF field is recorded;
# RENAINSF holds NULL for 'No Complication' only if the PUF stores it that way, so test both fields
A[, elig_renal := puf_year <= 2020 & !(DIALYSIS %in% "Yes") & !is.na(OPRENAFL) & !is.na(RENAINSF)]
A[, y_renal := fifelse(elig_renal, as.integer(RENAINSF %in% "Progressive Renal Insufficiency" |
                                                OPRENAFL %in% "Acute Renal Failure"), NA_integer_)]
A[, death30 := as.integer(!is.na(DOPERTOD))]

# LOS cut-off: 75th percentile (observed value, type 1) of the primary cohort, computed once
prim_los <- A[in_primary & elig_los, los_days]
LOS_CUT <- as.numeric(quantile(prim_los, 0.75, type = 1))
A[, y_los_p75 := fifelse(elig_los, as.integer(los_days >= LOS_CUT), NA_integer_)]
A[, y_los_ge9 := fifelse(elig_los, as.integer(los_days >= 9), NA_integer_)]
outcome_csv(data.table(cohort = "primary", n = length(prim_los), p75_days_type1 = LOS_CUT,
  p75_days_type7 = as.numeric(quantile(prim_los, 0.75, type = 7)), median_days = median(prim_los),
  prop_ge_cutoff = mean(prim_los >= LOS_CUT), prop_ge9 = mean(prim_los >= 9),
  note = "still-in-hospital or DOPTODIS>30 coded as 31 days (>30); in-hospital deaths excluded"),
  "cohort_los_cutoff.csv")

# diagnostic: how often outcome-specific field combos are odd
odd <- A[in_primary == TRUE, .(
  readm_unplanned_missing_flag = sum(READMISSION1 %in% "Yes" & is.na(UNPLANNEDREADMISSION1)),
  readm_days_gt30 = sum(READMISSION1 %in% "Yes" & !is.na(READMPODAYS1) & READMPODAYS1 > 30),
  readm_days_null = sum(READMISSION1 %in% "Yes" & is.na(READMPODAYS1)),
  readm_flag_null = sum(is.na(READMISSION1)),
  ssi_fields_null = sum(!ssi_known), renal_fields_null_legacy = sum(puf_year <= 2020 & !elig_renal & !(DIALYSIS %in% "Yes")),
  dischdest_null = sum(is.na(DISCHDEST)), stillinhosp = sum(still))]
outcome_csv(odd, "cohort_outcome_field_checks.csv")

# ---------------------------------------------------------------- outcome-specific flow
OUT_FLOW <- function(d, cname) {
  rows <- list(); add <- function(o, st, desc, rem, remaining)
    rows[[length(rows) + 1]] <<- data.table(section = "outcome", cohort = cname, outcome = o, step = st,
      description = desc, n_removed = rem, n_remaining = remaining)
  seqstep <- function(o, steps) {
    keep <- rep(TRUE, nrow(d)); add(o, 0L, "Cohort", NA_integer_, nrow(d)); i <- 0L
    for (s in steps) { i <- i + 1L; pass <- s[[2]](d); pass[is.na(pass)] <- FALSE
      add(o, i, s[[1]], sum(keep & !pass), sum(keep & pass)); keep <- keep & pass }
  }
  seqstep("nonhome", list(list("Died in hospital (Expired)", function(x) !x$expired),
    list("Still in hospital or stay > 30 days", function(x) !x$still),
    list("Destination missing, AMA or Unknown", function(x) x$DISCHDEST %in% c(NONHOME, HOME))))
  seqstep("los", list(list("Died in hospital (Expired)", function(x) !x$expired),
    list("Stay length missing", function(x) !is.na(x$los_days))))
  seqstep("readmission", list(list("Died in hospital (Expired)", function(x) !x$expired),
    list("Still in hospital or stay > 30 days", function(x) !x$still),
    list("Discharge destination missing", function(x) !is.na(x$DISCHDEST)),
    list("Readmission flag missing", function(x) !is.na(x$READMISSION1))))
  seqstep("ssi", list(list("SSI outcome fields missing", function(x) x$ssi_known),
    list("SSI present at time of surgery (PATOS)", function(x) !x$ssi_patos),
    list("Open wound (field exists 2012-2020 only)", function(x) !x$open_wound)))
  seqstep("cardiac", list(list("MI or arrest field missing", function(x) x$elig_cardiac)))
  seqstep("renal", list(list("2021 (mixed renal definition) excluded", function(x) x$puf_year != 2021),
    list("2022-2024 (new renal definition) excluded", function(x) x$puf_year <= 2020),
    list("Preoperative dialysis", function(x) !(x$DIALYSIS %in% "Yes")),
    list("Renal outcome fields missing", function(x) !is.na(x$OPRENAFL) & !is.na(x$RENAINSF))))
  rbindlist(rows)
}
COHORT_FLAG <- c(primary = "in_primary", broad = "in_broad", strict = "in_strict", primary_nm = "in_primary_nm",
                 dx_required = "in_dx_required")
of <- rbindlist(lapply(names(COHORT_FLAG), function(n) OUT_FLOW(A[get(COHORT_FLAG[n]) == TRUE], n)), fill = TRUE)
flow <- rbindlist(c(flow_rows, list(of)), fill = TRUE)
setcolorder(flow, c("section", "cohort", "outcome", "step", "description", "n_removed", "n_remaining"))
outcome_csv(flow, "cohort_flow.csv")

# ---------------------------------------------------------------- events tables
OUTS <- data.table(
  outcome = c("nonhome", "los_p75", "los_ge9", "readm_unplanned", "readm_any", "ssi", "dehiscence",
              "cardiac", "mi", "cardiac_arrest", "renal"),
  ycol = c("y_nonhome", "y_los_p75", "y_los_ge9", "y_readm_unplanned", "y_readm_any", "y_ssi",
           "y_dehiscence", "y_cardiac", "y_mi", "y_cardiac_arrest", "y_renal"),
  designated = c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE))
ev_long <- function() rbindlist(lapply(names(COHORT_FLAG), function(cn) {
  d <- A[get(COHORT_FLAG[cn]) == TRUE]
  rbindlist(lapply(seq_len(nrow(OUTS)), function(i) {
    y <- d[[OUTS$ycol[i]]]; ok <- !is.na(y)
    data.table(cohort = cn, outcome = OUTS$outcome[i], designated_outcome = OUTS$designated[i],
               puf_year = d$puf_year[ok], era = as.character(d$era[ok]), y = y[ok])
  }))
}))
EV <- ev_long()
summ <- function(by) EV[, .(n = .N, events = sum(y), prevalence = round(mean(y), 4)), by = by]
byo <- summ(c("cohort", "outcome", "designated_outcome"))
outcome_csv(byo[order(cohort, outcome)], "events_by_outcome.csv")
outcome_csv(summ(c("cohort", "outcome", "designated_outcome", "puf_year"))[order(cohort, outcome, puf_year)], "events_by_year.csv")
outcome_csv(summ(c("cohort", "outcome", "designated_outcome", "era"))[order(cohort, outcome, era)], "events_by_era.csv")

# cohort size by year / era and characteristics
cy <- rbindlist(lapply(names(COHORT_FLAG), function(cn) A[get(COHORT_FLAG[cn]) == TRUE, .(cohort = cn, n = .N), by = puf_year]))
outcome_csv(cy[order(cohort, puf_year)], "cohort_by_year.csv")
ce <- rbindlist(lapply(names(COHORT_FLAG), function(cn) A[get(COHORT_FLAG[cn]) == TRUE, .(cohort = cn, n = .N), by = era]))
outcome_csv(ce[order(cohort, era)], "cohort_by_era.csv")

# ---------------------------------------------------------------- missingness by year (primary cohort)
Pm <- A[in_primary == TRUE]
raw_na <- function(col) Pm[, .(variable = col, kind = "cleaned value missing", n = .N,
                               pct_missing = round(100 * mean(is.na(get(col))), 2)), by = puf_year]
miss_vars <- c("age", "bmi", LABS, "so_weight_loss", "so_open_wound")
mm <- rbindlist(lapply(miss_vars, raw_na))
unk <- function(col, lv) Pm[, .(variable = col, kind = "category Unknown/Not assigned", n = .N,
                               pct_missing = round(100 * mean(get(col) %in% lv), 2)), by = puf_year]
mm <- rbind(mm, unk("race_group", "Unknown"), unk("hispanic", "Unknown"), unk("functional_status", "Unknown"),
            unk("asa_class", "Not assigned"), unk("admission_origin", "Unknown"))
so2 <- Pm[, .(variable = "so_dyspnea", kind = "cleaned value missing", n = .N,
               pct_missing = round(100 * mean(is.na(so_dyspnea)), 2)), by = puf_year]
mm <- rbind(mm, so2)
mm[, set := fifelse(variable %in% c("so_weight_loss", "so_open_wound", "so_dyspnea"), "sensitivity_only", "common")]
outcome_csv(mm[order(variable, puf_year)], "cohort_missingness_by_year.csv")

# ---------------------------------------------------------------- candidate parameter count + pmsampsize
mm_pred <- function(vars, d) {
  f <- as.formula(paste("~", paste(vars, collapse = "+")))
  dd <- copy(d[, ..vars]); for (v in vars) if (is.numeric(dd[[v]])) dd[[v]][is.na(dd[[v]])] <- 0
  ncol(model.matrix(f, data = dd)) - 1L
}
P_FULL <- mm_pred(PRED_COMMON, Pm)
P_NOLAB <- mm_pred(setdiff(PRED_COMMON, c(LABS, MISS)), Pm)
P_PARS <- 10L

r2_for <- function(prev, frac) {
  maxr2 <- 1 - exp(2 * (prev * log(prev) + (1 - prev) * log(1 - prev)))
  list(max = maxr2, cs = frac * maxr2)
}
ps_run <- function(prev, p, frac) {
  r2 <- r2_for(prev, frac)
  out <- tryCatch(pmsampsize(type = "b", csrsquared = r2$cs, parameters = p, prevalence = prev, shrinkage = 0.9),
                  error = function(e) NULL)
  if (is.null(out)) return(list(n = NA_real_, events = NA_real_, epp = NA_real_, cs = r2$cs, max = r2$max))
  list(n = out$sample_size, events = out$events, epp = out$EPP, cs = r2$cs, max = r2$max)
}
R2_FRAC <- 0.15
prim <- byo[cohort == "primary"]
P_SETS <- c(full = P_FULL, no_labs = P_NOLAB, parsimonious = P_PARS)
make_ss <- function(cname) {
  pr <- byo[cohort == cname]
  rbindlist(lapply(seq_len(nrow(pr)), function(i) {
    r <- pr[i]; prev <- r$events / r$n
    if (r$events == 0 || prev >= 1) return(NULL)
    rbindlist(lapply(names(P_SETS), function(nm) {
      p <- P_SETS[[nm]]
      g <- ps_run(prev, p, R2_FRAC)
      data.table(cohort = cname, outcome = r$outcome, designated_outcome = r$designated_outcome, predictor_set = nm,
                 n_cohort = r$n, events = r$events, prevalence = round(prev, 4), candidate_parameters = p,
                 max_r2cs = round(g$max, 4), assumed_r2cs = round(g$cs, 4), r2cs_assumption = "15% of max Cox-Snell R2 (Riley 2019 default when no prior model), shrinkage 0.9",
                 required_n = g$n, required_events = g$events, required_epp = g$epp,
                 events_per_parameter = round(r$events / p, 2),
                 riley_met = !is.na(g$n) & r$n >= g$n)
    }))
  }))
}
ss_all <- rbindlist(lapply(c("primary", "dx_required"), make_ss))
ss <- ss_all[cohort == "primary"]
# robustness of assumption: 10% and 30% of max for the full set (required N only)
rob <- rbindlist(lapply(c(0.10, 0.30), function(fr) rbindlist(lapply(seq_len(nrow(prim)), function(i) {
  r <- prim[i]; prev <- r$events / r$n; if (r$events == 0) return(NULL)
  g <- ps_run(prev, P_FULL, fr)
  data.table(outcome = r$outcome, predictor_set = "full", r2_fraction_of_max = fr, required_n = g$n,
             events_per_parameter = round(r$events / P_FULL, 2), riley_met = !is.na(g$n) & r$n >= g$n)
}))))
outcome_csv(ss_all, "sample_size.csv")
outcome_csv(rob, "cohort_sample_size_robustness.csv")

# ---------------------------------------------------------------- tiers (BEFORE any model)
prov <- data.table(outcome = c("nonhome", "los_p75", "readm_unplanned", "ssi", "cardiac", "renal"),
                   provisional_tier = c(1, 1, 1, 2, 3, 3),
                   ycol = c("y_nonhome", "y_los_p75", "y_readm_unplanned", "y_ssi", "y_cardiac", "y_renal"))
make_tiers <- function(cname, ssx) {
  prm <- byo[cohort == cname]
  tt <- rbindlist(lapply(seq_len(nrow(prov)), function(i) {
  o <- prov$outcome[i]; r <- prm[outcome == o]
  s_full <- ssx[outcome == o & predictor_set == "full"]; s_par <- ssx[outcome == o & predictor_set == "parsimonious"]
  ev <- r$events; pt <- prov$provisional_tier[i]
  tier <- NA_character_; rule <- NA_character_; models <- NA_character_
  if (ev < 100) { tier <- "Descriptive only"; rule <- "<100 events: no model developed (Phase 0 stop rule)"; models <- "none; incidence and univariable description" }
  else if (pt == 1 && ev >= 500) { tier <- "Tier 1"; rule <- ">=500 events and provisional Tier 1"; models <- "LR, spline LR, CART, RF, XGBoost, MLP" }
  else if (pt == 1) { tier <- "Tier 2"; rule <- "provisional Tier 1 but <500 events: demoted one tier"; models <- "penalised LR vs XGBoost" }
  else if (pt == 2) { tier <- "Tier 2"; rule <- "provisional Tier 2 and >=100 events"; models <- "penalised LR vs XGBoost, restricted predictors, exploratory" }
  else if (pt == 3 && isTRUE(s_par$riley_met)) { tier <- "Tier 3"; rule <- "provisional Tier 3, >=100 events and Riley criteria met for 10 parameters"; models <- "prespecified penalised LR, <=10 parameters" }
  else { tier <- "Descriptive only"; rule <- "provisional Tier 3 and Riley criteria not met for 10 parameters"; models <- "none; incidence and univariable description" }
  data.table(outcome = o, cohort = cname, n = r$n, events = ev, prevalence = r$prevalence,
             provisional_tier_background = pt, candidate_parameters_full = P_FULL,
             events_per_parameter_full = round(ev / P_FULL, 2), riley_met_full = s_full$riley_met,
             riley_met_10_parameters = s_par$riley_met, tier = tier, models_permitted = models, rule_applied = rule)
  }))
  g9 <- byo[cohort == cname & outcome == "los_ge9"]
  rbind(tt, data.table(outcome = "los_ge9", cohort = cname, n = g9$n, events = g9$events, prevalence = g9$prevalence,
    provisional_tier_background = 1, candidate_parameters_full = P_FULL,
    events_per_parameter_full = round(g9$events / P_FULL, 2),
    riley_met_full = NA, riley_met_10_parameters = NA, tier = tt[outcome == "los_p75", tier],
    models_permitted = "comparator definition only; same tier as los_p75", rule_applied = "comparator of los_p75"), fill = TRUE)
}
# primary rows first (unchanged); dx_required rows are the D1 prespecified sensitivity cohort
tiers <- rbind(make_tiers("primary", ss_all[cohort == "primary"]), make_tiers("dx_required", ss_all[cohort == "dx_required"]))
outcome_csv(tiers, "tiers.csv")
tiers_p <- tiers[cohort == "primary"]

# ---------------------------------------------------------------- predictor dictionary + elective sensitivity
lab_lookup <- c(age = "Age (years)", bmi = "Body mass index (kg/m2)", lab_sodium = "Serum sodium (mEq/L)",
  lab_bun = "Blood urea nitrogen (mg/dL)", lab_creatinine = "Serum creatinine (mg/dL)", lab_albumin = "Serum albumin (g/dL)",
  lab_hematocrit = "Hematocrit (%)", lab_platelets = "Platelet count (x10^3/uL)", lab_wbc = "White cell count (x10^3/uL)",
  lab_inr = "INR", n_cpt_codes = "Number of CPT codes recorded", sex_male = "Male sex", inpatient = "Inpatient status",
  smoker = "Current smoker", copd = "COPD", heart_failure = "Heart failure", hypertension = "Hypertension requiring medication",
  dialysis_preop = "Preoperative dialysis", disseminated_cancer = "Disseminated cancer", steroid_immunosuppression = "Chronic steroid or immunosuppressant use",
  bleeding_disorder = "Bleeding disorder", preop_transfusion = "Preoperative transfusion", preop_sirs = "Preoperative SIRS",
  three_column_osteotomy = "Three-column osteotomy (22206/22207)", posterior_column_osteotomy = "Posterior-column osteotomy (22212/22214)",
  pelvic_fixation = "Pelvic fixation (22848)", interbody_fusion = "Interbody fusion (22558/22630/22633)", combined_ap = "Combined anterior-posterior approach",
  race_group = "Race", hispanic = "Hispanic ethnicity", functional_status = "Functional status", asa_class = "ASA class",
  admission_origin = "Admission origin", diabetes = "Diabetes", instrumentation = "Instrumented segments (CPT)",
  post_deformity_arthrodesis = "Posterior deformity arthrodesis (CPT)", ant_deformity_arthrodesis = "Anterior deformity arthrodesis (CPT)",
  so_weight_loss = "Weight loss >10% (to 2020 only)", so_open_wound = "Open wound (to 2020 only)",
  so_dyspnea = "Dyspnea (to 2020 only)", so_wound_class = "Wound class (to 2020 only)")
lab_lookup[MISS] <- paste0("Missing: ", sub("^lab_", "", LABS))
vtype <- function(v) if (v %in% PRED_CAT || v %in% c("so_dyspnea", "so_wound_class")) "categorical" else if (v %in% c(PRED_BIN, MISS, "so_weight_loss", "so_open_wound")) "binary" else "numeric"
dict <- data.table(variable = c(PRED_COMMON, PRED_SENS_ONLY))
dict[, `:=`(label = lab_lookup[variable], type = vapply(variable, vtype, ""),
            set = fifelse(variable %in% PRED_SENS_ONLY, "sensitivity_only", "common"),
            drift_caveat = fifelse(variable %in% c("race_group", "admission_origin", "heart_failure", "steroid_immunosuppression", "preop_transfusion"), "yes", "no"),
            cpt_derived = fifelse(variable %in% c("instrumentation", "post_deformity_arthrodesis", "ant_deformity_arthrodesis", "three_column_osteotomy", "posterior_column_osteotomy", "pelvic_fixation", "interbody_fusion", "combined_ap", "n_cpt_codes"), "yes", "no"))]
outcome_csv(dict, "cohort_predictors.csv")
base_no_elective <- copy(D)
keep0 <- Reduce(`&`, lapply(base_steps[c(1:3)], function(s) { r <- s[[2]](base_no_elective); r[is.na(r)] <- FALSE; r }))
keep_rest <- Reduce(`&`, lapply(base_steps[5:11], function(s) { r <- s[[2]](base_no_elective); r[is.na(r)] <- FALSE; r }))
elig_p <- elig$primary[[2]](base_no_elective); elig_p[is.na(elig_p)] <- FALSE
neu <- !base_no_elective$dx_neuromusc
core <- keep0 & keep_rest & elig_p & neu
outcome_csv(data.table(
  definition = c("elective (primary: ELECTSURG=Yes & not emergent, or CASETYPE=Elective)", "non-emergent (EMERGNCY=No, or CASETYPE not Emergent)", "difference"),
  primary_cohort_n = c(sum(core & base_no_elective$elective), sum(core & base_no_elective$non_emergent),
                       sum(core & base_no_elective$non_emergent) - sum(core & base_no_elective$elective))),
  "cohort_elective_sensitivity.csv")
# primary-cohort arm composition and ICD split
arm <- Pm[, .(n = .N), by = .(arm = fifelse(c_def_arth, "deformity arthrodesis code (diagnosis not required)", "22843/22844 plus deformity diagnosis"),
                              has_deformity_dx = dx_def, icd_version = icd)]
outcome_csv(arm[order(arm, has_deformity_dx, icd_version)], "cohort_primary_arms.csv")
dx_class <- function(d) fcase(
  d$dx_def, "Deformity (M40/M41/M43.8/Q76; 737)",
  d$icd == 10L & rx("^M48", d$dx_nd), "Spinal stenosis (M48)",
  d$icd == 10L & rx("^M51", d$dx_nd), "Disc disorder (M51)",
  d$icd == 10L & rx("^M96", d$dx_nd), "Post-procedural (M96)",
  d$icd == 10L & rx("^T84", d$dx_nd), "Hardware complication (T84)",
  d$icd == 9L & rx("^996", d$dx), "Hardware complication (996)",
  d$icd == 9L & rx("^(721|722|724)", d$dx), "ICD-9 degenerative (721/722/724)",
  d$icd == 10L & rx("^M43", d$dx_nd), "Other M43",
  default = "Other")
CLASS_LEVELS <- c("Deformity (M40/M41/M43.8/Q76; 737)", "Spinal stenosis (M48)", "Disc disorder (M51)", "Other M43",
                  "Post-procedural (M96)", "Hardware complication (T84)", "Hardware complication (996)",
                  "ICD-9 degenerative (721/722/724)", "Other")
cm <- Pm[, .(arm = fifelse(c_def_arth, "Arthrodesis code 22800-22812", "22843/22844 plus deformity diagnosis"),
             dx_class = factor(dx_class(.SD), levels = CLASS_LEVELS))]
cm_t <- cm[, .N, by = .(arm, dx_class)]
cm_t <- cm_t[CJ(arm = unique(cm$arm), dx_class = factor(CLASS_LEVELS, levels = CLASS_LEVELS)), on = c("arm", "dx_class")]
cm_t[is.na(N), N := 0L]
setnames(cm_t, "N", "n")
cm_t[, pct_of_arm := round(100 * n / sum(n), 2), by = arm]
cm_t[, pct_of_primary := round(100 * n / nrow(Pm), 2)]
cm_t[, dx_class := as.character(dx_class)]
stopifnot(sum(cm_t$n) == nrow(Pm))
outcome_csv(cm_t[order(arm, match(dx_class, CLASS_LEVELS))], "cohort_case_mix.csv")
outcome_csv(Pm[, .(n = .N, n_outpatient = sum(inpatient == 0), pct_outpatient = round(100 * mean(inpatient == 0), 2))], "cohort_inout.csv")

# ---------------------------------------------------------------- save derived (patient level, never committed)
keep_cols <- c("case_uid", "puf_year", "era", "in_primary", "in_broad", "in_strict", "in_primary_nm", "in_dx_required",
               PRED_COMMON, PRED_SENS_ONLY,
               "y_nonhome", "y_los_p75", "y_los_ge9", "y_readm_unplanned", "y_readm_any", "y_ssi", "y_dehiscence",
               "y_cardiac", "y_mi", "y_cardiac_arrest", "y_renal", "los_days", "death30", "dx_def", "dx_neuromusc", "icd")
AR <- A[, ..keep_cols]
setorder(AR, puf_year, case_uid)
write_parquet(AR, file.path(P$derived, "analysis_ready.parquet"), compression = "uncompressed")
saveRDS(list(common = PRED_COMMON, numeric = PRED_NUM, binary = PRED_BIN, categorical = PRED_CAT,
             missing_indicators = MISS, labs = LABS, sensitivity_only = PRED_SENS_ONLY,
             outcomes = OUTS$ycol, designated_outcomes = prov$ycol, los_cutoff_days = LOS_CUT,
             candidate_parameters = c(full = P_FULL, no_labs = P_NOLAB, parsimonious = P_PARS),
             cohort_flags = COHORT_FLAG, seed = SEED), file.path(P$derived, "predictor_sets.rds"), compress = FALSE)

# ---------------------------------------------------------------- summary
cat("Primary n =", sum(A$in_primary), "\n")
print(byo[cohort == "primary"][order(outcome)])
print(tiers[, .(cohort, outcome, n, events, prevalence, events_per_parameter_full, tier)])
cat("LOS cutoff:", LOS_CUT, " P_FULL:", P_FULL, " P_NOLAB:", P_NOLAB, "\n")
fs <- c(list.files(P$results, pattern = "^(cohort_.*\\.csv|events_.*\\.csv|tiers\\.csv|sample_size\\.csv)$", full.names = TRUE),
        file.path(P$derived, c("analysis_ready.parquet", "predictor_sets.rds")))
fs <- fs[basename(fs) != "cohort_checksums.csv"]
chk <- tools::md5sum(fs)
fwrite(data.table(file = basename(fs), md5 = unname(chk))[order(file)], file.path(P$results, "cohort_checksums.csv"))
