# Machine learning versus logistic regression for outcomes after elective adult spinal deformity fusion (ACS NSQIP 2012-2024)

Analysis code and aggregate results for a comparison of logistic regression with flexible models
(restricted-cubic-spline logistic regression, penalized regression, CART, random forest, XGBoost,
single-hidden-layer neural network) for non-home discharge, prolonged hospitalization and unplanned
30-day readmission, with surgical site infection, cardiac events and renal failure as lower-event outcomes.
Primary cohort: 7,699 elective deformity fusion operations.

## Contents

| Path | Content |
|---|---|
| `analysis/00_setup.R` | paths and read-only DuckDB helper |
| `analysis/01_build_cohort.R` | cohort, outcomes, predictors, modeling tiers, sample-size criteria |
| `analysis/02_models.R` | repeated stratified nested cross-validation (5 x 5 x 5) and 60/20/20 split fits |
| `analysis/03_metrics.R` | AUROC, AUPRC, Brier score, calibration, paired bootstrap differences versus logistic regression, decision curves |
| `analysis/04_sensitivity.R` | era and forward-chaining validation, cohort and predictor-set sensitivity, imputation, recalibration, instability |
| `lib/` | shared preprocessing, fitting, cross-validation and metric functions |
| `tests/` | test scripts for the cohort, model, metric and sensitivity steps |
| `results/` | aggregate outputs only (csv, md); no patient-level data |
| `refs/` | reference list in JSON and Markdown |
| `run_all.R` | driver for steps 00-04 |

## Data

The ACS NSQIP Adult Participant Use Files are not redistributed. Obtain them from ACS NSQIP under their
data use agreement, load them into DuckDB, and set `NSQIP_DUCKDB` to that file. Patient-level intermediate
files (`derived/`) are written locally by step 01 and are excluded from this repository.

## Run

```
Rscript run_all.R
```

Steps 01-04 are computationally heavy (nested cross-validation and bootstrap refits). Analyses used R 4.6.0
(glmnet, ranger, xgboost, nnet, rms, dcurves) with fixed seeds.

## License

MIT. See `LICENSE`.
