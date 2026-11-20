# Cohort summary (T1, Phase 0 gate)

Source: ACS-NSQIP Adult PUF 2012-2024 via `nsqip.duckdb`; script `analysis/01_build_cohort.R`; all numbers in `results/*.csv`.

## Cohort
- Candidate pool 20,215 -> primary **7,699** (elective, ASA < 5, no preop sepsis/vent, no tumour/infection/fracture principal dx, not cervical-only, deformity arthrodesis code or 22843/22844 plus deformity dx, M41.4x excluded). By year 360 (2012) to 709 (2019); era blocks 2,107 / 1,736 / 2,012 / 1,844.
- Sensitivity cohorts: broad 13,730; strict 4,058; primary + neuromuscular 7,799; **diagnosis required in every arm (dx_required) 4,604** (prespecified, D1).
- 3,095 primary cases (40.2%) qualify on the arthrodesis code without a deformity diagnosis; flagged author decision (promote dx_required to primary?). Case mix in `cohort_case_mix.csv`.
- dx_required events (n): non-home 1,686 (4,542); LOS >= 8 d 1,452 (4,592); readmission 345 (4,542); SSI 187 (4,570); cardiac 68 (4,604); renal 19 (3,093). Tiers there: non-home and LOS Tier 1, readmission and SSI Tier 2, cardiac and renal Descriptive only.

## Events, primary cohort (denominator, events, prevalence)
| Outcome | n | Events | % | Tier |
|---|---|---|---|---|
| Non-home discharge | 7,602 | 2,448 | 32.2 | 1 |
| Prolonged LOS (>= 8 d, cohort 75th percentile) | 7,671 | 1,975 | 25.7 | 1 |
| LOS >= 9 d (comparator) | 7,671 | 1,440 | 18.8 | with LOS |
| Unplanned readmission (30 d) | 7,597 | 557 | 7.3 | 1 (>= 500 events) |
| Any SSI | 7,639 | 273 | 3.6 | 2 |
| Cardiac (MI or arrest) | 7,699 | 107 | 1.4 | 3 |
| Renal (2012-2020 only) | 5,187 | 30 | 0.6 | Descriptive only |

Dehiscence 67, MI 71, arrest 45 (descriptive). Era-block events: readmission 136/134/148/139; SSI 61/56/76/80 (rising); cardiac 24/23/32/28.

## Sample size (pmsampsize, 62 candidate parameters, Cox-Snell R2 = 15% of max)
Readmission is Tier 1 by the >= 500 event rule but the full-set Riley criteria are NOT met (events per parameter 9.0 vs 10 needed; required N 8,804 vs 7,597). Non-home and LOS meet them. SSI (4.4 events per parameter) and cardiac (1.7) do not for the full set; both meet them for a 10-parameter model. Renal fails everything.

## Caveats
- LOS 75th percentile is 8 days (published 9 is the comparator).
- Readmission flag covers readmission 1 only; counted from surgery.
- Open-wound SSI exclusion impossible from 2021; renal definition changes July 2021 (legacy era only).
- Elective era-equivalence unconfirmed (+799 cases under a non-emergent rule).
- CPT-derived plan variables are assumed planned and come from performed codes; T4 must run drop-`n_cpt_codes` and drop-all-CPT-derived sensitivities.
- Pelvic/ASD-case mix not checked beyond CPT flags; strict cohort tests it.
