# Post hoc extension E10, 2026-10-10: stability of the 2026 dyslipidemia guideline's 30-year-risk statin pathway on remeasurement, and its contemporary US reach

Written before any E10 quantity was computed. Only data availability was checked:
- TRIGLY_L (NHANES August 2021–August 2023) has the fasting weight WTSAF2YR and LDL-C by NIH equation 2 (LBDLDLN; non-missing count only).
- The Study 1 30-year risks exist in the E9 code path.

E9's rule-v2 output (`ext3_statin_results.json`) has been seen, so this extension is post hoc and descriptive; it is not preregistered. No frozen estimate changes. Results go to a results file and the evidence ledger; the manuscript text stays frozen.

## Motivation

- The 2026 guideline's COR 2a pathway is: age 30–59, PREVENT 10-year ASCVD <3%, LDL-C 70–159 mg/dL and base PREVENT 30-year ASCVD ≥10%. It accounts for most newly eligible US adults (Essa 2026; Spitz 2026; see `literature_check_20261010.md`).
- Both eligibility studies are cross-sectional and cannot assess within-person change.
- In E9 (rule v2), changes in the COR 1-or-2a reading were dominated by 30-year-risk and LDL-C-160 crossings.

## Part A: Study 1 (NHANES III), remeasurement stability

**Population A1.** Study 1 analysis set (759), aged 30–59 at visit 1, not taking a statin (statin = 0). Counts are reported.

**30-year risk.** Base PREVENT 30-year ASCVD at each visit, computed as in E9:
- the same clamped visit inputs and fixed covariates as the 10-year risk, with age at visit 1;
- the pinned preventr `base_30yr` table, checked against pinned pyprevent (≤1e-9 percentage points; non-finite values fail the check).

**Primary estimand A (θ30).** Proportion of A1 whose 30-year ≥10% classification differs between visits. Interval: CR1 logit-t over the 98 pseudo-PSUs (Korn-Graubard if fewer than 10 events).

**Threshold-distance curve π30.** Same specification as the post hoc 7.5% total-CVD curve (`nh3_posthoc_cvd75_curve.py`):
- logistic regression of the flip indicator on a restricted cubic spline (three knots at the 10th, 50th and 90th percentiles) of |logit(r30 at visit 1 / 100) − logit(0.10)|;
- the same 4,000 cluster-bootstrap replicates (S.boot_weights, seed 20261009), with pointwise 95% percentile intervals;
- π30 at the threshold;
- instability zones: the contiguous 30-year-risk range around 10% where π30 ≥25% and where π30 ≥10%;
- a lookup at 30-year risks of 5, 7, 8, 9, 10, 11, 12, 13, 15 and 20%.

**Secondary A.**
1. **LDL-C 160 mg/dL** (Sampson, as E9), in A1 with triglycerides ≤800 at both visits:
   - proportion whose LDL-C ≥160 classification differs between visits;
   - the same curve form on |ln(LDL-C at visit 1) − ln(160)|, with zones in mg/dL;
   - if there are fewer than 15 crossings, only the proportion is reported.
2. **Pathway domain at visit 1**: A1 with triglycerides ≤800 at both visits; no risk-independent indication by rule v2; LDL-C 70–159; 10-year risk <3%. Within it:
   - the proportion whose 30-year ≥10% status differs at visit 2;
   - the proportion whose E9 rule-v2 COR 1-or-2a status (I, R, C, S or H) differs at visit 2;
   - each separately for visit-1 pathway-eligible (30-year ≥10%) and not eligible.

   These condition on visit 1 and so include regression to the mean.
3. **Single visit versus two-visit mean**: the 30-year risk from two-visit mean inputs, defined as in E9. Report the proportion whose single-visit 30-year ≥10% classification differs from the mean-based one, averaged over the two visits, with a paired interval (S.paired_ci).
4. **Clustering**: θ30 with visit-1 physician clustering (CR1, plus missing code 88888 as singletons).

## Part B: NHANES August 2021–August 2023, exploratory transport

**Population B.** As E2:
- examined; aged 30–59 (the 30-year model's age range); explicitly no CHF, CHD, angina, MI or stroke; not pregnant;
- in the fasting subsample with WTSAF2YR >0;
- not taking cholesterol medication (BPQ101D ≠ 1; statin proxy as in E2; missing excluded);
- complete base PREVENT inputs and LDL-C (LBDLDLN).

Inputs as E2:
- oscillometric SBP: mean of readings 2–3, otherwise of the available readings;
- diabetes by questionnaire or HbA1c ≥6.5%;
- smoking; BP treatment; statin term 0;
- eGFR by 2021 CKD-EPI;
- the same clamping.

**Domain B1 (30-year pathway domain).**
- No rule-v2 risk-independent indication: LDL-C <190; not age 40–75 with diabetes; not age 40–75 with eGFR <60 and LDL-C ≥70.
- LDL-C 70–159.
- PREVENT 10-year ASCVD <3%.

**Estimands.** Survey-weighted with WTSAF2YR and SDMVSTRA/SDMVPSU (Taylor linearization; the Svy class from E2). Millions are reported as complete-case weighted totals and also scaled as in E2.
1. Share and number of B1 with 30-year risk ≥10% (pathway-eligible). This is an external plausibility check against Essa 2026, not a formal comparison: different years, LDL equation and statin handling.
2. Share and number of B1, and of pathway-eligible B1, whose 30-year risk lies in the Part A π30 ≥25% zone and in the ≥10% zone. Zone limits are fixed at the Part A point estimates.
3. Transported expected proportion of B1 whose 30-year ≥10% classification would differ on remeasurement: π30 averaged over B1. The interval combines the 4,000 refitted curve replicates with the survey variance, as in E2.
4. If secondary A1 yields a curve: the share of B1 with LDL-C in the LDL-C-160 ≥25% zone.

**Assumptions to state.**
- 1988–1994 within-person variation and distance-only dependence apply to 2021–2023.
- Mercury versus oscillometric BP.
- Fasting LDL-C by NIH equation 2.
- Statin proxy from questionnaire.
- Risk enhancers, CAC, HIV and LDL-C goals not modelled.
- This is a projection, not an observed remeasurement.

## Reporting and limits

- Code: `scripts/nh3_ext4_30yr.py`. Output: `analysis_v1_posthoc/ext4_30yr_results.json`.
- Results enter the evidence ledger. Part B is labelled exploratory (PH-P, transport).
- Departures are recorded in an addendum to this file.
- An independent code audit follows the first run.

## Addendum after an independent audit of the first run (2026-10-10; protocol text above unchanged, SHA-256 of the text before this addendum 5375ED06...6F15; file time 20:33:21)

The audit reproduced every first-run number with independent code. The 2021-2023 codebooks confirm the definitions:
- BPQ101D: 1/2 = yes/no; any cholesterol medicine.
- LBDLDLN: NIH equation 2, matching Sampson to rounding.
- WTSAF2YR: the fasting MEC weight.
- MCQ160B-F codes, and the design df of 15.

Wording and reporting defects, fixed in the final run with first-run estimates unchanged (first run kept as ext4_30yr_results_run1.json):

1. "Eligible at visit 1 and not at visit 2" was mislabelled. The 14/70 count is a loss of the 30-year >=10% criterion only. The final run adds:
   - the strict count of visit-1 pathway-eligible participants who were not pathway-eligible at visit 2 (out of the domain, or 30-year <10%), with reasons;
   - for visit-1 not-eligible participants, the count who became pathway-eligible.

   All conditional counts are labelled as conditioning on visit 1, so they include regression to the mean.
2. Zone language: the pi30 zones are where remeasurement changes the 30-year classification, not eligibility as a whole. Eligibility can also change through the 10-year 3% cut, LDL-C 70/160 and the indications.
3. The Part B zone shares have survey-only intervals, because the zone limits are fixed at the Part A point estimates (as E2). The final run adds a version that carries the curve uncertainty: the zone is recomputed in each of the 4,000 curve replicates and the replicate variance is combined with the survey variance.
4. Additions not requested by the protocol:
   - counts and descriptive summaries: ge10 counts, up/down, median/IQR, in-zone counts, tier counts, the LDL-C lookup, B1 share and totals;
   - curve grids and a weighted 30-year-risk histogram, for a figure.

   These are descriptive extras.
5. Population numbers are complete-case totals; the scaled totals are reported beside them.
6. Assumptions list completed: mercury (1988-1994) versus oscillometric BP; LDL-C goals and risk enhancers not modelled.
7. Trivia:
   - CR1 intervals use the clusters present in each subset (97 in A1), not all 98.
   - Zones come from a 0.001-logit grid.
   - The transported-flip interval uses 1.96 (t with 15 df would give 4.80-10.06%).