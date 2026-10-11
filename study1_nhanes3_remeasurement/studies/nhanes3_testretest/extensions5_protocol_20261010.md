# Post hoc extension E11, 2026-10-10: clinical stakes of measurement-driven decision changes (net benefit lost by single-visit decisions)

Written before any E11 quantity was computed. Post hoc and model-based; not preregistered. No frozen estimate changes. Results go to a results file and the evidence ledger; the manuscript text stays frozen.

## Question and rival

- **Rival (Rahimi 2025; the expected reviewer objection).** Decisions that flip on remeasurement belong to people whose usual-level risk is close to the threshold, where treatment's expected net benefit is close to zero. So flips may be frequent yet clinically low-stakes.
- **E11 question.** What fraction of the net benefit achievable by deciding on usual-level risk is lost when the same risk-threshold rule is applied to a single visit's measurement? How much do two-day averaging, or a second day only near the threshold, recover?
- **Scope.** Risk-threshold components only. Rule components defined by single measurements (BP stage, LDL-C 160/190, HbA1c 6.5%) are not risk-based. Their stakes are not captured, and this is stated as a limit.

## Framework (decision-curve net benefit; the guideline threshold sets the exchange rate)

For threshold pt and usual-level risk r, treating rather than not treating has net benefit (r − pt)/(1 − pt), in true-positive-equivalent units over the risk horizon.
- **Regret of a single-visit decision:** |r − pt|/(1 − pt) if the single-visit decision differs from the usual-level decision; otherwise 0.
- **Achievable net benefit** (usual-level rule versus treating nobody): E[max(r − pt, 0)]/(1 − pt).
- **Primary metric:** loss fraction = E[regret] / achievable net benefit. It does not depend on any relative risk reduction, but it assumes the equation is calibrated and that the threshold reflects the benefit–harm exchange rate.
- **Secondary metrics:**
  - absolute regret per 1,000 people;
  - mean |r − pt| among people whose single-visit decision differs from the usual-level decision, in percentage points;
  - the single-visit decision disagreement itself.

## Measurement model (as E4)

- **Linear predictor.** Visit v value: L_v = T + s·δ/2 + ε_v, with ε ~ N(0, σ²), s = −1 at visit 1 and +1 at visit 2. δ and σ come from the two-visit differences in the analysed population.
- **Usual level.** The distribution G of T is the NPMLE: 400-point grid, EM run to the 1e-9 tolerance or 100,000 iterations, convergence reported. It is fitted to the two-visit means with SD σ/√2.
- **Expectations.** Computed over G; single-visit quantities are averaged over the two visits (s = ±1).
- **Usual level defined.** The person's average over the about 2-week two-visit process, as in E4. Not a lifetime mean.

## Thresholds and populations (Study 1, NHANES III)

1. **Base PREVENT 10-year ASCVD at 5%** (the 2026 statin COR 1 threshold within the risk pathway), statin-naive participants. **Primary.** Also reported at 3% and 10%, same population.
2. **PREVENT 10-year total CVD at 7.5%** (the 2025 BP drug threshold within stage 1), among participants not taking BP medication. Total-CVD risks come from the frozen Study 1 participant file.
3. **Base PREVENT 30-year ASCVD at 10%** (the 2026 COR 2a pathway), in the E10 Part A population (ages 30–59, statin-naive). Risks computed as in E10.

## Strategies compared (for each threshold)

- (a) A single visit.
- (b) The average of two independent days (error SD σ/√2).
- (c) A targeted second day, only when the first single-visit risk lies in that threshold's π ≥25% instability zone; otherwise decide from the first visit.
  - Zones: 5% from the frozen π_clin (K3), clipped at the category midpoints as in E2; 7.5% from `cvd75_curve.json`; 30-year 10% from E10.
  - For 3% and 10%, the frozen K3 zones.
  - Reported with the expected proportion remeasured.
  - Computed by Monte Carlo over G: 200,000 draws per threshold, fixed seed 20261011.

## Uncertainty

- **Percentile intervals** from the first 1,000 of the frozen 4,000 cluster-bootstrap replicates (S.boot_weights, seed 20261009), restricted to each population. EM is warm-started from the full-sample solution, with at most 400 iterations per replicate (as E4).
- **Monte Carlo strategy (c) in replicates:** 20,000 draws.

## Sensitivity (primary threshold 5% only)

- (i) Empirical symmetrized single-visit errors in place of the normal (as E4).
- (ii) Smooth g-model (gsmooth, as E4) in place of the NPMLE.
- (iii) A model-free empirical version that treats the two-visit mean as the usual level: regret of visit-1 and visit-2 decisions relative to the mean decision. This is expected to understate regret.

## Prespecified interpretation of the primary loss fraction (single visit, 5%)

| Loss fraction | Reading |
|---|---|
| <5% of achievable net benefit | Minor: flips are frequent but low-stakes for net benefit. The rival holds for outcomes. The clinical claim narrows to reproducibility of advice, communication and workload. |
| 5% to <15% | Modest. |
| ≥15% | Substantial: supports remeasurement before decisions near the threshold. |

If the 95% interval spans two readings, the result is reported as such, with no single reading chosen. The same reading scale is applied descriptively to the other thresholds.

## Reporting and limits

- **Code:** `scripts/nh3_ext5_stakes.py`.
- **Output:** `analysis_v1_posthoc/ext5_stakes_results.json`.
- **Limits to state:**
  - assumes calibration and the threshold exchange rate;
  - risk-based components only;
  - the usual level is defined over 2 weeks;
  - 1988–1994 measurement;
  - model-based.
- **Audit:** an independent code audit follows the first run.

## Addendum after an independent audit of the first run (2026-10-10; written after the first run, before the changes below were computed; protocol text above unchanged, SHA-256 of the text before this addendum 9601C8D6...9C70)

The audit reproduced every first-run point estimate exactly with independent code (5% loss fraction 0.8666%; 30-year 1.5650%; 10-year 10% 7.6877%; 3% 0.4356%; 7.5% 0.8854%). It confirmed the formulas, the averaging over visits, the two-day and targeted strategies, and the zone sources. Defects and post-audit changes:

1. Bootstrap starting point. Each replicate's EM was warm-started from the full-sample solution and limited to 400 iterations (as specified here and in E4). The audit found that this narrows the intervals relative to refitting each replicate as the point estimate was fitted (from a uniform start).
   - Final run: also refits every replicate from a uniform start to convergence, using SQUAREM-accelerated EM. SQUAREM is checked against plain EM on the full sample.
   - Both intervals are reported. The uniform-start interval is used for interpretation. The reading rule above is applied to it, and an interval spanning two readings is reported as such.
2. The 30-year bootstrap departed from this protocol: replicates were generated from the subgroup's own cluster labels (97 of 98 clusters present). Final run: the frozen replicates are built on the full cluster vector and restricted to the subgroup's columns. In the other two populations all 98 clusters are present, so their replicates were already the frozen ones.
3. Post hoc sensitivities prompted by the audit (not prespecified):
   - (a) error SD by age band (30-39, 40-59, 60-79; within-person SD falls with age), as E4's age sensitivity;
   - (b) error variance relative to a longer-term usual level inflated 1.5-fold and 2-fold. The extra variance is shared by both visits, so it is not removed by averaging two days;
   - (c) the 5% threshold within ages 30-59;
   - (d) loss relative to the better of treating all and treating none (relevant where treating all has positive net benefit, eg, 3%).
4. Reporting:
   - Regret per 1,000 is reported beside every loss fraction, and the fraction is described as specific to this population.
   - The 7.5% analysis covers all participants not taking BP medication, not only stage 1.
   - "Recovers about half" applies only at 5% and the 30-year threshold.
   - The single-visit disagreement with the usual level is not called kappa.
   - The 10-year 10% result is described as depending on the error model, not as a single reading.