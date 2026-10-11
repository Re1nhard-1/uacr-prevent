# Post hoc extensions, set 2, 2026-10-10 (written before running E3; see disclosure for E4-E8)

All analyses are post hoc and descriptive, specified after every Study 1 and Study 2 result was known, and not preregistered. No frozen estimate changes.
They respond to an independent review of revision 5 (clinical, methodological and literature reviewers).

**Disclosure.** Before this protocol was written, a methodological referee ran rough exploratory versions of E4-E8 (usual-level misclassification
about 7%; examiner-clustered interval about 7.4-14.4%; BP minus lipid Shapley difference 1.4 pp, -0.6 to 3.4; heteroscedastic standardization
about 8.4%; Study 2 shift-removed disagreement about 4.5%). E4-E8 below are therefore not blind. E3 (decision level) had not been computed by
anyone; only data availability was counted (DBP readings 2-3 at both visits for all 759; 623 untreated).

## E3. Antihypertensive drug-initiation decision under the 2025 AHA/ACC high blood pressure guideline (primary of this set)

- **Population.** Study 1 analysis set (759), untreated (bptreat = 0). Treated participants: descriptive secondary only.
- **BP at each visit.** SBP = mean of physician readings 2 and 3 (unclamped; as in Study 1). DBP = mean of K5 readings 2 and 3 (PEP6H3, PEP6I3;
  PER6H3, PER6I3); a value of 0 or a missing code is missing; if one of the two is missing, the other is used; if both are missing, the mean of
  the available readings 1-3; if none, the participant is excluded from E3 (count reported).
- **Categories** (inclusive lower bounds): stage 2 if SBP >=140 or DBP >=90; stage 1 if not stage 2 and (SBP >=130 or DBP >=80); otherwise below
  130/80.
- **High-risk flag at each visit**: diabetes (Study 1 visit definition dm_v), or CKD (eGFR_v <60 mL/min/1.73 m2, calibrated creatinine, or UACR_v
  >=30 mg/g where UACR is available), or base PREVENT 10-year total CVD risk >=7.5% (Study 1 visit risk, clamped inputs).
- **Decision at each visit**: D3 stage 2 (drug now; two drugs); D2 stage 1 with high-risk flag (drug now); D1 stage 1 without it (lifestyle first,
  drug if still >=130/80 after 3-6 months); D0 below 130/80 (no drug).
- **Primary estimand**: proportion of untreated participants whose drug-now status (D2 or D3 versus D0 or D1) differs between visits (theta_Rx).
- **Secondary**: 4x4 transition table; flips of the hypertension label (>=130/80) and of stage 2 (two-drug start); conditional flips
  P(not drug-now at visit 2 | drug-now at visit 1) and the converse; attribution with two input groups (BP: SBP and DBP; non-BP: lipids, eGFR,
  diabetes, UACR), as the proportions flipping when only one group is remeasured, averaged over both reference visits, and the two-player
  Shapley split of theta_Rx; flips caused by the 7.5% risk criterion alone (stage 1 at both visits with a change in the high-risk flag).
- **Guideline-concordant comparator**: decision from the average of both visits (SBP and DBP means; geometric means of total and HDL cholesterol,
  creatinine and UACR; diabetes if fixed diabetes or mean HbA1c >=6.5%; PREVENT-CVD from the averaged, clamped inputs). Reported: proportion whose
  single-visit decision differs from the two-visit decision (mean over visit 1 and visit 2).
- **Uncertainty**: as Study 1 (CR1 logit-t over the 98 pseudo-PSUs, Korn-Graubard if <10 events); E5 examiner clustering as a sensitivity.
- **Sensitivity**: (a) means of all three readings for SBP and DBP; (b) CKD by eGFR only; (c) confirmed diabetes fixed across visits (SA3a rule);
  (d) excluding the 54 referred or unknown-referral participants (SA5 rule).
- **Treated (descriptive)**: flips of 'above 130/80' between visits.
- **Interpretation limits**: office readings only; the guideline also asks for readings on two or more occasions and out-of-office confirmation;
  CKD from a single specimen; 1988-1994 mercury readings recorded to even numbers.

## E4. Misclassification relative to the usual level (kappa)

- **Model** (logit of base PREVENT-ASCVD risk, LP): LP_iv = mu_i + s_v delta/2 + e_iv (s_1 = -1, s_2 = +1), e ~ N(0, sigma^2), independent.
  delta = mean(LP2 - LP1); sigma^2 = var(LP2 - LP1)/2. Usual level mu_i = the person's mean over the two-visit process.
- **G** (distribution of mu): nonparametric maximum likelihood from m_i = (LP_i1 + LP_i2)/2 with known error N(0, sigma^2/2), on a 400-point grid,
  EM to convergence.
- **Estimands**: kappa_model = mean over v of P(category(mu + s_v delta/2 + e) != category(mu)) under G; theta_model = P(category at visit 1 !=
  category at visit 2) under the same model; calibrated kappa = kappa_model x theta_obs/theta_model (primary); kappa(r) = P(usual category !=
  category(r) | first-visit risk r), from the posterior of mu given LP_1, at the Table 3 risks; the model-free lower bound theta_obs/2.
- **Sensitivity**: error distribution equal to the symmetrized empirical distribution of (LP2 - LP1)/sqrt(2); error SD by age band.
- **Uncertainty**: 1,000 cluster-bootstrap replicates (frozen weights, seed 20261009), re-estimating delta, sigma, G and theta_obs.

## E5. Examiner clustering

- CR1 logit-t intervals for theta with clusters = visit-1 physician, visit-2 physician and physician pair; wild cluster bootstrap (Webb six-point,
  9,999 replicates, seed 20261010) percentile-t interval for visit-1 physician clustering. Theta by same versus different physician, with the
  number of changes expected from the frozen pi_clin model in each group.

## E6. Blood pressure versus lipids

- Psi_BP - Psi_LIPID (base model) with the frozen 4,000 cluster-bootstrap replicates (percentile 95% CI); the frozen Psi_BP interval must be
  reproduced.

## E7. Heteroscedastic standardization

- Cells: age band (30-39, 40-59, 60-79) x race and ethnicity (non-Hispanic White, non-Hispanic Black, Mexican American; other merged with
  non-Hispanic White). sigma_c = SD(LP2 - LP1) in cell c (sample). Distance d* = d x sigma_all/sigma_c. Refit pi_clin on d* with the frozen
  specification (three-knot spline at the 10/50/90th percentiles of d*, nearest-threshold indicator); apply to the 9,478 targets with their
  cell sigma_c. Joint bootstrap as K4 (2,000 replicates; cluster bootstrap refit with sigma_c re-estimated; Rao-Wu target). Report
  theta_std,het and its paired difference from the frozen theta_std.

## E8. Study 2: systematic shift versus random variation

- Weighted geometric mean ratio of first-morning to random UACR (examination weights). Rescale first-morning UACR by its inverse; recompute
  PREVENT-UACR; proportion changing category (Taylor-logit, 16 df), with down and up counts. Run only if the Study 2 inputs reproduce the locked
  risks exactly; otherwise reported as not run.

## Reporting

E3 in the main text (Results, one paragraph; Discussion); E4 as a column of Table 3 and one Results sentence; E5-E8 in the supplement with one
sentence each where they qualify a main-text claim. Departures are recorded in an addendum to this file.

## Addendum after the first run (2026-10-10; protocol text above unchanged, SHA-256 of the text before this addendum DCCD10E9...B11B)

1. E4 pointwise kappa(r): the NPMLE of G is discrete (atoms), so pointwise kappa(r) near thresholds depended on atom placement (final run:
   NPMLE 62.5%, 66.8% and 39.1% at 3%, 5% and 10%). A smooth g-model (log g polynomial of degree 7, ridge penalty 1) was added as a
   sensitivity after this was seen; it gave 50.2%, 51.2% and 53.6%. Overall kappa was stable across estimators (final run: NPMLE 6.9%, smooth
   7.4%, empirical errors 6.5%, age-specific SD 6.8%). kappa(r) is therefore not shown as a column of Table 3 (departure from the reporting
   plan); both versions are reported in the supplement and the main text gives overall kappa only. kappa(r) values are uncalibrated posteriors.
2. E4 bootstrap used the first 1,000 of the 4,000 frozen cluster-bootstrap replicates, with the NPMLE warm-started from the full-sample
   solution and limited to 400 EM iterations per replicate.
3. E8 ran: the Study 2 inputs reproduced the locked risks (maximum difference 4 x 10^-14 pp) and the locked primary estimate and interval.

## Second addendum after an independent audit of the first run (2026-10-10)

4. E4 grid and convergence: the first run used 200 grid points and at most 3,000 EM iterations (kappa 7.0%); the protocol specifies 400
   points and EM to convergence. The final run uses 400 points and up to 100,000 iterations; EM did not reach the 10^-9 tolerance (final
   maximum change in mass 8.8 x 10^-8), and kappa changed by about 0.1 pp (6.9%, 95% CI 5.5-8.8%). First-run outputs are preserved in
   ext2_results_run1.json.
5. E3 physician clustering, required by the E3 uncertainty clause, was omitted from the first run and added (CR1 over 13 visit-1 physicians and
   a wild cluster bootstrap).
6. E7 first run drew the Rao-Wu target replicates with a different seed from the frozen K4; the final run uses the frozen seed and reproduces
   the frozen theta_std interval exactly.
7. The disclosure above said DBP readings 2-3 were available at both visits for all 759; zero K5 values were counted as available. Three
   untreated participants had K5 = 0 at every reading of one visit and were excluded as the E3 rule specifies (620 analysed).
8. Language: kappa is described as model-based disagreement with the usual-level category, the usual level being the person's average over
   the two-visit (about 2-week) process, not as misclassification (frozen Study 1 protocol language rule).
9. Reporting beyond the plan: E3 also appears in the abstract, Key Points, Central Illustration and Limitations; E4 also in the Discussion;
   E5 in Results and Limitations; E8 in Results and Discussion. The protocol's '7.5% criterion alone' count (stage 1 at both visits with a
   change in the high-risk flag; 5) is reported, with a narrower count (no change in diabetes or kidney status; 1) added after the first run.
10. Timing statements in the manuscript were corrected: only E3 was computed after the protocol was saved without prior rough values.

## Third addendum (2026-10-11; written before computing; protocol text above unchanged)

E4 bootstrap starting values. The E11 audit (extensions5_protocol_20261010.md, addendum) found that warm-starting each bootstrap replicate's EM
from the full-sample NPMLE solution with at most 400 iterations narrows the interval relative to refitting each replicate as the point estimate
was fitted. E4 used the same scheme (first addendum, item 2). Post-audit addition, computed by `scripts/nh3_ext2_e4_uniform.py`:
- the same first 1,000 frozen cluster-bootstrap replicates, each refitted from a uniform start with SQUAREM-accelerated EM (20,000 EM
  evaluations; the full-sample SQUAREM fit is checked against the plain-EM E4 point estimate);
- percentile intervals for kappa, kappa_model, kappa/theta and kappa(r), reported beside the original intervals; the uniform-start interval is
  used for interpretation. The E4 point estimates are unchanged. Output: analysis_v1_posthoc/ext2_e4_uniform_ci.json.