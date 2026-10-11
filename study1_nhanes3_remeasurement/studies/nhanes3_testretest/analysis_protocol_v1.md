# Stability of PREVENT 10-year ASCVD risk categories when every routinely measured input is remeasured: frozen analysis protocol for the NHANES III second-examination test-retest study (v1, 2026-10-09; reviewer-corrected)

Status: FROZEN 2026-10-09 (see Pre-freeze resolutions at the end). Designed by judge-panel workflow wf_ef535660-fd4 (three independent designs, judge synthesis, adversarial data check with 22 corrections). No PREVENT, PCE or other risk value has been computed in NHANES III data.

## Question

Population: NHANES III adults aged 30-79 at visit 1 with no self-reported MI, stroke or heart failure. Each had a repeat mobile examination centre (MEC) examination 3-48 days after the first (median 16). The repeat used the same standardized protocol, the same MEC staff and the same examination session.

Primary question: how often does the base PREVENT 10-year ASCVD risk category (<3, 3-<5, 5-<10, >=10%) differ between the two visits when every routinely measured input is remeasured? The remeasured inputs are physician SBP, total and HDL cholesterol, serum creatinine/eGFR and laboratory diabetes status (HbA1c, exam-day insulin).

Key secondary questions:
(a) Does adding a remeasured urine ACR (UACR-enhanced model) change stability?
(b) How much of the instability comes from each input group?
(c) How does the chance of a change depend on distance from the nearest category threshold?
(d) How much of the laboratory part could analytical imprecision alone explain?

This is a test-retest (precision) study of a fixed equation. Neither visit is a reference standard. A category change is never called an error, a misclassification or a treatment change.

Within-person variation of the individual inputs in this second examination has been published before (Lacher 2005, Clin Chem; Selvin 2013, AJKD). The new contribution is propagating real, simultaneous remeasurement of all inputs through PREVENT.

## Population

Source: the frozen count gate (count_gate_v1/count_gate_results.json; variable_map.json SHA256 712e64a3...6c81). The reviewer reproduced every count below from the data files on 2026-10-09; all are non-outcome counts.

Eligibility flow (re-derived exactly in P3):
- 2,596 second-exam records
- 1,859 MEC second exams (examse PERTECH, cols 258-262, not blank)
- 1,859 linked to a MEC visit 1 (exam DMPSTAT=2 and MXPAXTMR present)
- 1,106 aged 30-79 at visit 1 (360 <= MXPAXTMR <= 959)
- 990 without self-reported MI, stroke or CHF (HAF10, HAC1D, HAC1C)
- 988 not pregnant at either visit (MAPF12R, MARF12R)
- 979 with determinable history (HAF10, HAC1D and HAC1C all = 2)
- 759 with complete base PREVENT inputs at both visits
- 725 who also have UACR at both visits

PRIMARY ANALYSIS SET: all 759 complete pairs. No one is excluded on input values; out-of-domain inputs are clamped (see primary_estimand).
- Female 420; current smokers 201; antihypertensive users 136 (none reached through a 'don't know' answer); statin users 10.
- Visit-1 age in completed years: 30-39, 190; 40-59, 311; 60-79, 258 (569 aged 40-79).
- Race-ethnicity (DMARETHN): 1 non-Hispanic white 362; 2 non-Hispanic black 205; 3 Mexican American 166; 4 other 26.
- Survey phase 1/2: 353/406.
- Interval (MXRRDAYS) known for all 759: 3-48 days, median 16.
- Same MEC session at both visits (MXPSESSR = MXRSESSR) for all 759.
- Same BP arm (PEP6B = PER6B) for 721; different arm for 38.
- Physician examiner codes: same 641, different 117, both blank-but-applicable 1 (13 distinct codes at visit 1, 11 at visit 2).
- Clustering: 98 masked pseudo-PSUs in 49 pseudo-strata. PSU size 1-14 (median 8); stratum size 3-24 (median 15). 8 pseudo-PSUs contain participants from both survey phases.

UACR SET: the 725 of the 759 with valid urine albumin and creatinine at both visits. Urine creatinine below detection (7.9) counts as incomplete. 97 pseudo-PSUs in 49 strata; one stratum has a single PSU. Across all second-exam pairs with both UACR values, UACR ranges from 0.18 to 16,064 mg/g.

Selection: NCHS second-exam documentation (EXAMSE-acc pp.10-11) says no statistical sampling design was applied. MEC staff selected a nonrandom ~5% at each location: mainly adults, about half aged 20-39 and half men. NCHS advises thorough investigation of selection bias. There are no weights. All primary estimates are unweighted and conditional on this sample, and Table 1 compares the sample with the eligible MEC population.

STANDARDIZATION TARGETS (K4 only):
- Base model: 9,478 MEC-examined NHANES III adults who meet the same visit-1 eligibility and have complete visit-1 base inputs under the same rules (8,173 with no second exam, 427 with a home second exam, 878 with a MEC second exam). All have WTPFEX6 > 0; 49 strata, 98 PSUs.
- UACR model: the 9,316 of these with a valid visit-1 UACR (8,025 / 423 / 868); 49 strata, 98 PSUs.

## Exposures and inputs

The exposure is repeat measurement under an identical protocol. Columns and fill rules follow variable_map.json (v1 plus the P1 v2 additions and superseding rules).

MEASURED AT EACH VISIT:
- SBP: mean of physician readings 2 and 3. Visit 1: PEP6H1 (cols 1403-1405) and PEP6I1 (1413-1415). Visit 2: PER6H1 (163-165) and PER6I1 (173-175). Never the composite PEPMNK/PERMNK variables.
- Lipids: TCP/TCR and HDP/HDR (mg/dL); non-HDL = TC - HDL inside PREVENT.
- Creatinine: CEP/CER, IDMS-calibrated at both visits as 0.960 x raw - 0.184 (applying it at visit 2 is an assumption). Then CKD-EPI 2021 (race-free, unrounded) with sex (HSSEX) and the visit-1 completed-year age at both visits.
- Diabetes at visit v = 1 if any of:
  - non-gestational self-reported diabetes (HAD1/HAD3/HAD4);
  - household insulin (HAD6) or diabetes pills (HAD10);
  - exam-day insulin at either visit (PHPINSU/PHRINSU = 1, treated as medication history);
  - HbA1c_v >= 6.5% (GHP/GHR). When HbA1c is missing at one visit, the other visit's value is used at both, so missingness cannot create a change.
  Otherwise 0, under the gate's self-report rules. This differs from the frozen gate rule, which is SA3c.
- UACR (mg/g) = 100 x UBP/URP at visit 1 and 100 x UBR/URR at visit 2. Albumin below detection is kept as 0.4. Urine creatinine of 7.9 counts as incomplete. URRSI is never used.

IDENTICAL AT BOTH VISITS: sex; age = floor(MXPAXTMR/12) (an integer, as in Paper 1, which also avoids pyprevent's NaN above 79); current smoking (HAR1/HAR3); antihypertensive use (HAE5A chain); statin use (pupremed HQRXCODE 1664/3064/3088).

DOMAIN CLAMP (primary): SBP 90-200, TC 130-320, HDL 20-100, eGFR 15-140; UACR 0.1-25,000. SA1b and SA1d use SBP 90-180. PCE clamp in S11: TC 130-320, HDL 20-100, SBP 90-200, age 40-79. Every model call goes through the guarded wrapper.

STRUCTURE VARIABLES (P1): SDPSTRA6, SDPPSU6, SDPPHASE, WTPFEX6 (targets only), MXRRDAYS, MXPSESSR/MXRSESSR, PEP6B/PER6B, PEPTECH/PERTECH, PEPLEVEL, GHPMETH/GHRMETH, DMARETHN (PCE and descriptive use only).

MODELS:
- PREVENT base 10-year ASCVD (primary).
- PREVENT UACR 10-year ASCVD (K1).
- PREVENT base 10-year total CVD (S10).
- Original 2013 PCE (P6 and S11 only).

## Primary estimand

theta = (1/759) x sum over i of 1[C(r_i2) != C(r_i1)]. This is the unweighted proportion of the 759 whose base PREVENT 10-year ASCVD category differs between visit 1 and visit 2.

Risk r_iv: pinned pyprevent prevent_base(...)['prevent_base_10yr_ASCVD'] (commit 69f34147), used only after P4 verification against the preventr base_10yr coefficient table (commit ebe2ac15). bmi is passed as NaN, which affects only the HF outputs. Risk = 100 x expit(LP), unrounded.

Categories: C(r) = searchsorted([3, 5, 10], r, side='right'), applied to unrounded risk, with inclusive lower bounds (the Paper 1 rule).

Inputs taken from their own visit:
- SBP: mean of physician readings 2 and 3.
- TC and HDL-C (mg/dL).
- eGFR: CKD-EPI 2021 (race-free), unrounded as in Paper 1. It uses creatinine calibrated as 0.960 x raw - 0.184 (calibrated value must be > 0; observed minimum 0.392), sex, and the visit-1 age in completed years.
- Diabetes status, under the rule in exposures_and_inputs.

Inputs identical at both visits: sex; age in completed years at visit 1, floor(MXPAXTMR/12); current smoking; antihypertensive use; statin use.

Domain rule (applied before scoring to each continuous input, at each visit and in every hybrid):
- Clamp to the developers' coded PREVENT domain: SBP 90-200 mmHg, TC 130-320 mg/dL, HDL-C 20-100 mg/dL. eGFR is clamped to 15-140, preventr's range; pyprevent itself rejects only eGFR <= 0.
- Sources: pinned pyprevent (AHAprevent parity, _core.py and _base.py) returns NA outside SBP 90-200, TC 130-320 and HDL 20-100. preventr estimate_risk.R (lines 410-419) states that the PREVENT equations technically permit SBP up to 200 mmHg and that preventr's 180 limit is the package author's clinical choice.
- UACR: 0.1-25,000 mg/g. No specimen falls outside this range.
- Non-outcome counts (reproduced): 49 of 759 have at least one clamped value (36 at visit 1, 33 at visit 2, 20 at both). Person-visits clamped (visit 1/visit 2): TC 20/18, HDL 7/8, SBP 7/6, eGFR 2/1. 710 are unclamped at both visits (676 in the UACR set, where 49 are clamped). 13 are in domain at visit 1 but out at visit 2.
- Rationale: clamping keeps the frozen denominator, treats both visits alike, and avoids excluding people on their remeasured values. It mirrors entering the nearest value a calculator accepts.
- Departures from Paper 1, which excluded people outside SBP 90-180: SA1a excludes under the primary bounds; SA1b clamps at SBP 90-180; SA1d reproduces Paper 1's exclusion domain exactly.

Input guards (asserted on every model call by a project wrapper):
- Integer age 30-79.
- Sex, diabetes, smoking, antihypertensive use and statin use in {0,1}.
- Continuous inputs finite and inside the clamp domain.
- On UACR-model calls, UACR finite and within [0.1, 25,000]. This guard is essential because pyprevent prevent_uacr silently substitutes a missing-UACR term when UACR is NaN.

Reporting: number of changers, theta, the upward and downward components, and the primary 95% CI (see uncertainty).

Pre-specified claim template: 'When every routinely measured PREVENT input was remeasured a median of 16 days later under the same standardized protocol, X of 759 adults (Y%, 95% CI a-b) changed base PREVENT-ASCVD category.'

Planning expectation only, not a hypothesis: theta of about 6-12%, with a 95% CI width of about 4-5 pp at a design effect of 1.0-1.3.

## Secondary estimands

- K1 (key) UACR-enhanced model, n=725. theta_U is the same estimand computed with the pinned prevent_uacr 10-year ASCVD model, with UACR also remeasured. That model is already verified against preventr uacr_10yr (maximum difference 3.55e-14 pp). The paired contrast is Delta_U = theta_U - theta_base, both computed in the same 725 people and reported with discordant-pair counts. It answers whether adding a remeasured, high-variability urine test makes the category less stable.
- K2 (key) Input attribution for the base model (4 input groups) and the UACR model (5 groups). Reported: the symmetrized Shapley split Psi_k of theta (headline); single-group swap proportions; hold-one-fixed proportions; LP-scale mean-squared-disagreement shares pi_k; and the {BP, LIPID} coalition (the Ulbricht analogue). Methods are in input_partition_method.
- K3 (key) Threshold-distance curves, for both models:
- pi_clin(d): probability of a category change given the visit-1 distance to the nearest threshold. This is the clinician's view and contains regression to the mean by design.
- pi_mech(d): the same given the distance of the visit-mean LP, which is free of regression to the mean.
- Instability zones: risk intervals around 3, 5 and 10% within which pi_clin >= 10% and >= 25%.
Methods are in threshold_conditional_and_standardization.
- K4 (key) theta standardized to the visit-1 risk distribution of the eligible NHANES III MEC population (WTPFEX6): n=9,478 for the base model and 9,316 for the UACR model. Reported with the ratio theta_std/theta. The transport assumption is labelled, and methods are in threshold_conditional_and_standardization.
- S5 Transition structure, base and UACR models:
- the full 4x4 table of visit-1 by visit-2 category (counts, row %, total %);
- P(up) and P(down), and the net shift (up minus down) with a paired CI;
- the number of non-adjacent moves;
- the category distribution at each visit;
- for t = 3, 5 and 10%, the proportion for whom t lies between r_i1 and r_i2.
For continuity with Paper 1, a fixed window of visit-1 risk within +/-0.5 pp of any threshold is also reported: the share of people in the window and the share of all changes that fall in it.
- S6 Continuous agreement on the LP (logit) scale, back-transformed where useful:
- mean Delta LP (systematic visit effect);
- SD(Delta LP) and the within-person SD, sigma_w = SD/sqrt(2);
- ICC(2,1) with absolute agreement;
- Bland-Altman bias and 95% limits of agreement on the logit scale, expressed as odds ratios;
- RCV of PREVENT risk: exp(+/-1.96 x SD(Delta ln r)) - 1, plus a robust version using 1.4826 x MAD;
- median and 90th percentile of |Delta r| (pp), and the proportions with |Delta r| >= 1 pp and >= 2 pp;
- heteroscedasticity: regression of |Delta LP| on the visit-mean LP, with cluster-robust SEs.
- S7 Input-level repeatability, inputs only (no risk). This replicates Lacher 2005 and Selvin 2013 in the analysis set. For SBP, TC, HDL-C, non-HDL-C, calibrated creatinine, eGFR, urine albumin, urine creatinine, UACR and HbA1c:
- the mean between-visit difference;
- total within-person CVw = sqrt(exp(s^2) - 1), with s = SD(Delta ln x)/sqrt(2); linear SDw for SBP and HbA1c;
- ICC and the 95% log-normal RCV.
Comparisons:
- Derived biological CVI = sqrt(max(CVw^2 - CVa^2, 0)), labelled model-based, is compared with the Lacher 2005 CVI for every analyte Lacher reports, and with EFLM CVI. Both sources are transcribed and hashed before scoring.
- Total CVw is compared with Selvin 2013's calibrated creatinine value of 7.6%, and is also computed after Selvin's exclusion of values more than 3 SD from the mean.
- Selvin's log-scale CVs (urine albumin 37.3%, urine creatinine 43.0%, ACR 11.3%) are shown beside ours, with the caveat that the definitions differ.
Also reported: the count of diabetes-status discordance between visits.
For SBP, the between-day variance is separated from the variance between readings within a visit: sigma^2_between-day = (Var(Delta mean23) - Var(reading2 - reading3)/2)/2, with the reading differences pooled across visits.
- S8 SBP measurement benchmarks. Other inputs are held at their visit-1 values. Each is computed in the base model, and in the UACR model on its 725:
(a) visit-1 single reading 1 versus single reading 2 (the Paper 1 analogue, 4.25%), excluding anyone lacking visit-1 reading 1 (5 or fewer);
(b) visit-1 reading 2 versus reading 3 (the Ogunwale analogue, 5.8%);
(c) visit-1 reading 2 versus visit-2 reading 2 (a single reading on two days);
(d) the between-visit SBP-only swap (mean of readings 2-3), taken from K2.
Paired contrasts are reported for (c - b) and (theta - d).
- S9 Systematic versus random variation. Each input's sample-mean visit-2 minus visit-1 shift is removed from the visit-2 values: additively for SBP and HbA1c, and on the log scale for laboratory analytes. HbA1c-defined diabetes is recomputed after the shift, and the data are rescored. theta minus theta_shift-removed is the part of the instability attributable to systematic visit effects, such as habituation or laboratory drift. In every bootstrap replicate the shifts are re-estimated and the data are rescored with a vectorized LP.
- S10 PREVENT base 10-year total CVD at the 7.5% threshold (2025 AHA/ACC hypertension guideline; scratch extraction gl2025/glance.txt). Reported: the proportion whose risk lies on opposite sides of 7.5% at the two visits, among all 759; the same in the UACR set, using the base model; and, descriptively, among the 84 untreated participants with visit-1 SBP 130-139. DBP and full blood-pressure staging are not analysed.
- S11 PCE bridge to Ye 2014, run only after P6 passes. Population: the 569 participants aged 40-79 at visit 1 (completed years). Equations: the original 2013 Pooled Cohort Equations; DMARETHN=2 uses the African-American equations and everyone else the white equations. Diabetes is visit-specific. PCE domain clamps: TC 130-320, HDL 20-100, SBP 90-200. Two outcomes are reported, crossing 7.5% and category change across the 2018 ACC/AHA categories (<5, 5-<7.5, 7.5-<20, >=20%), each under (a) SBP-only remeasurement, (b) full remeasurement of SBP, TC, HDL and diabetes, and the paired difference (b - a). Base PREVENT in the same 569 is shown alongside. No significance claims are made across equations.
- S12 Independent-error simulation benchmarks, always labelled model-based.
Set-up:
- Each person's pseudo-visits are centred on the visit mean: arithmetic for SBP and HbA1c, geometric for log-scale analytes.
- 1,000 replicates per person, seed 20261009, with common random numbers across bootstrap replicates. If the Monte Carlo SE is >= 0.05 pp, 5,000 replicates are used.
- Clamping is applied after simulation. Diabetes is recomputed from simulated HbA1c in people not fixed-diabetic.
(A) Observed-variation model: independent per-visit errors with SD = SD(Delta)/sqrt(2) for each input on its analysis scale, so that simulated between-visit differences reproduce the observed SD. Inputs are uncorrelated (the Cole-style assumption) and there is no rounding. Reported as theta_sim(A) and the ratio theta/theta_sim(A).
(B) Analytical-only model: laboratory analytes (TC, HDL, creatinine before calibration, HbA1c; plus urine albumin and creatinine in the UACR model) receive log-normal per-visit error with SD sqrt(ln(1 + CVa^2)), at the frozen CVa. SBP is fixed. Simulated values are then rounded to the files' reporting resolution: TC and HDL 1 mg/dL, creatinine 0.1 mg/dL, HbA1c 0.1%, urine albumin 0.1 mg/L, urine creatinine 0.1 mg/dL. Low and high CVa brackets are run. Reported as the ratio to the observed laboratory-only remeasurement proportion D_sym({LIPID, KIDNEY, GLYC}).
- S13 Conditional and descriptive rates. Category-change rates within each visit-1 category (Ye-style, affected by regression to the mean) are reported beside rates within each visit-mean category (not affected). Also reported: the regression-to-the-mean slopes of Delta LP on (LP1 minus its mean) and on (visit mean minus its mean). Ye-style conditional estimate: among those at 5-<10% at visit 1, the proportion changing category and the proportion crossing 5% or 10%. Descriptive only, with no claims: theta by sex, age band, race-ethnicity, antihypertensive use, interval tertile, survey phase, and same versus different examiner.
- E14 Exploratory (supplement only; labelled as an assumption-heavy transport across eras). pi_clin from K3 is applied to the base-model visit-1 PREVENT distribution of the Paper 1 NHANES 2009-2010 cohort (n=3,166; WTMEC2YR; inputs from final_and_external_v1/primary_participant_results.csv). Paper 1's SBP is the mean of available readings, which is flagged as a mismatch. This assumes 1988-94 within-person variability applies to 2009-10. It is not used for any claim about Paper 1.

## Input partition method

Input groups K (each is one routinely ordered test or procedure):
- BP: SBP, mean of physician readings 2-3.
- LIPID: TC and HDL-C together (one specimen and one run). TC-only and HDL-only splits are supplementary; they are additive on the LP scale.
- KIDNEY: calibrated serum creatinine, which gives eGFR (with visit-1 age).
- GLYC: visit-specific HbA1c >= 6.5%. Status differs between visits for 9 participants under the protocol rule.
- URINE: UACR, albumin and creatinine from the same specimen; UACR model only.
Age, sex, smoking, antihypertensive use, statin use and self-reported or treated diabetes are identical at both visits. They contribute exactly zero, which is reported as a design limitation, not a finding.

Hybrid risks. R_i(S) is the risk with the groups in S at their visit-2 values and the remaining groups at visit-1 values. Clamping is applied input by input, so every hybrid lies inside the domain. This needs 16 evaluations per person for the base model and 32 for the UACR model, all through the guarded wrapper.
- Reference visit 1: D1_i(S) = 1[C(R_i(S)) != C(R_i(empty set))].
- Reference visit 2: D2_i(S) = 1[C(R_i(K minus S)) != C(R_i(K))].

(i) Single swap: theta_k = mean of D1({k}). The visit-2-reference version and the average of the two (theta_k,sym) are also reported. BP-only is the Ye-type design; URINE-only is the Paper 1-type design.

(ii) Hold-one-fixed: theta_(-k) = mean of D1(K minus {k}), with both references.

(iii) HEADLINE: exact symmetrized Shapley split of the change indicator.
- phi_ik = sum over S not containing k of [|S|!(|K|-|S|-1)!/|K|!] x [D_i(S union {k}) - D_i(S)], computed for each reference visit; psi_ik is the average of the two.
- By efficiency, the psi_ik sum over k to 1[C_i1 != C_i2] for each person, so the Psi_k = mean over i of psi_ik sum exactly to theta.
- Psi_k is reported in pp with a CI, and as the share Psi_k/theta.
- Negative contributions (one input offsetting another's crossing) are kept and reported as offsetting.
- If theta has fewer than 20 events, only absolute Psi_k are reported and shares are not interpreted.

Paired quantities built from symmetrized indicators, such as 1[change] - theta_BP,sym,i, take values in [-1, 1] (including +/-0.5), not only {-1, 0, 1}.

(iv) Continuous partition (exact). In pinned pyprevent _base.py and _uacr.py (verified), every product term involves age, antihypertensive use or statin use: age x nonHDL, HDL, SBP_high, DM, smoking and eGFR_low; bptreat x SBP_high; statin x nonHDL. ln UACR is a main effect only. With those covariates fixed, the LP separates exactly across groups: Delta LP_i = sum over k of f_ik, with f_ik = LP_i({k}) - LP_i(empty set), identical in any substitution order. Reported:
- the mean signed contribution mu_k;
- the mean |f_k|;
- the mean-squared-disagreement share pi_k = E[f_k x Delta LP]/E[Delta LP^2], which sums to 1;
- the centred share nu_k = Cov(f_k, Delta LP)/Var(Delta LP);
- the correlation matrix of the f_k (for example, TC and HDL moving together).

Under SA6 (visit-specific age) and SA7 (visit-specific smoking), age or smoking becomes an extra group. Additivity then fails, and exact Shapley enumeration is used on both scales. Under SA6, hybrid eGFR is computed from the KIDNEY group's creatinine and the AGE group's age.

Wording rule: an input is called 'the largest contributor' only if it ranks first by Psi_k in >= 90% of cluster-bootstrap replicates; otherwise contributions are described as comparable. Single-swap, hold-one-fixed and Shapley results are always shown together, so that threshold interactions stay visible. The non-additivity gap, theta minus the sum of theta_k, is also reported.

## Threshold conditional and standardization

Distance to the nearest threshold:
- d_i1 = min over t in {3, 5, 10}% of |LP_i1 - logit(t/100)|.
- dbar_i uses the visit-mean LP, M_i = (LP_i1 + LP_i2)/2.
- Distances are also expressed in units of SD(Delta LP) and translated to the risk scale.

Curves (K3):
- Unweighted logistic regression of the change indicator on a restricted cubic spline of the distance (3 knots at its 10th, 50th and 90th percentiles, fixed across bootstrap replicates), plus a 3-level nearest-threshold indicator.
- With fewer than 40 events the indicator is dropped; with fewer than 20, only binned proportions are reported.
- Observed proportions in fixed bins (0-0.05, 0.05-0.10, 0.10-0.20, 0.20-0.40, >0.40 logit) are overlaid.
- Values are reported at d = 0, 0.05, 0.10, 0.20, 0.30 and 0.50, together with the instability zones.
- Diagnostic: the closed-form mechanism curve, the empirical survival function of |Delta LP| evaluated at 2d, is plotted against pi_mech.

Regression to the mean:
- theta, the partition, pi_mech and the standardization either condition on nothing or condition symmetrically, so regression to the mean does not bias them.
- pi_clin conditions on the first measurement, as a clinician must, and is labelled accordingly.
- A non-zero mean visit difference is a systematic visit effect (S9), not regression to the mean.
- Downward movement of high-risk people in tables conditioned on visit 1 is not presented as a finding.

Standardization (K4):
- theta_std = sum over j of w_j x pihat_clin(d_j1) / sum over j of w_j, over the target persons (9,478 base; 9,316 UACR), with w_j = WTPFEX6.
- d_j1 comes from the target person's single-occasion visit-1 PREVENT LP, computed with identical input rules: readings 2-3 SBP, clamping, completed-year age, and the gate's visit-1 diabetes rule. For the 759, this rule agrees with the protocol's visit-1 diabetes status in every case (0 differ), so the model inputs are consistent.
- Check estimator: post-stratification over 8 visit-1 risk bands (<2, 2-<3, 3-<4, 4-<5, 5-<7.5, 7.5-<10, 10-<15, >=15%). A band with fewer than 20 sample participants is merged with its neighbour toward 5%.
- Code identity check: applying the fitted model back to the analytic sample must reproduce theta to 1e-8. This holds by construction for a logistic MLE with an intercept.
- Untestable assumption, stated: given distance from the threshold, the probability of a change is the same in the target as in the convenience sample.
- The target includes the 878 MEC second-exam participants (868 for UACR), the analysis sample among them, so source and target overlap. The independent joint bootstrap is therefore approximate.
- Headline rule: if |theta_std - theta| > 2 pp, the abstract leads with the threshold-distance curve, not with crude theta as a generalizable number.
- The 2009-2010 transport is exploratory only (E14).

## Uncertainty

Inference is unweighted, with the person as the unit of analysis. The target of inference is the short-interval remeasurement process in persons exchangeable with this convenience sample, not the US population. Clusters are treated as the independent sampling units.

Clusters: the masked pseudo-PSU (SDPSTRA6, exam.dat cols 42-43, x SDPPSU6, col 41) serves only as a cluster identifier, standing for shared location, staff, season and laboratory runs. It is never used for weights or population inference, and strata have no design role in the primary analysis. Base set G=98; UACR set G=97. Caveat: 8 of the 98 pseudo-PSUs span both survey phases, so they are not pure location-time units; see SA15(d).

Primary CI for proportions (analytic):
- Cluster-robust (CR1) variance of the unweighted mean: V = [G/(G-1)] x sum over g of (sum over i in g of (y_i - theta))^2 / n^2.
- Design effect deff = V / [theta(1-theta)/n], floored at 1, so V* = max(V, theta(1-theta)/n).
- CI on the logit scale, using SE_logit = sqrt(V*)/(theta(1-theta)) and t with G-1 df (97 base, 96 UACR).
- With fewer than 10 events, or more than n-10: Korn-Graubard interval (Clopper-Pearson with effective n = n/max(deff, 1)).
- Always reported alongside: event count, Wilson CI and deff.

Paired contrasts (Delta_U, S8, S11, theta minus theta_BP-only or theta_BP,sym):
- Per-person differences d_i in [-1, 1]; CR1 variance of their mean, floored as V* = max(V_CR1, s_d^2/n); Wald CI with t(G-1).
- With fewer than 20 discordant pairs, an iid Tango score interval is reported as well, and the paired-contrast wording rule uses the wider interval.

Nonlinear statistics (Shapley and LP shares, curves, ICC and limits of agreement, RCV, CVw, S9, S12 ratios, K4):
- Nonparametric cluster bootstrap: resample the G clusters with replacement, unstratified.
- B=4,000, seed 20261009, 95% percentile intervals.
- Fixed scenarios (visits, hybrids, SA input definitions) are scored once.
- Everything that depends on sample-estimated quantities is recomputed in each replicate: categories and Shapley aggregates; spline refits with fixed knots; S9 shifts, with vectorized rescoring; SA10 propensity refits and weights; ratios. S12 uses common random numbers.
- B runs in 4 batches. If the batch-means MCSE of any primary-figure CI endpoint exceeds 0.1 pp, B rises to 10,000.

K4: joint bootstrap, B=2,000. The source side is the cluster bootstrap with refits. The target side is a Rao-Wu rescaling bootstrap within SDPSTRA6 (n_h - 1 PSUs per stratum, WTPFEX6 rescaled; WTPFEX6 has an embedded decimal). The two sides are treated as independent, which is approximate because of the overlap noted in K4.

SA15 (uncertainty methods, supplement):
(a) iid Wilson intervals.
(b) Unweighted stratified Taylor linearization, with df = PSUs minus strata after collapsing the UACR set's singleton-PSU stratum into its numerically adjacent stratum (49 for both sets).
(c) Percentile cluster bootstrap of theta.
(d) CR1 clustered on pseudo-stratum (G=49, sizes 3-24, t48). If its CI is more than 20% wider than the primary CI, the main text reports both.

Testing: no p-values and no 'significant' wording; estimates are given with CIs. There is one primary estimand and four key secondary families (K1-K4), whose CIs are unadjusted and labelled secondary. Everything else is supplementary or exploratory.

Precision label: 'adequately precise' if the primary CI width is <= 6 pp, otherwise 'imprecise'.

Paired-contrast wording: 'similar' only if the 95% CI lies within +/-2 pp; 'higher' or 'lower' only if the CI excludes 0; otherwise 'inconclusive'.

Not propagated, and stated as such: uncertainty in the creatinine calibration, the PREVENT and PCE coefficients, and the CVa values (S12 is bracketed instead).

Software (verified in the project .venv): Python 3.12.14, numpy 2.3.5, pandas 3.0.1, scipy 1.18.1, statsmodels 0.15.0, pyreadstat 1.3.6, rdata 1.1.0. Versions are frozen to environment_freeze.txt.

## Pipeline checks

- P0 Freeze.
- Save this protocol as studies/nhanes3_testretest/analysis_protocol_v1.md.
- Before any code reads data for analysis, record SHA256 of: the protocol; data/download_manifest.json, with all 17 data files re-hashed; method_sources/manifest.json; variable_map.json (v1 712e64a3...6c81 and the P1 v2); count_gate_results.json; the transcribed CVa, EFLM and Lacher 2005 tables.
- Outputs go to a new studies/nhanes3_testretest/analysis_v1/. The code refuses to overwrite it, and copies of the protocol and executed scripts are kept there.
- A deviation log records date, reason, and whether each deviation came before or after results were seen.
- Public registration of the protocol hash (OSF) is recommended; the user decides.
- P1 Variable map v2. Each column below was verified by the reviewer on 2026-10-09 against the .sas INPUT statements and codebook text. Each must be recorded with status 'verified', and the map re-hashed before freeze.
- exam.dat: DMARETHN 12 (1 non-Hispanic white, 2 non-Hispanic black, 3 Mexican American, 4 other); SDPPHASE 40; SDPPSU6 41; SDPSTRA6 42-43; WTPFEX6 59-67 (embedded decimal, e.g. '001737.58'); MXPSESSR 1234; PEP6B 1384 (arm selected); PEPLEVEL 1497 (1 = Level I, immediate attention; 2 = Level II, attention within two weeks; 3 = Level III; 8 = blank but applicable); PEPTECH 1498-1502.
- examse.dat: PER6B 144; PERLEVEL 257 (same codes).
- lab.dat: GHPMETH 1865.
- labse.dat: GHRMETH 594 (1-3 = Diamat instruments 1-3; 4 = affinity chromatography; 8 = blank but applicable); MXRSESSR 6.
Rules added in v2, superseding conflicting v1 rules:
- the protocol diabetes rule;
- the clamp rule with its sources, replacing the v1 'unknown' PREVENT-range rule;
- age: visit-1 completed-year age for PREVENT and eGFR at both visits; visit-2 age only in SA6;
- SBP: mean of readings 2-3, as executed in the gate, superseding v1's 'any valid reading' rule.
These positions were used only for protocol-stage non-outcome counts.
- P2 Independent reader.
- A second parser builds its column specifications automatically from the .sas INPUT statements, with no hand-typed positions. It must reproduce every analysis field exactly for all relevant SEQNs, including the K4 targets.
- Unit tests: CER '0001' = 1.0; UBR '00160' = 160; UBP '015700' = 15700; 8-fills = missing; URP/URR '007.9' = below detection; UBP/UBR 0.4 kept; WTPFEX6 '001737.58' = 1737.58.
- Full-file counts of fill codes must match the variable-map codebook counts.
- P3 Count reproduction (non-outcome). Any mismatch means STOP and a versioned reconciliation. Expected values:
- Flow: 2,596 / 1,859 / 1,859 / 1,106 / 990 / 988 / 979 / 759 / 725.
- Domain, primary bounds (SBP 90-200): 49 persons clamped at either visit (36 at visit 1, 33 at visit 2, 20 at both). Person-visits TC 20/18, HDL 7/8, SBP 7/6, eGFR 2/1. 710 unclamped at both (UACR 676; 49 clamped in the UACR set). 723 in domain at visit 1 (UACR 689). 13 in domain at visit 1 but out at visit 2.
- Domain, SBP 90-180 bounds: 56 / 41 / 36 / 21; SBP person-visits 13/9; 703 (670); 718 (685); 15; 55 clamped in the UACR set.
- Diabetes: 9 discordant under the protocol rule (10 under the frozen gate rule); HbA1c missing at exactly one visit 7; diabetes at both visits 90; visit-1 diabetes under the protocol rule differs from the visit-1-only gate rule in 0.
- Creatinine: minimum raw 0.6 mg/dL; minimum calibrated 0.392.
- Clustering: 49 strata / 98 PSUs for the 759 (PSU size 1-14, median 8; stratum size 3-24); 49 / 97 with one singleton stratum for the 725; 8 PSUs span both phases.
- Design features: same session 759; same arm 721; PEPLEVEL 3/2/1/missing = 705/50/2/2; visit-1 or visit-2 reading 1 missing 5; HbA1c method-class switch (Diamat versus affinity) 65; examiner same/different/both blank = 641/117/1.
- Subgroups: untreated with visit-1 SBP 130-139, 84; aged 40-79, 569.
- Sensitivity exclusions: Rose exertional pain 42; antianginal or coronary vasodilator users 9; women aged 30-59 with unknown pregnancy status 50.
- K4 targets: 9,478 base (8,173 / 427 / 878); 9,316 UACR (8,025 / 423 / 868).
- Ye cohort Y: 1,424 -> 1,190 -> 1,119 -> 1,078 -> 910; variants 920, 930 and 908 (CHF excluded); 24 aged 80+; none top-coded at 90+.
- P4 Equation verification, on synthetic inputs only, before any participant risk exists. Any failure means STOP.
(a) Base 10-year ASCVD and 10-year total CVD:
- Pinned pyprevent prevent_base (bmi=NaN) must agree with an independent implementation of the preventr base_10yr table to max |difference| < 1e-10 pp.
- The table is 23 rows; the columns used are female/male_ascvd and female/male_total_cvd; BMI rows 8, 9 and 21 are zero for these outcomes.
- Feature order is that of verify_prevent_uacr.py, without the ln UACR and missing-UACR rows.
- Use at least 1,200 profiles in both sexes, covering the spline knots (SBP 110/130; eGFR 60/90) and the domain bounds (SBP 90/200, TC 130/320, HDL 20/100, eGFR 15/140, age 30/79).
- BMI invariance: outputs identical for bmi = NaN, 5, 25 and 90.
(b) Re-hash and re-run the existing UACR verifications (method_verification.json; ascvd_input_verification_v2) without overwriting them.
(c) The CKD-EPI 2021 code (unrounded) is identical (< 1e-12) to the Paper 1 function (run_uacr_prevent_pilot.py line 125). The calibration 0.960 x raw - 0.184 is applied exactly once. Three hand-computed cases are stored.
(d) Unit tests for the clamp function at, inside and beyond every bound, including idempotence. Unit tests that the guarded wrapper refuses: NaN UACR on UACR-model calls (pyprevent would impute it), age > 79 or non-integer age, out-of-domain continuous inputs, and non-binary indicators.
(e) PCE: Goff 2013 coefficients are transcribed twice, independently, and the two transcriptions compared. The four published worked examples (age 55, TC 213, HDL 50, untreated SBP 120, nonsmoker, nondiabetic) must be reproduced within 0.05 pp; the values 2.1/3.0/5.3/6.1% are confirmed from the source before use. The pinned sysdata.rda holds no PCE tables (verified) and preventr's PCE source is not pinned, so no preventr comparison is claimed unless the user approves pinning it.
(f) On synthetic pairs: LP additivity (|Delta LP - sum of f_k| < 1e-10); Shapley efficiency (|sum of psi_ik - 1[change]| < 1e-12); the empty and full hybrids equal the visit-1 and visit-2 risks.
- P5 Input plausibility and linkage, before any scoring (inputs only).
- Linkage falsification: between-visit Pearson r >= 0.80 for TC, HDL and calibrated creatinine; >= 0.70 for SBP; >= 0.60 for ln UACR; >= 0.85 for HbA1c. Any failure: STOP and audit the SEQN joins.
- Flags, kept unless a parsing error is found and corrected with versioning: |Delta SBP| > 40 mmHg; |Delta TC| > 80 mg/dL; |Delta HDL| > 25 mg/dL; creatinine ratio outside 0.67-1.5; UACR ratio outside 0.1-10; |Delta HbA1c| > 1.5. The raw fixed-width text of flagged fields is listed.
- External consistency (triggers an audit, not a gate by itself): a ratio outside 0.5-2 triggers a parsing and unit audit for any of the following comparisons.
  - Calibrated-creatinine total CVw (ages >= 30, plus a version after Selvin's > 3 SD exclusion) versus Selvin's 7.6%.
  - Derived CVI for TC versus Lacher's 8.2%.
  - Derived CVI for HbA1c versus Lacher's 1.5%.
  - Urine albumin, urine creatinine and ACR are compared with Selvin descriptively only, because her log-scale CV definition differs.
- R0, Ye input audit in cohort Y (inputs only):
  - Median visit-1 and visit-2 SBP (mean of readings 2-3) versus Ye's 118 and 117 mmHg: a difference > 4 mmHg triggers the P6 audit list before any PCE is computed.
  - Median Delta SBP versus -1 mmHg: a difference > 2 mmHg triggers the same audit.
  - Descriptive only: mean age (46.1), % men (47.3), % non-Hispanic white (42.2), % smokers (29.2), % treated (14.5), mean TC (5.31 mmol/L) and HDL (1.31 mmol/L), and the 1st-99th percentiles of Delta SBP (-23 to 32).
- P6 BINDING gross-error screen: Ye et al. 2014 reproduction, PCE only. It runs after P0-P5 and before any PREVENT value exists. The scoring code refuses to run unless ye_reproduction.json records PASS, or FAIL with a user decision to continue.

Cohort Y, frozen now from non-outcome counts (n=910):
- MEC exam at both visits; visit-1 age >= 20 (MXPAXTMR >= 240), with no upper bound.
- Three valid physician SBP readings at each visit, in the same arm (PEP6B = PER6B).
- Complete visit-1 TCP, HDP, smoking, antihypertensive use, sex, and DMARETHN 1-4.
- Exclude HAF10=1 (MI, Ye's 'CHD'), HAC1D=1 (stroke) or HAD1=1 (diabetes).
- Require HAF10, HAC1D and HAD1 all = 2 (unknown history is treated as missing data).
Counts: 1,424 -> 1,190 -> 1,119 -> 1,078 -> 910. Variants are reported as counts only and never scored: 920 without the determinate-history rule; 930 excluding only non-gestational diabetes; 908 also excluding CHF.

Known non-reconstructability, stated in advance:
- Ye's equivalent step (three same-arm readings) had 956 people, against 1,119 here.
- Ye's '2,174 second visits' exceeds the 1,859 MEC second exams of all ages, so it must include home exams.
- Ye's coding is unpublished.
- Y overlaps heavily with Ye's 808, so the two-sample tolerances below are lenient. P6 detects gross pipeline errors; it does not establish exact reproduction.

Scoring:
- PCE at both visits, applied at all ages. This is inferred from Ye's inclusion of adults >= 20 (mean age 46.1, SD 17.2); Ye does not state it.
- SBP = mean of readings 2-3 at each visit; all other inputs from visit 1; diabetes = 0.
- DMARETHN=2 uses the African-American equations, everyone else the white equations.
- No clamping; unweighted; unrounded risk.

PASS requires all four. Each tolerance is |p - p_Ye| <= 1.96 x sqrt(pbar(1-pbar)(1/n + 1/n_Ye)):
- R1: crossing 7.5% (risk on opposite sides of 7.5% at the two visits) among those at 5-<10% at visit 1, versus 12/106;
- R2: crossing among the rest (<5% or >=10%), versus 2/702;
- R3: |Delta risk| >= 1 pp, versus 23% of 808;
- R4: share at 5-<10% at visit 1, versus 106/808.

Diagnostics (not binding, all reported):
- Ye-derived visit-1 counts: 58 at 5-<7.5% and 48 at 7.5-<10%.
- Transitions: 7 (12.1%) down and 6 (10.3%) up from 5-<7.5%; 6 (12.5%) down and 10 (20.8%) up from 7.5-<10%; 1.2% up from <5% (one into 7.5-<10%); 4.5% down from >=10% (one into 5-<7.5%).
- Median risk 2.5% (IQR 0.5-10.0) at visit 1 and 2.4% (0.5-10.0) at visit 2; median Delta risk 0.0 (1st-99th percentile -5.2 to +7.1).
- |Delta risk| >= 3 pp in 8% and >= 5 pp in 4%.
- Medians and IQRs of SBP at each visit; median Delta SBP.
- Y's full transition table in Ye's categories (<5, 5-<7.5, 7.5-<10, >=10%).

FAIL path:
1. Audit the shared components: BP fields, 888 handling, reading order, arm, lipid units, smoking and treatment coding, race mapping, PCE code.
2. Re-run with the P2 reader.
3. If an error is found: fix it with versioning and re-run P6.
4. If no error is found: record FAIL-UNEXPLAINED. PREVENT scoring proceeds only on the user's explicit decision in chat, entered in the deviation log, and the failure is reported in Methods, Results and Limitations.
The cohort definition may not be changed to approach Ye's numbers. The full reproduction table goes in the supplement whatever the result.
- P7 Single scoring run with exactness checks.
- All PREVENT and PCE estimands, including the K4 target LPs, are computed in one execution, through the guarded wrapper.
- A vectorized numpy LP equals scalar pyprevent for every visit, every hybrid and every target person (< 1e-10).
- All risks are finite and inside (0, 100); no UACR-model call receives a missing UACR.
- Per person: |Delta LP - sum of f_k| < 1e-10 and |sum of psi_ik - 1[change]| < 1e-12. The summed Psi_k equal theta.
- theta is unchanged when the visit labels are swapped.
- Export a participant-level CSV: SEQN, cluster, stratum, inputs and clamp flags for each visit, all hybrid risks, LPs, categories, f and psi; plus environment versions and hashes.
- P8 Independent recomputation.
- A separate standard-library Python script reads only the exported CSV. It must reproduce exactly: counts, theta, transition tables, single-swap, hold-one-fixed and Shapley means, the S10 and S11 indicators, and the CR1 variances (to 1e-10).
- Every exported risk is re-checked against the independent preventr-coefficient implementation (< 1e-10).
- Sanity checks: swapping no group gives 0 changes, and swapping every group reproduces theta.

## Sensitivity analyses

- SA1 [KEY: a] Domain handling.
(a) Exclude anyone with an out-of-domain input at either visit under the primary bounds: n=710 base, 676 UACR.
(b) Clamp with SBP 90-180 (the preventr and Paper 1 bound): 759/725; 56 persons clamped (55 in the UACR set).
(c) Clinical sequence: in domain at visit 1, visit-2 values clamped: 723 base, 689 UACR.
(d) Paper 1's convention exactly: exclude anyone outside SBP 90-180, TC 130-320, HDL 20-100 or eGFR 15-140 at either visit: 703 base, 670 UACR.
- SA2 [KEY: a] BP protocol, in the same 759.
(a) A single physician reading (reading 2) at each visit, representing office practice.
(b) Mean of all available readings 1-3 at each visit.
(c) Exclude the 38 pairs measured in a different arm.
- SA3 [KEY: a] Diabetes.
(a) Confirmed diabetes held fixed: 1 at both visits if self-reported, on medication, on exam-day insulin at either visit, or HbA1c >= 6.5% at both visits; otherwise 0 at both. GLYC is constant by construction.
(b) Exclude the 65 pairs whose HbA1c method class changed between visits (Diamat instruments 1-3 versus affinity 4).
(c) The frozen gate rule (missing HbA1c means no lab diabetes; insulin is visit-specific): 10 discordant.
- SA4 [KEY: a] Creatinine calibration.
(a) Raw Jaffe creatinine at both visits, which tests the assumption that the calibration also applies at visit 2.
(b) The Coresh 2002 correction (raw - 0.23 mg/dL) at both visits.
- SA5 [KEY] Referral-feedback guard. Exclude participants with visit-1 physician referral Level I or II, or unknown (PEPLEVEL 1, 2 or missing; 54 people, leaving n=705), because treatment may have started before visit 2. Always reported with the input-level mean Delta SBP (visit 2 minus visit 1) in the 54 and in the 705.
- SA6 Visit-specific age: completed years at each visit. Visit-2 age = floor((MXPAXTMR + MXRRDAYS/30.4375)/12), clamped at 79 (1 person reaches 80 at visit 2). Used in PREVENT and eGFR. Age becomes a Shapley player with exact enumeration; hybrid eGFR uses the KIDNEY group's creatinine and the AGE group's age.
- SA7 Visit-specific smoking: MAPA1 (visit 1, exam cols 5072-5074) and MARA1 (visit 2, examse cols 3308-3310) > 0 cigarettes in the past 5 days. Where these are 888, blank or DK, household status is used. Smoking becomes a Shapley player.
- SA8 Exclude pairs with any P5 implausible-change flag.
- SA9 Clinical eligibility. Exclude, separately and then together: Rose exertional chest pain without the site criterion (42); antianginal or coronary vasodilator users, prescription classes 0503/0508 (9); women aged 30-59 with unknown pregnancy status at either visit (50).
- SA10 [KEY] Completeness weighting from 979 to 759. A logistic model for having a complete pair among the 979 eligible, using age (spline), sex, race-ethnicity, phase, current smoking, antihypertensive use and self-reported diabetes; all are complete for the 979. Stabilized weights are truncated at the 1st and 99th percentiles, and the model is refitted in every bootstrap replicate. Visit-1 characteristics of included and excluded participants are always shown.
- SA11 Lipid-lowering proxy (HAE9D=1 or prescription class 0912) in place of the statin code list.
- SA12 Observer effect: same physician examiner at both visits only (641 pairs).
- SA13 Risk rounded to 0.1 pp before categorization, as the calculator displays it.
- SA14 UACR model only.
(a) Exclude pairs with albumin below detection at either visit.
(b) Keep urine-creatinine below-detection values (7.9 mg/dL), using the gate's sensitivity count of 732 pairs.
- SA15 Uncertainty methods: (a) Wilson; (b) stratified Taylor, df 49 for both sets; (c) cluster-bootstrap percentile; (d) CR1 clustered on pseudo-stratum (G=49, t48).
- Scope and reporting. Every SA targets theta. SA1, SA3, SA4 and SA5 are also run for theta_U and for the S10 total-CVD crossing. All SA results appear in one supplementary forest plot or table, whatever their direction. No SA replaces the primary estimate.

## Comparisons to prior work

All comparisons are descriptive. The studies differ in equation, era, thresholds and population, so no tests are made across them. Table 3 lists design, equation, inputs varied, interval, n and % changed for each.
- Ye 2014 (PCE, SBP only): screened at P6 and extended in S11 by SBP-only versus full remeasurement in the same equation and the same people, and by the PREVENT BP-only single swap (K2).
- Paper 1 (UACR random to first-morning, 5.34%; within-visit SBP reading 1 to 2, 4.25%): set beside the random-to-random URINE single swap (UACR model) and S8(a). The specimen contrasts, era and assay differ, so nothing is pooled.
- Ogunwale 2016 (ARIC, within-visit reading swap, 5.8%): set beside S8(b) and contrasted with the between-visit S8(c).
- Ulbricht 2018 (SCORE, TC and SBP remeasured): set beside the {BP, LIPID} coalition.
- Cole 2023 (simulated lipid biological variation, 6.3%): set beside the empirical LIPID single swap and S12(A), which tests whether an independent-error simulation reproduces the observed theta.
- Lacher 2005 (Clin Chem; biological variation of laboratory analytes in this same second exam) and Selvin 2013 (AJKD; within-person variation of kidney measures in this same second exam): input-level comparisons in S7 and P5, presented as replication.

## Displays

- Figure 1: STROBE flow from 2,596 to 759/725, with clamp counts under the primary bounds (49 persons), the K4 targets (9,478 base; 9,316 UACR) and Ye cohort Y (910, with the step comparison against Ye's 956/808).
- Figure 2: visit-1 versus visit-2 base PREVENT-ASCVD on log axes, with 3/5/10% gridlines and changers highlighted; inset 4x4 transition table; UACR panel.
- Figure 3 (key): attribution. Symmetrized Shapley Psi_k stacked to theta with CIs; dots for single-swap and hold-one-fixed proportions, with the Ye, Paper 1, Ogunwale and Cole benchmarks marked; LP-share panel. Base and UACR models.
- Figure 4: pi_clin and pi_mech against distance to the nearest threshold, with bootstrap bands, binned observed proportions and instability zones around 3, 5 and 10%.
- Figure 5: analyte variation. Observed total CVw and derived CVI versus frozen CVa, Lacher 2005 CVI and EFLM CVI for each analyte; SBP SD within versus between visits; S12 simulated versus observed theta.
- Table 1: visit-1 characteristics of the 759 beside the weighted 9,478 target and the 220 eligible non-completers (standardized differences); inputs at each visit with mean difference, CVw and RCV.
- Table 2: primary and key secondary estimates (theta, theta_U, Delta_U, theta_std, S10, S11), with CR1 and Wilson intervals, deff and n.
- Table 3: comparison with prior studies, including Lacher 2005 and Selvin 2013.
- Supplement: P4 verification report; P6 Ye screen table with variant counts and diagnostics; SA forest plot; full Shapley subset tables; TC/HDL split; S5-S9 and S13 tables; CVa/EFLM/Lacher source table; uncertainty-method comparison (including stratum clustering); E14; deviation log.

## Reporting and stopping rules

EXECUTION ORDER (each stage must pass before the next):
- Stage 1, no risk computed: P0 -> P1 -> P2 -> P3 -> P5 (including R0).
- Stage 2, PCE in cohort Y only: P4(e) -> P6.
- Stage 3: P4(a-d, f) -> P7 single scoring run -> primary estimate -> K1-K4 -> S5-S13 -> E14.
- Stage 4: sensitivity analyses -> P8 independent recomputation.
No PREVENT value exists before Stage 3.

STOP if:
- any of P0-P5 fails and cannot be fixed with a documented, versioned correction;
- P6 fails and the user has not decided, explicitly in chat, to continue after a FAIL-UNEXPLAINED audit;
- any P4 or P7 exactness or guard assertion fails.
A small, null or large theta never stops reporting.

FREEZE RULES:
- One scoring run.
- After scoring, no estimand, model, category, threshold, domain rule, input definition, cluster unit, B or seed may change.
- Post hoc analyses, including reviewer requests, are labelled post hoc with date and reason and shown beside the pre-specified results.
- Errors found after scoring are fixed in a new versioned directory, and both versions are kept. A corrected version becomes primary only for documented coding or parsing errors, never for changed analytic choices.

REPORT EVERYTHING:
- Every estimand, check and SA in this protocol is reported in the main text or supplement, in this order: primary, K1-K4, S5-S13, E14, SA1-SA15.
- Reporting standards: STROBE plus GRRAS. The overlap with Paper 1 (same program and author) is disclosed.

ABSTRACT must include:
- theta with its CI and precision label;
- up and down counts;
- theta_U and Delta_U under the +/-2 pp rule;
- the attribution statement under the 90% rank rule;
- the RCV of risk;
- theta_std;
- Ye screen status.
It must also state any KEY sensitivity analysis (SA1a, SA2a, SA3a, SA4a, SA5, SA10) that moves theta by >= 1.5 pp.

INTERPRETIVE RULES, fixed now:
(a) 'Laboratory remeasurement adds instability beyond BP' only if the CI of theta minus theta_BP,sym excludes 0. Otherwise: 'compatible with BP accounting for most of the instability', with the interval.
(b) S12(B) ratio (model-based analytical-only proportion divided by the observed D_sym({LIPID, KIDNEY, GLYC})):
- if its upper 95% limit is < 0.5, laboratory-driven change is described as mainly within-person biological or pre-analytical;
- if its CI includes 1, as compatible with analytical imprecision alone;
- otherwise, as partly analytical.
A definitive wording requires the low and high CVa brackets to agree; otherwise the wording is 'depends on the assumed analytical CV'.
(c) If S8(c) is below S8(b) by more than the half-width of its CI, the BP data are flagged (examiner effects, digit preference), and SBP instability is not attributed to day-to-day biological variation.
(d) Treatment or feedback between visits is invoked only if both hold: SA5 changes theta by >= 1.5 pp, and in the 54 excluded people the input-level mean visit-2 minus visit-1 SBP is negative, with a 95% CI excluding 0, and lower than in the other 705. Otherwise any SA5 change is described as a composition effect.
(e) Psi_GLYC is always described as the effect of re-applying a single HbA1c >= 6.5% criterion at each visit, not as a change in clinical diagnosis. SA3a holds GLYC constant by construction and is read through theta: if |theta(SA3a) - theta| >= 1.5 pp, the abstract reports both.
(f) The URINE contribution is described as sensitive to the handling of below-detection albumin in the 1988-94 immunoassay only if, under SA14a, Psi_URINE changes by >= 50% relative AND >= 0.5 pp absolute.
(g) The headline rule for theta_std is in threshold_conditional_and_standardization.

LANGUAGE:
- Say 'category change' or 'test-retest disagreement'; never 'misclassification', 'error' or 'treatment change'.
- S12 results are always called model-based.
- Subgroup results are descriptive. Proportions failing NCHS presentation standards (Korn-Graubard width and df rules, n < 30) are flagged.
- S7 is described as a replication of Lacher 2005 and Selvin 2013.

PAPER 1 BOUNDARY: the random-to-random URINE single swap and the S8(a) benchmark are shown beside Paper 1's 5.34% and 4.25% only as described comparisons. They are not pooled and are not used to revise Paper 1 before its submission.

NOVELTY CLAIM: immediately before submission, a dated PubMed/Scholar search is run: (PREVENT OR 'pooled cohort') AND (repeat OR remeasure* OR test-retest OR variability) AND (reclassif* OR category). The 'first rescoring after real remeasurement of all inputs' claim is softened if anything matches. Input-level biological variation in this dataset is never claimed as new.

RELEASE: code and derived data go to GitHub/Zenodo only with the user's explicit authorization.

## Threats and limitations

- Selection. The second exam is a nonrandom ~5% subsample chosen by staff, with no sampling design and no weights, and 220 of the 979 eligible lack complete pairs. Crude theta depends on how many people sit near thresholds, so it is not a prevalence. K4, SA10 and the band-specific curves rest on conditional-exchangeability assumptions that cannot be tested.
- Direction relative to routine care is mixed. Pushing toward less instability than in practice: the same staff and session, standardized physician mercury BP (mean of readings 2-3, recorded to even numbers), a short interval, and presumably the same laboratories. Pushing toward more: 1988-94 assays (kinetic Jaffe creatinine reported to 0.1 mg/dL; urine-albumin immunoassay with inter-assay CV up to about 16%). The BP component is therefore probably an underestimate for office practice, while the laboratory component may overstate contemporary analytical noise.
- Carried inputs. Smoking, antihypertensive use, statin use, self-reported diabetes and history come from a single household interview before visit 1. They contribute zero by construction, so questionnaire instability and medication changes are omitted (SA7 covers part of smoking). 52 participants referred at Level I/II may have started treatment between visits (SA5).
- Historical assays:
- the creatinine calibration was derived in adults aged 60+, is extrapolated to low values, and is assumed to apply at visit 2 (SA4);
- the urine creatinine instrument is disputed (Selvin reports a Beckman Synchron AS/Astra);
- HbA1c was measured by Diamat HPLC, with 65 pairs switching to or from the affinity method.
The analytical share may therefore exceed what contemporary assays would give, whereas biological variation is not specific to the era.
- Total within-person variation only. No duplicate assays exist, so the analytical/biological split rests on documented CVa and an independent, constant-CV log-normal model (S12, S7), with bracketed CVa sources.
- Prior input-level work. Lacher 2005 and Selvin 2013 already estimated within-person variation of these analytes in this dataset; S7 is a replication, and only the propagation to risk categories is new.
- Domain handling is not neutral. Clamping attenuates variation for 49 people under the primary bounds (56 under SBP 90-180), and excluding them (SA1a, SA1d) selects on remeasured values; both are shown. The primary SBP bound follows the developers' coded range (200), not Paper 1's 180.
- A single HbA1c >= 6.5% is not a clinical diagnosis, yet it carries a large coefficient. Status flips for 9 people (SA3).
- Masked pseudo-PSUs only approximate survey locations: 8 span both phases, and examiners cut across locations. Using design variables as cluster identifiers departs from the letter of NCHS advice, though not from its purpose (weighting). Wilson, stratified and stratum-clustered intervals are shown.
- Few events (likely tens) limit the precision of the category-scale Shapley shares, the curves and the paired contrasts. 'Inconclusive' is a likely wording for Delta_U. The LP-scale partition is precise.
- Representation: statins are rare (10 of 759; atorvastatin and rosuvastatin are absent), so the statin interaction is hardly exercised. The race-ethnic mix reflects NHANES III oversampling.
- There are no outcomes. Instability does not show that either visit predicts better. PREVENT was derived from single measurements, so regression dilution is partly built into its calibration.
- Guideline mapping is partial. Decisions also depend on LDL-C, risk enhancers, CAC, DBP and shared decision-making. Thirty-year, HF and HbA1c/SDI/full PREVENT models are not analysed.
- The Ye check cannot be an exact reproduction. Ye's same-arm step had 956 people versus 1,119 here; Y is 910 versus 808; Ye's coding is unpublished; the samples overlap, so the tolerances are lenient; and the PCE is extrapolated to ages 20-39 and 80+ (reproduction only).
- Age is held at visit 1 by design (SA6 relaxes this). The estimand excludes ageing, which is negligible over 3-48 days.

## Data facts relied on

- Frozen gate (count_gate_results.json), reproduced by the reviewer on 2026-10-09:
- Flow: 2,596 / 1,859 / 1,859 / 1,106 / 990 / 988 / 979 / 759 / 725.
- Interval 3-48 days, median 16, known for all 759.
- Missing among the 979 (visit 1/visit 2): TC 39/93, HDL 45/100, creatinine 52/99, SBP readings 2-3 41/104.
- Sensitivity counts: any-reading SBP 839; UACR keeping 7.9 values 732.
- Domain counts (judge_counts.py re-run and review_counts.py; eGFR used only for domain status; no risk computed):
- 1 of 759 has MXPAXTMR >= 949, and 1 reaches 80 at visit 2 on decimal age.
- SBP 90-200: 49 / 36 / 33 / 20; person-visits TC 20/18, HDL 7/8, SBP 7/6, eGFR 2/1; 710 (676); 723 (689); 13 in at visit 1 but out at visit 2.
- SBP 90-180: 56 / 41 / 36 / 21; SBP person-visits 13/9; 703 (670); 718 (685); 15.
- UACR among all second-exam pairs with both values: 0.18-16,064 mg/g.
- Minimum raw creatinine in the 759: 0.6 mg/dL (0.392 calibrated).
- Diabetes: 10 discordant under the gate rule (1 caused by HbA1c missing at one visit); 9 under the protocol rule; HbA1c missing at exactly one visit 7; diabetes at both visits 90. Visit-1 diabetes under the protocol rule equals the visit-1-only gate rule for all 759.
- Design structure:
- 759: 49 strata / 98 PSUs, no singleton strata; PSU size 1-14 (median 8); stratum size 3-24 (median 15); 8 PSUs span both phases.
- 725: 49 strata / 97 PSUs, one singleton stratum.
- Phase 353/406. Same session 759/759. Same arm 721, different 38. Examiner codes same 641, different 117, both blank-applicable 1.
- PEPLEVEL 3:705, 2:50, 1:2, missing 2. Reading 1 missing at either visit: 5. HbA1c method-class switches: 65.
- Composition and subgroups:
- Untreated with visit-1 SBP 130-139: 84. Aged 40-79: 569. Age 30-39 / 40-59 / 60-79: 190 / 311 / 258.
- Female 420; smokers 201; statin users 10; BP-treated 136; antihypertensive 'don't know' path 0.
- DMARETHN 362/205/166/26.
- Rose exertional pain 42; antianginal classes 0503/0508, 9; women aged 30-59 with unknown pregnancy status, 50.
- IPW covariates complete for all 979.
- Standardization targets and Ye cohorts:
- Base target 9,478 (8,173 / 427 / 878), all WTPFEX6 > 0, 49 strata / 98 PSUs.
- UACR target 9,316 (8,025 / 423 / 868), 49 / 98.
- Ye-style flow 1,424 -> 1,190 -> 1,119 -> 1,078. Variant A (excluding HAF10/HAC1D/HAD1=1) 920: 43.3% non-Hispanic white, 47.4% male. Non-gestational diabetes only: 930. A with CHF also excluded: 908. A with determinate history: 910 (frozen cohort Y). 24 aged 80+; none top-coded.
- Column positions and codes verified against the .sas INPUT statements and codebooks:
- exam: DMARETHN 12 (1-4 = NH white, NH black, Mexican American, other); SDPPHASE 40; SDPPSU6 41; SDPSTRA6 42-43; WTPFEX6 59-67 (embedded decimal); MXPSESSR 1234; PEP6B 1384; PEPLEVEL 1497 (1/2/3 = Level I/II/III, 8 = blank but applicable); PEPTECH 1498-1502; PEP6G1/H1/I1 1393-1395/1403-1405/1413-1415.
- examse: PER6B 144; PERLEVEL 257; PERTECH 258-262; PER6G1/H1/I1 153-155/163-165/173-175; MXRRDAYS 3727-3729.
- lab: GHPMETH 1865.
- labse: GHRMETH 594 (1-3 Diamat, 4 affinity); MXRSESSR 6.
- pyprevent (commit 69f34147):
- validate_common_inputs returns NaN if age < 30 or > 79, SBP < 90 or > 200, or eGFR <= 0. There is no eGFR upper bound and no lower bound at 15.
- prevent_base and prevent_uacr also return NaN outside TC 130-320 or HDL 20-100, or when statin or bptreat are missing.
- UACR in [0, 0.1) is floored at 0.1. A NaN UACR does NOT give NaN: the missing-UACR term is used.
- Both functions require bmi; BMI affects only the HF outputs.
- Every product term in the base and UACR 10-year ASCVD and total-CVD models involves age, bptreat or statin, and ln UACR is a main effect only, so the additive LP partition is exact.
- preventr (commit ebe2ac15):
- helpers.R valid ranges: SBP 90-180, TC 130-320, HDL 20-100, eGFR 15-140, UACR 0.1-25,000, HbA1c 4.5-15.
- estimate_risk.R lines 410-419: the PREVENT equations technically permit SBP up to 200; the 180 limit is the package author's choice.
- calc_egfr rounds eGFR to an integer (not used here).
- sysdata.rda contains base_10yr (23 x 11: female/male total_cvd, ascvd, heart_failure, chd, stroke; BMI rows zero for total_cvd and ascvd), uacr_10yr (25 rows), the hba1c/sdi/full tables, the 30-year tables, sdi_decile_dat and valid_zips. There are no PCE coefficient tables, and the PCE R source is not pinned.
- Paper 1 code (run_uacr_prevent_pilot.py): integer age (RIDAGEYR) in PREVENT and CKD-EPI 2021; unrounded eGFR (line 125); SBP = mean of available readings; domain SBP 90-180, TC 130-320, HDL 20-100, BMI 18.5-39.9, eGFR 15-140, age 30-79 (exclusion); categories via searchsorted([3, 5, 10], side='right').
- Existing verification: UACR 10-year ASCVD agrees with the preventr coefficients to a maximum of 3.55e-14 pp, with 2 vignette rows within 0.05 pp; v2 proved the UACR model invariant to BMI. Base ASCVD and total CVD are not yet verified (P4a).
- The 2025 AHA/ACC hypertension guideline uses 10-year PREVENT CVD >= 7.5% (scratch gl2025/glance.txt, lines 129-281).
- Variable map: applying the CER calibration at visit 2 is an assumption (no LABSE note). PERTECH not blank defines a MEC second exam. MAPA1 (exam 5072-5074) and MARA1 (examse 3308-3310) are verified. Statin codes 1664/3064/3088. Composite BP variables must not be used.
- Literature values verified from PubMed full text on 2026-10-09:
- Ye 2014 (PMC3959282): n 808; 12/106 versus 2/702; |Delta risk| >= 1/3/5 pp in 23%/8%/4%; median risk 2.5% and 2.4% (IQR 0.5-10.0); mean age 46.1 (SD 17.2); 47.3% men; 42.2% non-Hispanic white; 29.2% smokers; 14.5% treated; TC 5.31 and HDL 1.31 mmol/L; SBP 118 (108-131) and 117 (107-129); Delta SBP median -1 (-23 to 32); interval median 17 days (1-48); flow 2,174 / 956 / -32 / -116 / 808.
- Selvin 2013 (PMC3628297): total CVw for standardized creatinine 7.6% (ages >= 18, > 3 SD extremes excluded); log-scale CV for ACR 11.3%, urine albumin 37.3%, urine creatinine 43.0%; inter-assay CVs serum creatinine 2.0-2.7%, urine albumin 4.8-16.1%, urine creatinine 1.5-7.7%; Lacher biological CVs TC 8.2% and HbA1c 1.5% (analytical removed).
- Lacher 2005: Clin Chem 51:450-2, doi 10.1373/clinchem.2004.039354; table to be transcribed.
- Unverified descriptive comparators: Paper 1 5.34% (n=3,166) and 4.25%; Ogunwale 5.8%; Cole 6.3%.

## Open questions

- Before freeze: confirm the SBP bound against the primary source (Khan 2024 supplement or the AHA PREVENT calculator). The protocol now uses 90-200, on the evidence of the AHAprevent-parity code and preventr's own documentation. Revert to 90-180 (swap with SA1b and restore the expected counts 56/703) only if the primary source explicitly documents 90-180.
- Before freeze: confirm the category cut-points in the 2026 dyslipidemia guideline (DOI 10.1161/CIR.0000000000001423). Paper 1's <3/3-<5/5-<10/>=10% stay primary; if the guideline differs, add its categories as a secondary analysis before freeze.
- Before P4(e): confirm the Goff 2013 coefficients and the worked examples (2.1/3.0/5.3/6.1%) from the primary source. Optionally, with the user's approval, download and pin preventr's PCE R source at commit ebe2ac15 for an independent check.
- Before freeze: transcribe and hash the Lacher 2005 CVI table (Clin Chem 51:450-2) and read Selvin's cited CVw method, to fix the log-scale definitions used in the P5 and S7 comparisons.
- Before freeze: build the frozen CVa table and hash it: HbA1c from the LABSE documentation (Diamat inter-assay CV about 2.0%); serum and urine values from Selvin 2013; TC and HDL from the NHANES III Laboratory Procedures Manual, or the NCEP 3%/4% goals, labelled as such. Transcribe EFLM CVI entries with an access date.
- Before P6, optional: Ye's tables are not in the retrieved PMC text. All binding R1-R4 values and the diagnostics come from the verified text.
- Before submission (not before freeze): verify Ogunwale 2016 (5.8%), Cole 2023 (6.3%) and Ulbricht 2018 from the primary papers.
- The user decides whether to register the protocol hash publicly (OSF) before scoring, and which target journal to aim for; the journal affects display limits, not the analysis.
- Reviewer acceptability of pseudo-PSU clustering given NCHS advice on design variables: check the target journal's statistical guidance. The primary interval floors deff at 1, and Wilson, stratified and stratum-clustered intervals are always shown.

## Adversarial check corrections applied

- **Domain rule: SBP bound (primary_estimand, open_questions, P3, SA1, threats, displays)**: The protocol kept SBP 90-180 as the primary clamp until the developers' range was checked, but the pinned sources already settle it. AHAprevent-parity code accepts SBP 90-200. preventr's 180 limit is its author's own clinical choice, and preventr's documentation says so. Every downstream count (56 clamped, 703 unclamped) was tied to the 180 bound. Fix: Following the protocol's own swap rule: the primary clamp is now SBP 90-200 (TC 130-320, HDL 20-100, eGFR 15-140 unchanged). SA1b clamps at 90-180. A new SA1d reproduces Paper 1's exclusion domain exactly. All counts are updated. A residual open question allows a revert before freeze only if Khan 2024 or the AHA calculator explicitly documents 90-180.
- **pyprevent input handling (primary_estimand guards, P4d, P7)**: The protocol assumed out-of-domain or missing inputs give NaN, which the 'all risks finite' check would catch. In fact prevent_uacr does not return NaN for a missing UACR: it silently swaps in the missing-UACR model term, so a NaN UACR in a hybrid, SA or target row yields a plausible wrong risk. pyprevent has no eGFR 15-140 check (it rejects only eGFR<=0). Both functions also require a bmi argument. Fix: A wrapper now asserts integer age 30-79, binary indicators in {0,1}, finite in-domain continuous inputs, and a finite UACR in [0.1, 25,000] on every UACR-model call. P4(d) tests that the wrapper refuses NaN UACR and out-of-domain values. BMI is passed as NaN, and P4(a) proves base ASCVD and total CVD invariant to BMI (NaN, 5, 25, 90), as v2 already did for the UACR model.
- **P6 Ye 2014 reproduction: cohort reconstruction and tolerance**: (i) Ye's published flow (956 with three same-arm readings, then -32 missing, then -116 CHD/diabetes/stroke, giving 808) does not match cohort Y at the equivalent step (1,119 vs 956), so Y cannot reconstruct Ye's cohort. (ii) The two-sample tolerances treat heavily overlapping samples as independent, so P6 is a lenient gross-error screen, not a reproduction test. (iii) 'PCE applied at all ages' is inferred (adults >=20, mean age 46.1, SD 17.2), not stated. (iv) Many checkable figures in Ye's text were not used as diagnostics. (v) Ye's 2,174 'second visits' exceeds the 1,859 MEC second exams of all ages, so it must include home exams. Fix: R1-R4 are kept as binding but relabelled a gross-error screen, with the overlap caveat. Added: the step-by-step count comparison; derived category counts (58 at 5-<7.5%, 48 at 7.5-<10%); transition proportions; the text diagnostics, with numeric values; a non-binding Ye-category transition table; an input-only audit (R0) in P5. The inference about age is stated as an inference.
- **P5 external-consistency benchmarks (Selvin/Lacher)**: Lacher's TC 8.2% and HbA1c 1.5% are biological CVs with the analytical component removed, yet P5 compared them with total CVw. Selvin's ACR 11.3% (urine albumin 37.3%, urine creatinine 43.0%) are 'log-scale' CVs, not comparable with sqrt(exp(s^2)-1), so the 0.5-2 ratio trigger would fire spuriously. Selvin's creatinine 7.6% excluded values >3 SD from the mean and used ages >=18. Fix: TC and HbA1c are now compared as derived CVI = sqrt(max(CVw^2 - CVa^2, 0)) against Lacher. The audit trigger keeps calibrated creatinine (also computed after Selvin's >3 SD exclusion), TC CVI and HbA1c CVI. Urine albumin, urine creatinine and ACR are reported descriptively only. Selvin's CVw method reference is read before freeze.
- **Missing primary prior work: Lacher 2005**: Lacher, Hughes and Carroll estimated the biological variation of laboratory analytes from this same NHANES III second examination, in Clinical Chemistry (a target journal). The protocol cites it only second-hand. S7 would then be presented without its closest precedent. Fix: Cite Lacher 2005 directly. Transcribe and hash its CVI table before freeze. Frame S7 as replication. Add it to comparisons_to_prior_work and Table 3. The novelty claim stays limited to rescoring the risk equation.
- **Interpretive rule (e): GLYC and SA3a**: SA3a holds diabetes fixed across visits, so GLYC's contribution is zero by construction. The rule 'if GLYC falls by >50% under SA3a' therefore always fires and carries no information. Fix: Psi_GLYC is always described as the effect of re-applying a single HbA1c >=6.5% criterion at each visit. SA3a is read through theta: if |theta(SA3a) - theta| >= 1.5 pp, the abstract reports both values.
- **Interpretive rule (b): S12(B)**: The two conditions can hold together: a ratio point estimate <0.5 with a CI that includes 1. No precedence is set, and the result is not required to hold across the low and high CVa brackets. Fix: Precedence is now set by CI limits: upper limit <0.5 gives 'mainly biological/pre-analytical'; a CI including 1 gives 'compatible with analytical imprecision'; otherwise 'partly analytical'. A definitive wording requires both CVa brackets to agree.
- **Interpretive rule (d): SA5 attribution**: Removing 54 people can move theta by >=1.5 pp through composition alone, because referred people have high risk and high SBP. Such a move does not show treatment between visits. Fix: Treatment or feedback is invoked only if theta moves by >=1.5 pp AND the input-level mean visit-2 minus visit-1 SBP in the 54 is negative (CI excluding 0) and below that of the other 705. Otherwise the change is described as a composition effect.
- **Interpretive rule (f): SA14a**: The link from 'contribution falls under SA14a' to 'specific to the assay era' is not justified, and a 50% relative threshold on a small Psi is noisy. Fix: Reworded: the URINE contribution is described as sensitive to the handling of below-detection albumin in the 1988-94 immunoassay only if Psi_URINE changes by >=50% relative AND >=0.5 pp absolute.
- **Uncertainty: rescoring inside bootstrap replicates**: 'Person-level risks are fixed' contradicts S9 (shifts re-estimated per replicate, which requires rescoring), SA10 (propensity model refitted) and S12 (simulation). Fix: Risks are fixed for fixed scenarios (visits, hybrids, SA input definitions). Statistics that depend on sample-estimated quantities are recomputed per replicate, using vectorized LP rescoring for S9.
- **Paired contrasts**: theta - theta_BP,sym has per-person values in {-1, -0.5, 0, 0.5, 1}, not {-1, 0, 1}. No design-effect floor was set, and Wald intervals are poor with few discordant pairs. Fix: Per-person differences can lie anywhere in [-1, 1]. V* = max(V_CR1, s_d^2/n). With <20 discordant pairs, an iid Tango score interval is reported beside the cluster-robust interval, and the wording rule uses the wider of the two.
- **SA15 Taylor degrees of freedom**: After the singleton-PSU stratum is collapsed in the UACR set, df = 97 PSUs - 48 strata = 49, not 48. Fix: df is now defined as PSUs minus strata after collapsing: 49 for both sets.
- **Cluster unit validity**: Pseudo-PSUs are masked units, not single location-time stands: 8 of the 98 pseudo-PSUs in the 759 contain participants from both survey phases, and large certainty locations may be split across PSUs. Correlation between PSUs within a stratum would then be ignored. Fix: Added SA15(d): CR1 clustered on pseudo-stratum (G=49, sizes 3-24, t48). If its CI is >20% wider than the primary CI, the main text reports both.
- **S12 simulation specification**: The error SD was ambiguous: drawing each pseudo-visit with SD equal to SD(Delta) doubles the variance. The analytical-only model ignored reporting resolution (creatinine at 0.1 mg/dL is about 3% SD at 1 mg/dL). Diabetes recomputation and clamping order were unspecified. Fix: Per-visit SD is now SD(Delta)/sqrt(2) on the analysis scale for model A, and sqrt(ln(1+CVa^2)) for model B. Model B rounds simulated values to reporting resolution (model A does not, since resolution is already in the observed SD). Diabetes is recomputed from simulated HbA1c in people not fixed-diabetic. Clamping is applied after simulation.
- **SA6 hybrid eGFR**: Once age becomes a Shapley player, eGFR depends on both the creatinine visit and the age visit, so hybrid eGFR is undefined. Fix: Hybrid eGFR uses the KIDNEY group's creatinine and the AGE group's age.
- **K4 for the UACR model: target undefined**: K4 is specified 'for the base and UACR models', but the 9,478 target requires only base inputs. A UACR-model standardization needs visit-1 UACR in the target, and that population had not been defined or counted. Fix: UACR target defined as the 9,478 with a valid visit-1 UACR (urine creatinine not below detection): n=9,316 (8,025 no second exam, 868 MEC, 423 home), 49 strata, 98 PSUs. Added to P3.
- **K4 visit-1 diabetes rule and source-target overlap**: The sample's visit-1 diabetes rule can use visit-2 insulin and carry HbA1c over from visit 2, while the target uses visit-1 data only. The 878 MEC second-exam members of the target include the 759, so source and target are not independent. Fix: Rule agreement is verified now (0 of 759 differ) and fixed as a P3 expected value. The overlap is stated, and the joint bootstrap is described as approximate.
- **eGFR computation details**: Rounding was unspecified: preventr's calc_egfr rounds eGFR to an integer, while Paper 1 does not. No guard existed against non-positive calibrated creatinine. Fix: eGFR is unrounded, as in Paper 1. Calibrated creatinine must be >0 (the observed minimum is 0.392, from raw 0.6).
- **Threats: 'lower bound' claim**: Calling the result 'probably a lower bound for real-world instability' contradicts the stated point that 1988-94 assays (Jaffe creatinine at 0.1 mg/dL resolution, older urine immunoassay) may add more analytical noise than contemporary assays. Fix: Reworded as a mixed direction: the BP component is probably an underestimate for office practice, while the laboratory component may overstate contemporary analytical noise.
- **P1 codebook facts and variable-map rule conflicts**: The PEPLEVEL code meanings were left open. GHPMETH/GHRMETH, WTPFEX6 and DMARETHN codes were not stated. Three variable-map rules conflict with the protocol: SBP averaging ('any valid reading'), age for eGFR (visit-2 age) and the diabetes rule. The map's PREVENT-range rule still has status 'unknown'. Fix: All P1 columns and codes are verified now and recorded. Map v2 marks the protocol rules as superseding the conflicting map rules. The PEPLEVEL open question is closed.
- **Literature values in data_facts and open_questions**: The Ye and Selvin values were listed as 'read through a summarising fetch, to be re-read'. They are now verified from full text. Ogunwale, Cole and Ulbricht remain unverified. Fix: Ye and Selvin are marked verified (PMC full text, 2026-10-09). An open question is added to verify Ogunwale 5.8%, Cole 6.3% and Ulbricht from the primary papers before submission; they are descriptive only, so this is not needed before freeze.
- **Verification of remaining protocol facts (no change needed)**: None. Listed so the reviewer's checks are traceable. Fix: No change.

## Pre-freeze resolutions (2026-10-09, workflow wf_29a307eb-edb; facts in prefreeze_facts_v1.json)

- SBP domain: 90-200 mmHg confirmed by the developer paper (Khan et al., Circulation 2024, Methods, Study Population: exclusion SBP <90 or >200, TC <130 or >320, HDL-C <20 or >100) and the official AHA PREVENT calculator code (SBP validation <90 or >200; eGFR 15-140; age 30-79). Primary clamp stays SBP 90-200; SA1b (90-180 clamp) and SA1d (Paper 1 exclusion domain) unchanged.
- Risk categories: <3%, 3-<5%, 5-<10%, >=10% per the 2026 dyslipidemia guideline section 4.2.3.2 (verified in an earlier session; PROJECT_BRIEF). No secondary category set is added.
- Lacher 2005 (Clin Chem 51:450-452, PMID 15590751): Table 1 transcribed from the archived full text. CVw (biological, analytical removed, original scale, Tukey 3-IQR outliers removed) and method CVi: total cholesterol 8.2% / 2.3%; HDL-C 12.4% / 2.5%; serum creatinine 6.8% / 1.0%; HbA1c 1.5% / 3.1%; plasma glucose 8.3% / 1.7%; urine creatinine 43.0% / 2.2%. No urine albumin or ACR row. S7 is framed as a replication of Lacher for these analytes; comparisons use CVI = sqrt(max(CVw^2 - CVa^2, 0)) on the original scale.
- CVa table (frozen): low and high brackets are the minimum and maximum of Lacher method CVi and the NHANES III Laboratory Procedures Manual QC tables (TC 1.25-2.3%; HDL-C 2.49-2.5%; serum creatinine 1.0-2.7%; HbA1c 1.57-3.1%; glucose 1.6-2.2%; urine creatinine 2.2%). Urine albumin CVa has no NHANES-specific value; urine measures are reported descriptively only, as already specified. Exact source locations are in prefreeze_facts_v1.json.
- Selvin 2013 CVw are log-scale (via Juraschek 2012); they are compared descriptively only, as already specified.
- Pooled Cohort Equations (Stage 2 Ye 2014 reproduction): coefficients, mean coefficient x value and baseline survival for all four sex-race groups from Goff et al., Circulation 2014;129(25 Suppl 2):S49-S73, Appendix 7 Table A (identical to PMC4700825 Appendix 4 Table A), verified by two independent extractions; values and the worked example are in prefreeze_facts_v1.json.
- Still open, not blocking: user decision on public registration of this protocol hash (e.g., OSF) and target journal; verification of Ogunwale 2016, Cole 2023 and Ulbricht 2018 values before submission.
