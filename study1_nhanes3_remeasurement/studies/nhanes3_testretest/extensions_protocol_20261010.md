# Post hoc extensions, 2026-10-10 (written before running either analysis)

Both analyses were specified after all Study 1 and Study 2 results were known. They are post hoc, descriptive and not preregistered.
No Study 1 or Study 2 estimate is changed.

## E1. Laboratory-only remeasurement in NHANES 2001-2002 (Study 2 external sample)

- **Sample**: the 278 participants of the Study 2 external analysis (repeat MEC examination, aged 30-69 years), restricted to those with
  total and HDL cholesterol (L13_2_B), serum creatinine (L40_2_B) and HbA1c (L10_2_B) at the second examination. Flow is reported.
- **Inputs**: visit 1 values as in the Study 2 participant file (blood pressure, age, sex, smoking, treatment and statin use fixed);
  visit 2: LB2TC, LB2HDL, LB2SCR (eGFR, 2021 CKD-EPI, first-visit age), LB2GH (diabetes = questionnaire-positive [DIQ010, DIQ050 or DIQ070
  equal to 1] or HbA1c >= 6.5%), and the bridged repeat UACR for the UACR model. Visit 1 inputs are recomputed from the source files and must
  reproduce the Study 2 participant file (TC, HDL, eGFR, diabetes, first-urine UACR-model risk). Inputs are clamped to the PREVENT domain
  as in Study 1 (no participant is expected to need it).
- **Primary estimand**: proportion whose base PREVENT-ASCVD category (<3, 3-<5, 5-<10, >=10%) changes when the laboratory inputs are
  remeasured with blood pressure fixed (theta_lab,2001). Secondary: the same for the UACR model with laboratory inputs and urine
  remeasured; up/down counts; within-person CVw (SD of the log difference / sqrt 2) of TC, HDL, creatinine; SD of the HbA1c difference / sqrt 2;
  interval in days.
- **Comparator (descriptive)**: Study 1 laboratory-only swap with visit 1 as reference (base hybrid mask 14) and laboratory+urine swap
  (UACR hybrid mask 30), from the frozen participant export.
- **Uncertainty**: equal weights; stratified Taylor variance over the masked strata and PSUs; logit-scale t interval with df = PSUs - strata
  (as in the Study 2 external analysis). With fewer than 10 events, Korn-Graubard intervals.
- **Interpretation**: a different era, assays and population; no blood pressure remeasurement; small numbers. Described, not tested.

## E2. Contemporary US adults: NHANES August 2021-August 2023

- **Population**: examined participants aged 30-79 years with positive MEC weight; each of MCQ160B-F explicitly 'no'; not pregnant
  (RIDEXPRG != 1); complete inputs.
- **Inputs**: SBP and DBP = mean of oscillometric readings 2 and 3 (else mean of available readings); total cholesterol (LBXTC), HDL
  (LBDHDD), creatinine (LBXSCR; 2021 CKD-EPI), diabetes (DIQ010/DIQ050/DIQ070 = 1 or HbA1c >= 6.5%; negative requires observed HbA1c),
  current smoking (SMQ020 = 1 and SMQ040 in 1-2), antihypertensive treatment (BPQ150 = 1; not asked if BPQ020 = 2, coded 0), statin
  proxy = taking prescribed cholesterol-lowering medication (BPQ101D = 1; drug-level prescription data are not released for this cycle).
  Clamped to the PREVENT domain as in Study 1.
- **Estimands** (MEC-weighted, Taylor linearization over SDMVSTRA/SDMVPSU, df = PSUs - strata):
  1. Distribution of base PREVENT-ASCVD categories; weighted number of adults.
  2. Proportion and number with first-measurement risk inside the frozen base pi_clin instability zones (>=25% and >=10%, truncated).
  3. Transported expected proportion changing category on remeasurement = weighted mean of pi_clin(risk), interval combining the
     4,000 refitted curve replicates (population fixed) and the survey variance (curve fixed), assumed independent.
  4. Untreated stage 1 hypertension (no antihypertensive treatment; SBP 130-139 or DBP 80-89 with SBP < 140 and DBP < 90): number,
     proportion with PREVENT 10-year total CVD >= 7.5%, proportion inside the >=10% zone of the 7.5% crossing curve (post hoc curve,
     cvd75_curve.json), and the transported expected proportion crossing 7.5%.
- **Assumptions to state**: 1988-1994 within-person variation, physician mercury blood pressure and research protocols apply to 2021-2023
  oscillometric measurements and laboratory methods; change probability depends only on the distance from a threshold; the statin proxy
  includes non-statin lipid-lowering drugs. Exploratory; supplementary only, with at most one sentence in the main text.

## Departures and corrections recorded after the first run (2026-10-10; review wf_70e4f249-73d)

The text above is unchanged. After the first outputs were seen:

1. E1 intervals: the first run floored the Taylor variance at the binomial variance, which the protocol did not specify; the floor was removed
   (base 4.4%, 2.7-7.1%; UACR model 10.3%, 7.2-14.6%).
2. E1 clamping: contrary to the expectation above, 9 participants had a visit-2 input outside the PREVENT domain (8 total cholesterol
   <130 mg/dL, 1 HDL >100 mg/dL; none at visit 1). Checks: excluding them, 12/262 (4.6%); without clamping, 13/271 (4.8%).
3. E2 weights: blood analytes come from the phlebotomy subsample, so the phlebotomy weight (WTPH2YR) replaced the MEC weight named above.
4. E2 population numbers: complete-case weighted totals are additionally reported scaled by the ratio of eligible to complete-case weights,
   which assumes missing inputs are unrelated to risk; this scaling was added after the first output.
5. E2 expected crossing of 7.5%: the interval now combines the refitted 4,000 curve replicates with the survey variance.
6. Correction of a factual error above: drug-level prescription data were not collected in the 2021-2023 cycle (not merely unreleased).
7. Reporting: the first draft placed E2 in the abstract, Results and Discussion, beyond the 'at most one sentence in the main text' rule above;
   this was reverted to one Results sentence, with the rest in Supplementary Table S30.

## Reporting note added with the condensed manuscript (revision 5, 2026-10-10)

8. In the condensed main text, E2 has one Results sentence; the Methods only name it within the sentence that introduces both extensions and
   point to the supplement. Its methods and the assumptions to state (including that the statin proxy includes non-statin lipid-lowering
   drugs) are in the Supplementary Methods and the Supplementary Table S30 note. No estimate changed.
