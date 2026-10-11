# Post hoc extension E9, 2026-10-10: statin-recommendation tier under the 2026 ACC/AHA dyslipidemia guideline at two visits

Written before any E9 quantity was computed. Only data availability was checked (triglycerides TGP/TGR, LDL LCP/LCR and fasting hours
PHPFAST/PHRFAST exist in the NHANES III lab files at both visits). Post hoc, descriptive, not preregistered; no frozen estimate changes.
Under the user's evidence-first workflow, results go to the results file and the evidence ledger, not to the frozen manuscript text.

## Guideline rules used and their verification status

The primary guideline text (Circulation 2026;153:e1154-e1276; JACC doi 10.1016/j.jacc.2025.11.016) could not be retrieved automatically
(publisher sites return 403; not in PMC). Secondary summaries (ACC, Guideline Central, DynaMed, Healio, a JACC at-a-glance excerpt) agree on
the thresholds below but differ on the class of recommendation for the 5-<10% tier ('recommended' vs 'should be considered'). E9 therefore
compares tiers defined by thresholds, not by class labels:

- PREVENT-ASCVD 10-year risk tiers for adults 30-79 years without ASCVD and LDL-C 70-189 mg/dL: <3% (no drug therapy, counselling),
  3-<5% (therapy may be considered), 5-<10% (statin, at least moderate intensity), >=10% (high-intensity statin).
- Risk-independent statin indications: LDL-C >=190 mg/dL; age 40-75 years with diabetes; age 40-75 years with CKD stage 3-4 (eGFR 15-59).
- Not implemented (data unavailable or outside scope): the 30-year pathway at 30-59 years, HIV, risk enhancers (Lp(a), apoB, hsCRP, family
  history), CAC, and LDL-C goals during treatment.

If the primary text later shows a different rule, E9 is rerun with the corrected rule and both versions are kept.

## Population and inputs

- Study 1 analysis set (759) not taking a statin (statin = 0). Statin users are counted and excluded.
- At each visit: total and HDL cholesterol and triglycerides from that visit; LDL-C by the Sampson (NIH equation 2) formula,
  LDL = TC/0.948 - HDL/0.971 - (TG/8.56 + TG x nonHDL/2140 - TG^2/16100) - 9.44 (valid to TG 800 mg/dL; participants with TG >800 at either
  visit excluded and counted); base PREVENT-ASCVD 10-year risk (Study 1 visit risk); diabetes (Study 1 visit definition); eGFR (calibrated
  creatinine, 2021 CKD-EPI); age at visit 1.

## Tier at each visit (hierarchical)

1. Indicated, risk-independent (I): LDL-C >=190, or age 40-75 with diabetes, or age 40-75 with eGFR 15-59.
2. Otherwise, LDL-C <70: no risk-based pathway (N).
3. Otherwise (LDL-C 70-189), by PREVENT-ASCVD: <3% low (L), 3-<5% consider (C), 5-<10% statin (S), >=10% high-intensity statin (H).

## Estimands

- **Primary**: proportion whose statin-initiation status differs between visits, with 'statin recommended' = I, S or H (the 5% threshold
  reading of the guideline).
- **Secondary**: (a) the same with 'statin recommended or considered' = I, C, S or H (the 3% reading); (b) any change among the four
  action levels {N or L} / C / {S or I} / H; (c) high-intensity status (H) change; (d) mechanisms: changes in which the PREVENT category
  changed, LDL-C crossed 70 or 190, or the diabetes or eGFR flag changed (not mutually exclusive), and changes driven by the PREVENT category
  alone (no LDL-C, diabetes or eGFR flag change); (e) recommendation from the two-visit means (geometric means of TC, HDL, TG and creatinine;
  diabetes if fixed or mean HbA1c >=6.5%; PREVENT-ASCVD from averaged clamped inputs) compared with each single visit.
- **Uncertainty**: CR1 logit-t over the 98 pseudo-PSUs (Korn-Graubard if <10 events); visit-1 physician clustering as a sensitivity.
- **Sensitivity**: (i) NHANES-reported Friedewald LDL-C (LCP/LCR) in participants with it at both visits; (ii) Friedewald computed for all
  (TG <400); (iii) confirmed diabetes fixed across visits (SA3a rule).

## Interpretation limits

Office measurements from 1988-1994; lipid levels and treatment patterns differ from today's; secondary-source verification of the rules;
LDL-C estimated, not measured directly; risk enhancers and CAC unavailable. A tier change is a change in the guideline recommendation computed
from one visit, not a prescribing decision.

## Addendum after an independent code audit of the first run (2026-10-10; protocol text above unchanged, SHA-256 of the text before this addendum 844BE4EB...9BBB)

The audit reproduced every first-run number with independent code and found no material defect. Changes and disclosures:

1. Exclusions: the first run reported 7 exclusions under one label. They are 6 participants with TG >800 mg/dL at a visit (protocol rule)
   and 1 with triglycerides missing at visit 2 (not named in the protocol; a departure). The final run reports them separately.
2. Physician clustering: the missing physician code (88888) formed one cluster, as in E3 and E5. A version with those participants as
   singleton clusters is added.
3. Fasting hours were read but entered no estimand. A post-protocol descriptive check is added: fasting status (>=9 h) at each visit, and
   the primary estimand among participants with the same fasting status at both visits.
4. Mechanism counts are co-occurrence, not attribution. The final run adds the overlap pattern of the four flags and the number of PREVENT
   category changes in which the risk-based pathway applied at both visits (tier L, C, S or H at both).
5. The first run's JSON label was truncated; corrected. First-run outputs are preserved in ext3_statin_results_run1.json.
6. Design limits noted by the audit (not departures): the high-intensity estimand counts tier H only (the LDL-C >=190 and diabetes pathways
   are not split by intensity), so high-intensity changes are undercounted; CKD stage 3-4 is defined from one eGFR per visit, not a
   persistent reduction.
## Second addendum: rules verified against the primary guideline text; corrected-rule rerun (2026-10-10, written before any corrected-rule quantity was computed)

The full guideline text was read on jacc.org (JACC 2026;87(19); doi 10.1016/j.jacc.2025.11.016; free access). As planned above ("rerun with
the corrected rule and both versions kept"), E9 is rerun. Verified recommendations:

- Section 4.2.3.7 (adults 30-79 years, LDL-C 70-189 mg/dL), recommendation 5, COR 1, LOE A: intermediate risk (PREVENT-ASCVD 5-<10%), at
  least a moderate-intensity statin is recommended. The primary 5% reading is therefore the COR 1 reading. The secondary summaries' "should
  be considered" wording comes from the Top 10 Take-Home Messages (item 3), not from the recommendation table.
- 4.2.3.7 recommendation 4, COR 2a: borderline (3-<5%), moderate-intensity statin reasonable if a decision to start is made (tier C, as before).
- 4.2.3.7 recommendation 7, COR 1: high (>=10%), high-intensity statin recommended when lipid-lowering therapy is started (tier H, as before).
- 4.2.3.7 recommendations 2-3: low (<3%) at age 30-59, counselling if LDL-C <160 and 30-year risk <10% (COR 1); moderate-intensity statin
  reasonable (COR 2a) if LDL-C 160-189 or 30-year ASCVD risk >=10%. Not implemented in run 1.
- Adults with diabetes without ASCVD, recommendation 1, COR 1, LOE A: age 40-75 with diabetes, moderate-intensity statin indicated, with no
  LDL-C condition (as before).
- Severe hypercholesterolemia, recommendation 2, COR 1: LDL-C >=190, maximally tolerated statin (as before).
- Adults with CKD stage 3 or higher, recommendation 1, COR 1, LOE B-R: age 40-75 with CKD stage >=3 and LDL-C 70-189. Run 1 used eGFR
  15-59 without the LDL-C condition.
- Adults >75 years (COR 2b, older adults and diabetes sections) are not modelled separately; 4.2.3.7 covers ages 30-79.

Corrected tier at each visit (hierarchical):
1. I: LDL-C >=190; or age 40-75 with diabetes; or age 40-75 with eGFR <60 (one measurement per visit; stage >=3 includes <15; dialysis
   status unknown) and LDL-C >=70.
2. Otherwise LDL-C <70: N.
3. Otherwise (LDL-C 70-189), by PREVENT-ASCVD 10-year risk: >=10% H; 5-<10% S; 3-<5% C; <3%: L2 (statin reasonable, COR 2a) if age 30-59
   and (LDL-C >=160 or base PREVENT 30-year ASCVD risk >=10%), else L.

Readings: primary 'recommended' (COR 1) = I, S or H (definition unchanged); 'recommended or reasonable' (COR 1 or 2a) = I, S, H, C or L2,
which replaces the 3% reading; action levels {N, L} / {C, L2} / {S, I} / H; high intensity = H.
30-year risk: base PREVENT 30-year ASCVD from the same clamped visit inputs and fixed covariates as the 10-year risk (age at visit 1),
computed with the pinned pyprevent and independently from the pinned preventr base_30yr table; the two must agree within 1e-9 percentage
points. It is used only at ages 30-59.
All E9 estimands, intervals, sensitivity analyses and the two-visit-mean comparator are recomputed under the corrected rule, which becomes
the reference E9 result; run-1-rule results are kept. Mechanism flags add LDL-C crossing 160 and 30-year risk crossing 10% (ages 30-59).
## Third addendum after an independent audit of the rule-v2 run (2026-10-10; written after computing)

The audit reproduced every rule-v2 number and the 30-year risk (table versus pyprevent, 20,000 random inputs, maximum difference 4e-14 pp).
1. Mechanism reporting fix: the LDL-C 160 and 30-year 10% crossings cannot change the primary COR 1 status (tiers L and R are both 'not
   recommended'). The primary's mechanism table therefore keeps the run-1 flags (PREVENT category, LDL-C 70/190, diabetes, eGFR). The two new
   flags are reported for changes in the COR 1-or-2a reading, restricted to ages 30-59 with 10-year risk <3% at either visit (the
   restriction to 30-59 had not been applied to the LDL-C 160 flag).
2. Missing-value guards: the 30-year check now fails on a non-finite value, and the tier function asserts a finite 30-year risk at ages
   30-59 instead of treating a missing value as <10%. Neither guard is triggered by these data.
No estimate changes.