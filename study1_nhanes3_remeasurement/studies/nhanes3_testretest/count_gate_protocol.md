# PREVENT test-retest under full remeasurement: count-only feasibility gate

Date frozen: 2026-10-08, before any NHANES III file was downloaded. Status: awaiting user permission to download public NHANES III files.

## Why this gate

The candidate second full paper asks how stable PREVENT-ASCVD categories are when all inputs (systolic blood pressure, total and HDL cholesterol, serum creatinine/eGFR, urine ACR, diabetes markers) are remeasured at a short-interval repeat examination, and how much each input contributes. Prior work changed one input (Ye et al. 2014, SBP only, Pooled Cohort Equations; Ulbricht 2018, SCORE) or simulated variation (Cole et al. 2023). Feasibility depends on how many NHANES III second-examination adults have complete inputs at both visits; nothing else is decided until that count exists.

## Files (public; sizes from HTTP HEAD, 2026-10-08)

exam.dat 195,090,426 B and exam.sas; lab.dat 57,945,628 B and lab.sas; adult.dat 67,107,314 B and adult.sas; examse.dat 11,609,312 B and examse.sas; labse.dat 1,804,220 B and labse.sas; pupremed.xpt (prescription medicines) with documentation. Source: https://wwwn.cdc.gov/nchs/data/nhanes3/ (1a, 2a, 3a). About 335 MB in total. Store under `studies/nhanes3_testretest/data/` with SHA256 manifest.

## Counting rules (fixed now)

- Second examination performed in the mobile examination centre (home examinations have no physician BP or urine).
- Age 30–79 at the first examination; no self-reported myocardial infarction, coronary heart disease/angina, stroke or congestive heart failure; not pregnant.
- SBP at each visit: mean of available physician readings 2 and 3, at least two readings required at each visit (sensitivity: any available physician readings).
- Complete base-model inputs at both visits: SBP, total and HDL cholesterol, serum creatinine (published NHANES III calibration applied before CKD-EPI 2021), diabetes status, smoking, antihypertensive and statin use (statin from the prescription file; status as of the first examination, carried to the second).
- UACR at both visits counted separately (LABSE UBR/URR, primary UBP/URP; note that in NHANES III URP is urinary creatinine).

## Decision

- Proceed to a frozen analysis protocol only if at least 600 participants have complete base-model inputs at both visits AND at least 400 also have UACR at both visits.
- Otherwise: narrow to base PREVENT with an SBP/lipid partition if at least 600 base pairs exist, or stop.
- Before any PREVENT scoring, reproduce Ye et al. 2014 (808 adults; 12/106 vs 2/702 crossing 7.5% with the Pooled Cohort Equations, other inputs fixed) as a pipeline check.

Only counts of eligibility and completeness are computed at this gate; no risk equation, category or change is calculated. A by-product, the NHANES III random-to-random urine arm, is held for a possible reviewer response for the UACR–PREVENT paper and not added before submission.
