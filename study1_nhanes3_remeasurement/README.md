# Study 1: PREVENT risk categories on remeasurement (NHANES III)

Code and inputs for Study 1 of *Reproducibility of PREVENT Risk Categories and Treatment Recommendations on Remeasurement: Sources, Threshold Zones, and Expected Net Benefit* (Yunhao Jiang) and its post hoc analyses E1–E12. Study 2 (NHANES 2009–2010 urine
specimens) is the version 2.0.0 package in the parent folder; its numbers enter here through `studies/nhanes3_testretest/manuscript_evidence_v1/study2_evidence.json`.

## Run

```powershell
python -m pip install -r requirements.txt     # Python 3.12 (the frozen run: 3.12.14, Windows 11)
python download_nhanes3.py                    # 17 NHANES III files, 335 MB, SHA-256 checked against the original download
.\run_all.ps1 -SkipLong                       # about 10 minutes
```

Run in a fresh copy of this folder: the scripts refuse to overwrite a finished run. `run_all.ps1` runs the scripts from `scripts/`
in order and writes everything under `studies/nhanes3_testretest/`. Without `-SkipLong` it also recomputes E11, the uniform-start
κ interval and E12 (about 6 hours with 12 processes; `Pool(12)` is hard-coded); with it, their shipped results in
`analysis_v1_posthoc/` are used. `-SkipDocs` leaves out the Word builds. Other shells: run the same scripts in the order listed in
`run_all.ps1` with `PYTHONUTF8=1` and copy `records/*` into `analysis_v1/` after the first script.

## Outputs

| Folder | Content |
|---|---|
| `analysis_v1/` | Stage 1–4 under the frozen protocol: readers and counts, first PREVENT scoring (primary estimand 79/759 = 10.4%, 95% CI 8.2–13.1%), independent recomputation, sensitivity analyses |
| `analysis_v1_posthoc/`, `analysis_v1_corrections/`, `manuscript_evidence_v1/` | post hoc analyses and E1–E12, the S9 correction, Table 1 |
| `figures_v1/`, `figures_v7/`, `submission_ajpc/figures/` | figure panels; Figures 1–4 of the paper; submission exports and the Central Illustration |
| `EVIDENCE_LEDGER.md` | every number the paper reports, with the file it is read from |
| `manuscript_v7/`, `submission_ajpc/manuscript/` | the manuscript and Supplementary Material (Word), built from the result files |

## Contents

- `scripts/`: the analysis, figure and document scripts of the run that produced the paper, with comments and docstrings removed (the code is otherwise identical); `vendor/figure_tools.py` exports the figures.
- `studies/nhanes3_testretest/`: the frozen protocol (`analysis_protocol_v1.md`, hashes in `freeze_record_v1.json`, checked by Stage 1), the dated post hoc protocols with their addenda, the count-gate result, the NHANES III download manifest, the NHANES 2001–2002 and 2021–2023 files with manifests, `records/` (two files written once and annotated by the authors), the Study 2 numbers and the bibliography data, and the results of the three long steps.
- `studies/uacr_prevent/`: pinned `preventr` coefficient tables and `pyprevent` sources (MIT; `pyprevent` is imported from here, not from PyPI), the Study 2 participant tables (same values as the v2.0.0 `data/` files on every shared column), `DEMO_F.xpt`, and the NHANES 2001–2002 first-exam files used by E1.

## Notes

- `nh3_stats.py` differs from the copy that ran Stage 3 by one line (the Tango score-interval variance, corrected 2026-10-09 after the results); no frozen number depends on it.
- Result files carry a completion timestamp; everything else reproduced byte for byte on the frozen environment (checked 2026-10-11). Last floating-point digits can differ across platforms.
- The protocol was not registered publicly; the freeze record and its hash checks are the evidence of prespecification.

## Licenses

Author's code: MIT (`LICENSE`). Author-created data tables: CC BY 4.0. `preventr` and `pyprevent`: MIT (their license files are in `studies/uacr_prevent/method_sources/`). NHANES public-use files keep NCHS terms.
