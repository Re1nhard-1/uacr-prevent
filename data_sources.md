# Sources and inputs

[NHANES public-use data](https://wwwn.cdc.gov/nchs/nhanes/Default.aspx): 2009–2010 random/first-morning urine pairs and 2001–2002 selected repeat-examination pairs. Original component links are below. SEQN is the public survey identifier. Cohort selection and derived-variable rules are in the manuscript supplement.

| Processed CSV | Contents |
|---|---|
| primary_participant_results | 3,166 pairs, equation inputs, risks/categories, survey design and sensitivity-subset fields |
| external_bridged_participant_results | 278 pairs, bridged UACRs, equation inputs, risks/categories and pregnancy flag |
| external_native | Same 278 pairs with native assay measurements |
| alternative_bp | 3,160 participants with the alternative blood-pressure inputs |
| response_target | 3,441 eligible first specimens; 275 missing usable repeats |
| primary_design / external_design | 31 / 30 original positive-weight stratum/PSU units; placeholders retain zero-domain variance contributions, not additional participants |

UACR: mg/g; cholesterol/LDL: mg/dL; SBP: mmHg; BMI: kg/m²; eGFR: mL/min/1.73 m². Risk is percent. Derived sex: male=0, female=1; diabetes/smoking/treatment indicators: 0/1. Categories 0–3: <3%, 3–<5%, 5–<10%, ≥10%, using unrounded risks. Primary weights: WTMEC2YR; fasting subset: WTSAF2YR. External unit weights describe selected repeaters, not national prevalence.

Calculation: [pyprevent 0.2.0, commit 69f341472daf9cb049a399b3e38db60fd8896e56](https://github.com/kingrc15/pyprevent/tree/69f341472daf9cb049a399b3e38db60fd8896e56), UACR/core modules, MIT license in `third_party/`. Comments/docstrings were removed and the package entry point narrowed; calculation logic is unchanged. Equation: [Khan et al., Circulation 2024](https://doi.org/10.1161/CIRCULATIONAHA.123.067626). External inputs use the [CDC urine-creatinine bridge](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2007/DataFiles/ALB_CR_E.htm); bridge-coefficient uncertainty is not propagated.

## Original public files

### NHANES 2009–2010

- [DEMO_F](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/DEMO_F.htm)
- [ALB_CR_F](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/ALB_CR_F.htm)
- [BIOPRO_F](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/BIOPRO_F.htm)
- [TCHOL_F](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/TCHOL_F.htm)
- [HDL_F](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/HDL_F.htm)
- [BMX_F](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/BMX_F.htm)
- [BPX_F](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/BPX_F.htm)
- [BPQ_F](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/BPQ_F.htm)
- [DIQ_F](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/DIQ_F.htm)
- [GHB_F](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/GHB_F.htm)
- [SMQ_F](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/SMQ_F.htm)
- [MCQ_F](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/MCQ_F.htm)
- [RXQ_RX_F](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/RXQ_RX_F.htm)
- [TRIGLY_F](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/TRIGLY_F.htm)
### NHANES 2001–2002

- [DEMO_B](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/DEMO_B.htm)
- [L16_B](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/L16_B.htm)
- [L16_2_B](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/L16_2_B.htm)
- [L40_B](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/L40_B.htm)
- [L13_B](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/L13_B.htm)
- [BPX_B](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/BPX_B.htm)
- [BPQ_B](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/BPQ_B.htm)
- [DIQ_B](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/DIQ_B.htm)
- [L10_B](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/L10_B.htm)
- [SMQ_B](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/SMQ_B.htm)
- [MCQ_B](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/MCQ_B.htm)
- [RXQ_RX_B](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/RXQ_RX_B.htm)
- [BMX_B](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/BMX_B.htm)
