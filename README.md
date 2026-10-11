# UACR–PREVENT

Recalculate paired-specimen PREVENT-ASCVD risk-category analyses from processed NHANES data.

## Run

Python 3.12.14. From this folder, using a new output directory:

```sh
python -m pip install -r requirements.txt
python -B code/reproduce.py --output ../results
```

`data/` contains seven processed input tables. `code/reproduce.py` recalculates risks, survey estimates, sensitivities and response models, saving CSV/JSON outputs. Sources and variable conventions are in `data_sources.md`.

Raw downloads, cleaning code, plotting code and figures are not included. This package supports analysis from processed inputs; it does not reconstruct raw-data preparation. Cohort and variable definitions are in the manuscript supplement. Code: MIT (see LICENSE and third_party/LICENSE-pyprevent.txt). Author-created processing and organization of data/: CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/), attributed to YunhaoJiang. Original NHANES data retain their applicable terms.

## Version 3.0.0

`study1_nhanes3_remeasurement/` adds Study 1 (NHANES III remeasurement) and the post hoc analyses E1–E12 of *Reproducibility of PREVENT Risk Categories and Treatment Recommendations on Remeasurement: Sources, Threshold Zones, and Expected Net Benefit*: scripts, protocols, inputs and a runner. See its README. The files above are unchanged.

Archived: https://doi.org/10.5281/zenodo.23299698 (version 3.0.0).
