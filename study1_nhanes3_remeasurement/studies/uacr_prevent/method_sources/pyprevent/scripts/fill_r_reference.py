#!/usr/bin/env python3
"""
Write tests/fixtures/r_reference.csv from the current Python implementation.

Use scripts/generate_r_reference.R when Rscript + AHAprevent are available for
true upstream golden values. This script locks regression coverage in CI when R
is not installed.
"""

from __future__ import annotations

import pathlib
import sys

import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from prevent import PREVENT_OUTPUT_COLUMNS, compute_prevent  # noqa: E402

CASES_PATH = ROOT / "tests" / "fixtures" / "r_cases.csv"
OUT_PATH = ROOT / "tests" / "fixtures" / "r_reference.csv"


def cases_to_prevent_df(cases: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "PAT_ID": cases["case_id"],
            "AGE": cases["age"],
            "SEX": cases["sex"],
            "TCHOL": cases["tc"],
            "HDL": cases["hdl"],
            "SBP": cases["sbp"],
            "BMI": cases["bmi"],
            "EGFR": cases["egfr"],
            "T2DM": cases["dm"],
            "RECENT_SMOKING": cases["smoking"],
            "SMOKING_CURR": cases["smoking"],
            "UACR": cases["uacr"],
            "HBA1C": cases["hba1c"],
            "BPTREAT": cases["bptreat"],
            "STATIN": cases["statin"],
            "ZIP": "00000",
        }
    )


def main() -> int:
    print(
        "WARNING: This writes Python outputs, not AHAprevent. "
        "Prefer: bash scripts/r-env/run_generate_reference.sh",
        file=sys.stderr,
    )
    cases = pd.read_csv(CASES_PATH, na_values=["", "NA"])
    df = cases_to_prevent_df(cases)
    out = compute_prevent(
        df,
        sdi_series=cases["sdi"],
        smoking_preference="SMOKING_CURR",
    )
    ref = pd.DataFrame({"case_id": cases["case_id"]})
    for col in PREVENT_OUTPUT_COLUMNS:
        ref[col] = out[col].to_numpy()
    ref.to_csv(OUT_PATH, index=False, na_rep="")
    print(f"Wrote {len(ref)} rows × {len(ref.columns)} columns to {OUT_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
