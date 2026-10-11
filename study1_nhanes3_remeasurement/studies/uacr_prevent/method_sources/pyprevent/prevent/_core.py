from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import pandas as pd


REQUIRED_COLUMNS = [
    "PAT_ID",
    "AGE",
    "SEX",
    "TCHOL",
    "HDL",
    "SBP",
    "BMI",
    "EGFR",
    "T2DM",
    "RECENT_SMOKING",
    "SMOKING_CURR",
    "UACR",
    "HBA1C",
    "ZIP",
]


def mmol_conversion(cholesterol_mgdl: float) -> float:
    return 0.02586 * cholesterol_mgdl


def sdicat(sdi: float) -> float:
    if pd.isna(sdi):
        return np.nan
    if 0 < sdi < 4:
        return 0.0
    if 4 <= sdi < 7:
        return 1.0
    if 7 <= sdi <= 10:
        return 2.0
    return np.nan


def sigmoid_pct(x: float) -> float:
    return 100.0 / (1.0 + math.exp(-x))


def adjust(uacr: float):
    # Mirrors AHAprevent::adjust(): floor UACR in [0, 0.1) to 0.1 before log().
    if pd.isna(uacr):
        return None
    if uacr < 0:
        return None
    if uacr >= 0.1:
        return uacr
    if 0 <= uacr < 0.1:
        return 0.1
    return None


def invalid_uacr(uacr: float) -> bool:
    """R parity: pred_risk_* post-check ``!is.na(uacr) & uacr < 0``."""
    return not pd.isna(uacr) and uacr < 0


def invalid_hba1c(hba1c: float) -> bool:
    """R parity: pred_risk_* post-check ``!is.na(hba1c) & hba1c <= 0``."""
    return not pd.isna(hba1c) and hba1c <= 0


def invalid_sdi_decile(sdi: float) -> bool:
    """R parity: SDI must be missing or an integer decile in [1, 10]."""
    if pd.isna(sdi):
        return False
    try:
        decile = float(sdi)
    except (TypeError, ValueError):
        return True
    if decile != int(decile):
        return True
    return decile < 1 or decile > 10


def to_float(x) -> float:
    try:
        if pd.isna(x):
            return np.nan
        return float(x)
    except Exception:
        return np.nan


def to_binary01(x) -> float:
    try:
        if pd.isna(x):
            return np.nan
        xv = int(x)
        return float(xv) if xv in (0, 1) else np.nan
    except Exception:
        return np.nan


def to_binary01_series(values: pd.Series) -> pd.Series:
    return values.map(to_binary01)


def to_float_series(values: pd.Series) -> pd.Series:
    return values.map(to_float)


def normalize_sex(x) -> float:
    """
    PREVENT source expects:
      0 = male
      1 = female
    """
    if pd.isna(x):
        return np.nan
    s = str(x).strip().lower()
    if s in {"m", "male", "0"}:
        return 0.0
    if s in {"f", "female", "1"}:
        return 1.0
    try:
        xv = float(x)
        if xv in (0.0, 1.0):
            return xv
    except Exception:
        pass
    return np.nan


@dataclass(frozen=True)
class CoercedRow:
    pat_id: object
    age: float
    sex: float
    tchol: float
    hdl: float
    sbp: float
    bmi: float
    egfr: float
    t2dm: float
    recent_smoking: float
    smoking_curr: float
    uacr: float
    hba1c: float
    bptreat: float
    statin: float
    zip5: str | None


def coerce_zip5(zip_raw) -> str | None:
    if pd.isna(zip_raw):
        return None
    digits = "".join(c for c in str(zip_raw).strip() if c.isdigit())[:5]
    return digits.zfill(5) if digits else None


def coerce_zip5_series(zips: pd.Series) -> pd.Series:
    """Vectorized 5-digit ZCTA strings for SDI lookup (nullable)."""
    return zips.map(coerce_zip5)


def coerce_dataframe(df: pd.DataFrame) -> list[CoercedRow]:
    """Coerce an entire input frame to ``CoercedRow`` list (column-wise maps, one pass)."""
    n = len(df)
    has_bptreat = "BPTREAT" in df.columns
    has_statin = "STATIN" in df.columns
    zip5 = coerce_zip5_series(df["ZIP"]).tolist()
    ages = to_float_series(df["AGE"]).tolist()
    sexes = df["SEX"].map(normalize_sex).tolist()
    tchol = to_float_series(df["TCHOL"]).tolist()
    hdl = to_float_series(df["HDL"]).tolist()
    sbp = to_float_series(df["SBP"]).tolist()
    bmi = to_float_series(df["BMI"]).tolist()
    egfr = to_float_series(df["EGFR"]).tolist()
    t2dm = to_binary01_series(df["T2DM"]).tolist()
    recent_smoking = to_binary01_series(df["RECENT_SMOKING"]).tolist()
    smoking_curr = to_binary01_series(df["SMOKING_CURR"]).tolist()
    uacr = to_float_series(df["UACR"]).tolist()
    hba1c = to_float_series(df["HBA1C"]).tolist()
    bptreat = to_binary01_series(df["BPTREAT"]).tolist() if has_bptreat else [np.nan] * n
    statin = to_binary01_series(df["STATIN"]).tolist() if has_statin else [np.nan] * n
    pat_ids = df["PAT_ID"].tolist()
    return [
        CoercedRow(
            pat_id=pat_ids[i],
            age=ages[i],
            sex=sexes[i],
            tchol=tchol[i],
            hdl=hdl[i],
            sbp=sbp[i],
            bmi=bmi[i],
            egfr=egfr[i],
            t2dm=t2dm[i],
            recent_smoking=recent_smoking[i],
            smoking_curr=smoking_curr[i],
            uacr=uacr[i],
            hba1c=hba1c[i],
            bptreat=bptreat[i],
            statin=statin[i],
            zip5=zip5[i],
        )
        for i in range(n)
    ]


def coerce_inputs(row: pd.Series) -> CoercedRow:
    """
    Coerce types and normalize coding but do NOT clip ranges.

    Range handling must match the upstream R functions: out-of-range inputs
    lead to NA (NaN) risk outputs rather than silent clipping.
    """
    zip5 = coerce_zip5(row.get("ZIP"))
    return CoercedRow(
        pat_id=row.get("PAT_ID"),
        age=to_float(row.get("AGE")),
        sex=normalize_sex(row.get("SEX")),
        tchol=to_float(row.get("TCHOL")),
        hdl=to_float(row.get("HDL")),
        sbp=to_float(row.get("SBP")),
        bmi=to_float(row.get("BMI")),
        egfr=to_float(row.get("EGFR")),
        t2dm=to_binary01(row.get("T2DM")),
        recent_smoking=to_binary01(row.get("RECENT_SMOKING")),
        smoking_curr=to_binary01(row.get("SMOKING_CURR")),
        uacr=to_float(row.get("UACR")),
        hba1c=to_float(row.get("HBA1C")),
        bptreat=to_binary01(row.get("BPTREAT")) if "BPTREAT" in row.index else np.nan,
        statin=to_binary01(row.get("STATIN")) if "STATIN" in row.index else np.nan,
        zip5=zip5,
    )


def validate_common_inputs(age, sex, sbp, dm, smoking, egfr) -> bool:
    # Mirrors AHAprevent required-variable behavior (no clipping):
    # missing or out-of-range required inputs => all outcomes NA.
    if pd.isna(age) or age < 30 or age > 79:
        return False
    if pd.isna(sex) or sex not in (0.0, 1.0):
        return False
    if pd.isna(sbp) or sbp < 90 or sbp > 200:
        return False
    if pd.isna(dm) or dm not in (0.0, 1.0):
        return False
    if pd.isna(smoking) or smoking not in (0.0, 1.0):
        return False
    if pd.isna(egfr) or egfr <= 0:
        return False
    return True


def apply_age_horizon_masks(age: float, results: dict[str, float]) -> dict[str, float]:
    """
    Apply the R post-processing age rules (per model).

    - If 59 < age <= 79: all 30-year risks are NA
    - If age < 30 or age > 79: all risks (10-year and 30-year) are NA
    """
    if pd.isna(age):
        return results

    out = dict(results)

    if 59 < age <= 79:
        for k in list(out.keys()):
            if "_30yr_" in k:
                out[k] = np.nan
    elif age < 30 or age > 79:
        for k in list(out.keys()):
            out[k] = np.nan

    return out

