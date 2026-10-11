from __future__ import annotations

from typing import Iterable, Optional

import numpy as np
import pandas as pd

from ._base import prevent_base
from ._core import (
    REQUIRED_COLUMNS,
    CoercedRow,
    apply_age_horizon_masks,
    coerce_dataframe,
    coerce_inputs,
    coerce_zip5_series,
    to_binary01,
    to_float,
)
from ._core import adjust as adjust  # re-export
from ._core import normalize_sex as _normalize_sex  # backward-compatible test imports
from ._core import sdicat as _sdicat  # backward-compatible test imports
from ._core import sigmoid_pct as _sigmoid_pct  # backward-compatible test imports
from ._core import to_binary01 as _to_binary01  # backward-compatible test imports
from ._hba1c import prevent_hba1c
from ._full import prevent_full
from ._sdi import prevent_sdi
from ._uacr import prevent_uacr
from ._zip_sdi import lookup_sdi_decile_from_zip, lookup_sdi_deciles_series, resolve_sdi, sdi_score_to_decile

# R internal key prefix -> public column label (PREVENT{horizon}_{outcome}_{label}_PCT)
_MODEL_LABELS = {
    "prevent_base": "BASE",
    "prevent_uacr": "UACR",
    "prevent_hba1c": "HBA1C",
    "prevent_sdi": "SDI",
    "prevent_full": "FULL",
}

# Backward-compatible 10-year Basic column names (BASE -> BASIC)
_PREVENT10_BASIC_RENAME = {
    "PREVENT10_CVD_BASE_PCT": "PREVENT10_CVD_BASIC_PCT",
    "PREVENT10_ASCVD_BASE_PCT": "PREVENT10_ASCVD_BASIC_PCT",
    "PREVENT10_HF_BASE_PCT": "PREVENT10_HF_BASIC_PCT",
}

_PREVENT10_COLUMNS = [
    *_PREVENT10_BASIC_RENAME.values(),
    "PREVENT10_CVD_FULL_PCT",
    "PREVENT10_ASCVD_FULL_PCT",
    "PREVENT10_HF_FULL_PCT",
]

_PREVENT10_MODEL_PREFIXES = ("prevent_base", "prevent_full")

OPTIONAL_COLUMNS = ("BPTREAT", "STATIN")

PREVENT_OUTPUT_COLUMNS = [
    f"PREVENT{horizon}_{outcome}_{model}_PCT"
    for horizon in ("10", "30")
    for outcome in ("CVD", "ASCVD", "HF")
    for model in ("BASE", "UACR", "HBA1C", "SDI", "FULL")
]

__all__ = [
    "REQUIRED_COLUMNS",
    "OPTIONAL_COLUMNS",
    "PREVENT_OUTPUT_COLUMNS",
    "adjust",
    "compute_prevent",
    "compute_prevent10",
    "lookup_sdi_decile_from_zip",
    "lookup_sdi_deciles_series",
    "resolve_sdi",
    "sdi_score_to_decile",
    "coerce_inputs",
    "coerce_dataframe",
]


def _common_kwargs(r: CoercedRow, smoking: float, bptreat: float, statin: float) -> dict:
    return {
        "sex": r.sex,
        "age": r.age,
        "tc": r.tchol,
        "hdl": r.hdl,
        "sbp": r.sbp,
        "dm": r.t2dm,
        "smoking": smoking,
        "bmi": r.bmi,
        "egfr": r.egfr,
        "bptreat": bptreat,
        "statin": statin,
    }


def _score_row(
    r: CoercedRow,
    smoking: float,
    bptreat: float,
    statin: float,
    sdi: float,
    *,
    model_prefixes: Iterable[str] = _MODEL_LABELS,
) -> dict[str, float]:
    kw = _common_kwargs(r, smoking, bptreat, statin)
    prefixes = tuple(model_prefixes)
    merged: dict[str, float] = {}
    if "prevent_base" in prefixes:
        merged.update(prevent_base(**kw))
    if "prevent_uacr" in prefixes:
        merged.update(prevent_uacr(**kw, uacr=r.uacr))
    if "prevent_hba1c" in prefixes:
        merged.update(prevent_hba1c(**kw, hba1c=r.hba1c))
    if "prevent_sdi" in prefixes:
        merged.update(prevent_sdi(**kw, sdi=sdi))
    if "prevent_full" in prefixes:
        merged.update(prevent_full(**kw, uacr=r.uacr, hba1c=r.hba1c, sdi=sdi))
    return apply_age_horizon_masks(r.age, merged)


def _resolve_treatment(
    row_value: float,
    *,
    column_present: bool,
    call_default: float,
) -> float:
    if column_present:
        return row_value
    return call_default


def _to_public_columns(
    r_like: dict[str, float],
    *,
    model_prefixes: Iterable[str] = _MODEL_LABELS,
) -> dict[str, float]:
    prefixes = set(model_prefixes)
    out: dict[str, float] = {}
    for r_prefix, label in _MODEL_LABELS.items():
        if r_prefix not in prefixes:
            continue
        for horizon_key, horizon_num in (("10yr", "10"), ("30yr", "30")):
            for outcome in ("CVD", "ASCVD", "HF"):
                r_key = f"{r_prefix}_{horizon_key}_{outcome}"
                out[f"PREVENT{horizon_num}_{outcome}_{label}_PCT"] = r_like[r_key]
    return out


def _align_optional_series(df: pd.DataFrame, series: pd.Series, *, name: str) -> pd.Series:
    if not isinstance(series, pd.Series):
        raise TypeError(f"{name} must be a pandas Series")
    if series.index.equals(df.index):
        return series
    if len(series) == len(df):
        return pd.Series(series.to_numpy(), index=df.index, dtype=series.dtype)
    raise ValueError(
        f"{name} length ({len(series)}) does not match df length ({len(df)}); "
        "reindex to df.index or pass values in row order"
    )


def _resolve_sdi_per_row(
    zip5: str | None,
    explicit: float,
    *,
    explicit_series: bool,
    zip_lookup: float,
) -> float:
    if explicit_series:
        return explicit
    return zip_lookup


def _compute_scores(
    df: pd.DataFrame,
    *,
    bp_treat_default: Optional[int],
    statin_default: Optional[int],
    smoking_preference: str,
    sdi_series: Optional[pd.Series],
    model_prefixes: Iterable[str],
) -> pd.DataFrame:
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    if smoking_preference not in {"SMOKING_CURR", "RECENT_SMOKING"}:
        raise ValueError("smoking_preference must be 'SMOKING_CURR' or 'RECENT_SMOKING'")

    use_curr = smoking_preference == "SMOKING_CURR"
    bptreat_global = to_binary01(bp_treat_default) if bp_treat_default is not None else np.nan
    statin_global = to_binary01(statin_default) if statin_default is not None else np.nan

    has_bptreat = "BPTREAT" in df.columns
    has_statin = "STATIN" in df.columns
    explicit_sdi = sdi_series is not None
    if explicit_sdi:
        sdi_aligned = _align_optional_series(df, sdi_series, name="sdi_series")
        sdi_values = sdi_aligned.map(to_float).tolist()
    else:
        sdi_values = lookup_sdi_deciles_series(coerce_zip5_series(df["ZIP"])).tolist()

    coerced = coerce_dataframe(df)
    prefixes = tuple(model_prefixes)
    rows: list[dict[str, float]] = []
    for i, r in enumerate(coerced):
        smoking = r.smoking_curr if use_curr else r.recent_smoking
        bptreat = _resolve_treatment(r.bptreat, column_present=has_bptreat, call_default=bptreat_global)
        statin = _resolve_treatment(r.statin, column_present=has_statin, call_default=statin_global)
        if explicit_sdi:
            explicit = sdi_values[i]
            zip_lookup = np.nan
        else:
            explicit = np.nan
            zip_lookup = sdi_values[i]
        sdi = _resolve_sdi_per_row(
            r.zip5,
            explicit,
            explicit_series=explicit_sdi,
            zip_lookup=zip_lookup,
        )
        masked = _score_row(r, smoking, bptreat, statin, sdi, model_prefixes=prefixes)
        rows.append(_to_public_columns(masked, model_prefixes=prefixes))

    return pd.DataFrame(rows, index=df.index)


def compute_prevent(
    df: pd.DataFrame,
    bp_treat_default: Optional[int] = 0,
    statin_default: Optional[int] = 0,
    smoking_preference: str = "SMOKING_CURR",
    sdi_series: Optional[pd.Series] = None,
) -> pd.DataFrame:
    """
    Compute PREVENT risks for all supported models/horizons.

    Returns 30 columns: 5 models × 2 horizons × 3 outcomes.

    SDI is taken from ``sdi_series`` when provided (aligned to ``df.index``, or
  same length in row order); otherwise from ``ZIP`` via the bundled Robert
    Graham Center ZCTA crosswalk (2015–2019).

    BP treatment and statin use per-row ``BPTREAT`` / ``STATIN`` when those
    columns are present; otherwise ``bp_treat_default`` / ``statin_default``.
    """
    out = df.copy()
    scores = _compute_scores(
        out,
        bp_treat_default=bp_treat_default,
        statin_default=statin_default,
        smoking_preference=smoking_preference,
        sdi_series=sdi_series,
        model_prefixes=_MODEL_LABELS,
    )
    return pd.concat([out, scores], axis=1)


def compute_prevent10(
    df: pd.DataFrame,
    bp_treat_default: Optional[int] = 0,
    statin_default: Optional[int] = 0,
    smoking_preference: str = "SMOKING_CURR",
    sdi_series: Optional[pd.Series] = None,
) -> pd.DataFrame:
    """Backward-compatible wrapper: Base + Full 10-year risks only (12 internal equations)."""
    out = df.copy()
    scores = _compute_scores(
        out,
        bp_treat_default=bp_treat_default,
        statin_default=statin_default,
        smoking_preference=smoking_preference,
        sdi_series=sdi_series,
        model_prefixes=_PREVENT10_MODEL_PREFIXES,
    ).rename(columns=_PREVENT10_BASIC_RENAME)
    return pd.concat([out, scores[_PREVENT10_COLUMNS]], axis=1)
