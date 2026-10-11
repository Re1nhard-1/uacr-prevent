from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

# Robert Graham Center SDI by ZCTA (2015–2019); see package README for citation.
_ZCTA_SDI_CSV = Path(__file__).resolve().parent / "data" / "rgc_sdi_zcta2015_2019.csv"


def sdi_score_to_decile(score: float) -> float:
    """
    Map zip-code SDI percentile score (1–100) to PREVENT decile input (1–10).

    Crosswalk from AHA PREVENT documentation (same as AHAprevent / preventr):
      1–10 -> 1, 11–20 -> 2, …, 91–100 -> 10.
    """
    if pd.isna(score):
        return np.nan
    try:
        s = int(score)
    except (TypeError, ValueError):
        return np.nan
    if s < 1 or s > 100:
        return np.nan
    return float(min(10, (s - 1) // 10 + 1))


@lru_cache(maxsize=1)
def _zcta_decile_table() -> dict[str, float]:
    df = pd.read_csv(_ZCTA_SDI_CSV, usecols=["ZCTA5_FIPS", "SDI_score"], dtype={"ZCTA5_FIPS": str})
    zcta = df["ZCTA5_FIPS"].str.zfill(5)
    deciles = df["SDI_score"].map(sdi_score_to_decile)
    return dict(zip(zcta.tolist(), deciles.tolist()))


def lookup_sdi_decile_from_zip(zip5: str | None) -> float:
    """Return PREVENT SDI decile (1–10) for a 5-digit ZCTA, or NaN if unknown."""
    if zip5 is None:
        return np.nan
    key = str(zip5).zfill(5)[:5]
    return _zcta_decile_table().get(key, np.nan)


def lookup_sdi_deciles_series(zip5: pd.Series) -> pd.Series:
    """Map a series of 5-digit ZCTA strings to PREVENT SDI deciles (1–10)."""
    table = _zcta_decile_table()

    def _one(z) -> float:
        if z is None or (isinstance(z, float) and np.isnan(z)):
            return np.nan
        return table.get(str(z).zfill(5)[:5], np.nan)

    return zip5.map(_one)


def resolve_sdi(zip5: str | None, explicit: float, *, explicit_series: bool) -> float:
    """
    SDI resolution order:
      1. ``sdi_series`` when passed to ``compute_prevent*`` (may be NaN)
      2. ZIP → RGC ZCTA crosswalk
      3. NaN (missing-SDI coefficients in SDI / Full models)
    """
    if explicit_series:
        return explicit
    return lookup_sdi_decile_from_zip(zip5)
