from __future__ import annotations

import numpy as np
import pandas as pd

from ._core import invalid_sdi_decile, mmol_conversion, sdicat, sigmoid_pct, validate_common_inputs


def _sdi_nan() -> dict[str, float]:
    return {
        "prevent_sdi_10yr_CVD": np.nan,
        "prevent_sdi_10yr_ASCVD": np.nan,
        "prevent_sdi_10yr_HF": np.nan,
        "prevent_sdi_30yr_CVD": np.nan,
        "prevent_sdi_30yr_ASCVD": np.nan,
        "prevent_sdi_30yr_HF": np.nan,
    }


def prevent_sdi(sex, age, tc, hdl, sbp, dm, smoking, bmi, egfr, bptreat, statin, sdi) -> dict[str, float]:
    """
    R parity: AHAprevent::pred_risk_sdi (10yr + 30yr).
    """
    if not validate_common_inputs(age, sex, sbp, dm, smoking, egfr):
        return _sdi_nan()
    if invalid_sdi_decile(sdi):
        return _sdi_nan()

    can_cvd_ascvd = not (
        pd.isna(tc)
        or tc < 130
        or tc > 320
        or pd.isna(hdl)
        or hdl < 20
        or hdl > 100
        or pd.isna(statin)
        or statin not in (0.0, 1.0)
        or pd.isna(bptreat)
        or bptreat not in (0.0, 1.0)
    )
    can_hf = not (
        pd.isna(bmi)
        or bmi < 18.5
        or bmi >= 40
        or pd.isna(bptreat)
        or bptreat not in (0.0, 1.0)
    )

    def sdi_contrast(a, b):
        if not pd.isna(sdi):
            s = sdicat(sdi)
            return a * ((2 - s) * s) + b * ((s - 1) * (0.5 * s))
        return np.nan

    logor_10yr_CVD = np.nan
    logor_10yr_ASCVD = np.nan
    logor_10yr_HF = np.nan
    logor_30yr_CVD = np.nan
    logor_30yr_ASCVD = np.nan
    logor_30yr_HF = np.nan

    if sex == 1.0:  # female
        if can_cvd_ascvd:
            logor_10yr_CVD = (
                -3.461564
                + 0.7754083 * ((age - 55) / 10)
                + 0.0221756 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.1650828 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.2180808 * (min(sbp, 110) - 110) / 20
                + 0.3381188 * (max(sbp, 110) - 130) / 20
                + 0.8624372 * dm
                + 0.4663953 * smoking
                + 0.5919004 * (min(egfr, 60) - 60) / (-15)
                + 0.0516821 * (max(egfr, 60) - 90) / (-15)
                + 0.3182166 * bptreat
                - 0.1460816 * statin
                - 0.0574455 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.1302287 * statin * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.083509 * (age - 55) / 10 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                + 0.0282181 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.0952647 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.2718966 * (age - 55) / 10 * dm
                - 0.0641738 * (age - 55) / 10 * smoking
                - 0.1717026 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (sdi_contrast(0.1442776, 0.2421409) if not pd.isna(sdi) else 0.1885076)
            )

            logor_10yr_ASCVD = (
                -3.955898
                + 0.7028123 * ((age - 55) / 10)
                + 0.1056078 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.1502263 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.0488757 * (min(sbp, 110) - 110) / 20
                + 0.3402681 * (max(sbp, 110) - 130) / 20
                + 0.838022 * dm
                + 0.4064592 * smoking
                + 0.4838394 * (min(egfr, 60) - 60) / (-15)
                + 0.0480415 * (max(egfr, 60) - 90) / (-15)
                + 0.2270648 * bptreat
                - 0.0585626 * statin
                - 0.0349485 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.1017299 * statin * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.062389 * (age - 55) / 10 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                + 0.0285106 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.1033711 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.2477845 * (age - 55) / 10 * dm
                - 0.0544326 * (age - 55) / 10 * smoking
                - 0.1735372 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (sdi_contrast(0.1473705, 0.2451878) if not pd.isna(sdi) else 0.1691593)
            )

            logor_30yr_CVD = (
                -1.493211
                + 0.5124233 * ((age - 55) / 10)
                - 0.0978159 * (((age - 55) / 10) ** 2)
                + 0.0322131 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.1717884 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.1364536 * (min(sbp, 110) - 110) / 20
                + 0.3074443 * (max(sbp, 110) - 130) / 20
                + 0.6709275 * dm
                + 0.2897728 * smoking
                + 0.1670658 * (min(egfr, 60) - 60) / (-15)
                + 0.0618439 * (max(egfr, 60) - 90) / (-15)
                + 0.2969806 * bptreat
                - 0.0665514 * statin
                - 0.0458917 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.1168505 * statin * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0770419 * (age - 55) / 10 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                + 0.027634 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.0992045 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.3208137 * (age - 55) / 10 * dm
                - 0.134847 * (age - 55) / 10 * smoking
                - 0.1399842 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (sdi_contrast(0.1129725, 0.1975843) if not pd.isna(sdi) else 0.1627381)
            )

            logor_30yr_ASCVD = (
                -2.116951
                + 0.4396545 * ((age - 55) / 10)
                - 0.0918489 * (((age - 55) / 10) ** 2)
                + 0.1132729 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.1544977 * (mmol_conversion(hdl) - 1.3) / 0.3
                + 0.036315 * (min(sbp, 110) - 110) / 20
                + 0.3049229 * (max(sbp, 110) - 130) / 20
                + 0.6344794 * dm
                + 0.234514 * smoking
                + 0.0898312 * (min(egfr, 60) - 60) / (-15)
                + 0.0564502 * (max(egfr, 60) - 90) / (-15)
                + 0.1933487 * bptreat
                + 0.0220467 * statin
                - 0.0229229 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.0903326 * statin * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0579383 * (age - 55) / 10 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                + 0.0274011 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.1039749 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.2785102 * (age - 55) / 10 * dm
                - 0.1167267 * (age - 55) / 10 * smoking
                - 0.1269382 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (sdi_contrast(0.1149139, 0.1976537) if not pd.isna(sdi) else 0.1391241)
            )

        if can_hf:
            logor_10yr_HF = (
                -4.409382
                + 0.8819156 * ((age - 55) / 10)
                - 0.4495491 * (min(sbp, 110) - 110) / 20
                + 0.3457405 * (max(sbp, 110) - 130) / 20
                + 1.02632 * dm
                + 0.5371646 * smoking
                - 0.0168447 * (min(bmi, 30) - 25) / 5
                + 0.2805126 * (max(bmi, 30) - 30) / 5
                + 0.7315223 * (min(egfr, 60) - 60) / (-15)
                + 0.0651679 * (max(egfr, 60) - 90) / (-15)
                + 0.3491487 * bptreat
                - 0.0890335 * bptreat * (max(sbp, 110) - 130) / 20
                - 0.0971028 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.3528078 * (age - 55) / 10 * dm
                - 0.106216 * (age - 55) / 10 * smoking
                + 0.0064998 * (age - 55) / 10 * (max(bmi, 30) - 30) / 5
                - 0.1899413 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (sdi_contrast(0.1343318, 0.2496522) if not pd.isna(sdi) else 0.1915023)
            )

            logor_30yr_HF = (
                -2.317899
                + 0.5919097 * ((age - 55) / 10)
                - 0.1023133 * (((age - 55) / 10) ** 2)
                - 0.3864727 * (min(sbp, 110) - 110) / 20
                + 0.301876 * (max(sbp, 110) - 130) / 20
                + 0.8162909 * dm
                + 0.3449647 * smoking
                + 0.0574975 * (min(bmi, 30) - 25) / 5
                + 0.2367826 * (max(bmi, 30) - 30) / 5
                + 0.2790347 * (min(egfr, 60) - 60) / (-15)
                + 0.0742645 * (max(egfr, 60) - 90) / (-15)
                + 0.3352935 * bptreat
                - 0.0772532 * bptreat * (max(sbp, 110) - 130) / 20
                - 0.0995144 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.4000423 * (age - 55) / 10 * dm
                - 0.1770335 * (age - 55) / 10 * smoking
                + 0.0083046 * (age - 55) / 10 * (max(bmi, 30) - 30) / 5
                - 0.149585 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (sdi_contrast(0.0960646, 0.1987543) if not pd.isna(sdi) else 0.1562214)
            )
    else:  # male
        if can_cvd_ascvd:
            logor_10yr_CVD = (
                -3.159572
                + 0.7756377 * ((age - 55) / 10)
                + 0.0715325 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0976775 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.5186614 * (min(sbp, 110) - 110) / 20
                + 0.3235653 * (max(sbp, 110) - 130) / 20
                + 0.7722496 * dm
                + 0.3761129 * smoking
                + 0.5180893 * (min(egfr, 60) - 60) / (-15)
                + 0.0118451 * (max(egfr, 60) - 90) / (-15)
                + 0.2634094 * bptreat
                - 0.1455263 * statin
                - 0.0367013 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.1617785 * statin * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0507669 * (age - 55) / 10 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                + 0.0178356 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.1059337 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.2236755 * (age - 55) / 10 * dm
                - 0.0723216 * (age - 55) / 10 * smoking
                - 0.1548205 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (sdi_contrast(0.0889119, 0.291897) if not pd.isna(sdi) else 0.1508151)
            )

            logor_10yr_ASCVD = (
                -3.624712
                + 0.7150087 * ((age - 55) / 10)
                + 0.1627339 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.1194988 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.363659 * (min(sbp, 110) - 110) / 20
                + 0.3179476 * (max(sbp, 110) - 130) / 20
                + 0.7156422 * dm
                + 0.3404477 * smoking
                + 0.3545754 * (min(egfr, 60) - 60) / (-15)
                + 0.0157875 * (max(egfr, 60) - 90) / (-15)
                + 0.1786233 * bptreat
                - 0.1018269 * statin
                - 0.028313 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.1209467 * statin * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0285806 * (age - 55) / 10 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                + 0.0247348 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.0919494 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.1981491 * (age - 55) / 10 * dm
                - 0.0776891 * (age - 55) / 10 * smoking
                - 0.1284899 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (sdi_contrast(0.0728242, 0.2824453) if not pd.isna(sdi) else 0.1437348)
            )

            logor_30yr_CVD = (
                -1.251031
                + 0.437377 * ((age - 55) / 10)
                - 0.104443 * (((age - 55) / 10) ** 2)
                + 0.0812573 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.1069199 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.2786727 * (min(sbp, 110) - 110) / 20
                + 0.2729256 * (max(sbp, 110) - 130) / 20
                + 0.5279006 * dm
                + 0.1878949 * smoking
                + 0.0866569 * (min(egfr, 60) - 60) / (-15)
                + 0.0594948 * (max(egfr, 60) - 90) / (-15)
                + 0.2028246 * bptreat
                - 0.0308404 * statin
                - 0.0283679 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.1439353 * statin * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0510854 * (age - 55) / 10 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                + 0.0150236 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.1095448 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.2561109 * (age - 55) / 10 * dm
                - 0.1282945 * (age - 55) / 10 * smoking
                - 0.1011023 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (sdi_contrast(0.0314626, 0.2003953) if not pd.isna(sdi) else 0.0927451)
            )

            logor_30yr_ASCVD = (
                -1.836632
                + 0.3749788 * ((age - 55) / 10)
                - 0.0990063 * (((age - 55) / 10) ** 2)
                + 0.1708505 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.1272841 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.1275555 * (min(sbp, 110) - 110) / 20
                + 0.2659339 * (max(sbp, 110) - 130) / 20
                + 0.4676531 * dm
                + 0.1610104 * smoking
                - 0.0465144 * (min(egfr, 60) - 60) / (-15)
                + 0.0596996 * (max(egfr, 60) - 90) / (-15)
                + 0.1147096 * bptreat
                + 0.0052906 * statin
                - 0.0186687 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.1063151 * statin * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0307797 * (age - 55) / 10 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                + 0.0218126 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.0898242 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.210054 * (age - 55) / 10 * dm
                - 0.1246327 * (age - 55) / 10 * smoking
                - 0.0629358 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (sdi_contrast(0.0199201, 0.194949) if not pd.isna(sdi) else 0.0863835)
            )

        if can_hf:
            logor_10yr_HF = (
                -4.058977
                + 0.894179 * ((age - 55) / 10)
                - 0.7067398 * (min(sbp, 110) - 110) / 20
                + 0.350241 * (max(sbp, 110) - 130) / 20
                + 0.9252453 * dm
                + 0.4364765 * smoking
                - 0.0866297 * (min(bmi, 30) - 25) / 5
                + 0.3706765 * (max(bmi, 30) - 30) / 5
                + 0.6696768 * (min(egfr, 60) - 60) / (-15)
                + 0.0237374 * (max(egfr, 60) - 90) / (-15)
                + 0.2688352 * bptreat
                - 0.0434892 * bptreat * (max(sbp, 110) - 130) / 20
                - 0.1297155 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.299086 * (age - 55) / 10 * dm
                - 0.1079522 * (age - 55) / 10 * smoking
                + 0.0130483 * (age - 55) / 10 * (max(bmi, 30) - 30) / 5
                - 0.1797791 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (sdi_contrast(0.1235632, 0.3592212) if not pd.isna(sdi) else 0.17924)
            )

            logor_30yr_HF = (
                -2.060187
                + 0.5387527 * ((age - 55) / 10)
                - 0.1090333 * (((age - 55) / 10) ** 2)
                - 0.4829094 * (min(sbp, 110) - 110) / 20
                + 0.2843569 * (max(sbp, 110) - 130) / 20
                + 0.6827667 * dm
                + 0.2406677 * smoking
                + 0.0618028 * (min(bmi, 30) - 25) / 5
                + 0.2705615 * (max(bmi, 30) - 30) / 5
                + 0.2255837 * (min(egfr, 60) - 60) / (-15)
                + 0.0653632 * (max(egfr, 60) - 90) / (-15)
                + 0.2263243 * bptreat
                - 0.0316851 * bptreat * (max(sbp, 110) - 130) / 20
                - 0.1258716 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.3243709 * (age - 55) / 10 * dm
                - 0.1596172 * (age - 55) / 10 * smoking
                - 0.0103092 * (age - 55) / 10 * (max(bmi, 30) - 30) / 5
                - 0.1204785 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (sdi_contrast(0.0680528, 0.2619865) if not pd.isna(sdi) else 0.1151424)
            )

    return {
        "prevent_sdi_10yr_CVD": sigmoid_pct(logor_10yr_CVD) if not pd.isna(logor_10yr_CVD) else np.nan,
        "prevent_sdi_10yr_ASCVD": sigmoid_pct(logor_10yr_ASCVD) if not pd.isna(logor_10yr_ASCVD) else np.nan,
        "prevent_sdi_10yr_HF": sigmoid_pct(logor_10yr_HF) if not pd.isna(logor_10yr_HF) else np.nan,
        "prevent_sdi_30yr_CVD": sigmoid_pct(logor_30yr_CVD) if not pd.isna(logor_30yr_CVD) else np.nan,
        "prevent_sdi_30yr_ASCVD": sigmoid_pct(logor_30yr_ASCVD) if not pd.isna(logor_30yr_ASCVD) else np.nan,
        "prevent_sdi_30yr_HF": sigmoid_pct(logor_30yr_HF) if not pd.isna(logor_30yr_HF) else np.nan,
    }

