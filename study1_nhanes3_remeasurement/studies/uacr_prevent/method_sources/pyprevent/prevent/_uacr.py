from __future__ import annotations

import math

import numpy as np
import pandas as pd

from ._core import adjust, invalid_uacr, mmol_conversion, sigmoid_pct, validate_common_inputs


def _uacr_nan() -> dict[str, float]:
    return {
        "prevent_uacr_10yr_CVD": np.nan,
        "prevent_uacr_10yr_ASCVD": np.nan,
        "prevent_uacr_10yr_HF": np.nan,
        "prevent_uacr_30yr_CVD": np.nan,
        "prevent_uacr_30yr_ASCVD": np.nan,
        "prevent_uacr_30yr_HF": np.nan,
    }


def prevent_uacr(sex, age, tc, hdl, sbp, dm, smoking, bmi, egfr, bptreat, statin, uacr) -> dict[str, float]:
    """
    R parity: AHAprevent::pred_risk_uacr (10yr + 30yr).
    """
    if not validate_common_inputs(age, sex, sbp, dm, smoking, egfr):
        return _uacr_nan()
    if invalid_uacr(uacr):
        return _uacr_nan()

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

    logor_10yr_CVD = np.nan
    logor_10yr_ASCVD = np.nan
    logor_10yr_HF = np.nan
    logor_30yr_CVD = np.nan
    logor_30yr_ASCVD = np.nan
    logor_30yr_HF = np.nan

    if sex == 1.0:  # female
        if can_cvd_ascvd:
            logor_10yr_CVD = (
                -3.738341
                + 0.7969249 * ((age - 55) / 10)
                + 0.0256635 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.1588107 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.2255701 * (min(sbp, 110) - 110) / 20
                + 0.3396649 * (max(sbp, 110) - 130) / 20
                + 0.8047515 * dm
                + 0.5285338 * smoking
                + 0.4803511 * (min(egfr, 60) - 60) / (-15)
                + 0.0434472 * (max(egfr, 60) - 90) / (-15)
                + 0.2985207 * bptreat
                - 0.1497787 * statin
                - 0.0742889 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.106756 * statin * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0778126 * (age - 55) / 10 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                + 0.0306768 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.0907168 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.2705122 * (age - 55) / 10 * dm
                - 0.0830564 * (age - 55) / 10 * smoking
                - 0.1389249 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (0.1793037 * math.log(adjust(uacr)) if not pd.isna(uacr) else 0.0132073)
            )

            logor_10yr_ASCVD = (
                -4.174614
                + 0.7201999 * ((age - 55) / 10)
                + 0.1135771 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.1493506 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.0726677 * (min(sbp, 110) - 110) / 20
                + 0.3436259 * (max(sbp, 110) - 130) / 20
                + 0.7773094 * dm
                + 0.4746662 * smoking
                + 0.3824646 * (min(egfr, 60) - 60) / (-15)
                + 0.0394178 * (max(egfr, 60) - 90) / (-15)
                + 0.2125182 * bptreat
                - 0.0603046 * statin
                - 0.0466053 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.0733118 * statin * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0534262 * (age - 55) / 10 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                + 0.0325689 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.0999887 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.2411762 * (age - 55) / 10 * dm
                - 0.0826941 * (age - 55) / 10 * smoking
                - 0.1444737 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (0.1501217 * math.log(adjust(uacr)) if not pd.isna(uacr) else 0.0050257)
            )

            logor_30yr_CVD = (
                -1.583738
                + 0.5491768 * ((age - 55) / 10)
                - 0.0937311 * (((age - 55) / 10) ** 2)
                + 0.0359847 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.1642965 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.1483404 * (min(sbp, 110) - 110) / 20
                + 0.313353 * (max(sbp, 110) - 130) / 20
                + 0.6253766 * dm
                + 0.3147172 * smoking
                + 0.1094663 * (min(egfr, 60) - 60) / (-15)
                + 0.0550705 * (max(egfr, 60) - 90) / (-15)
                + 0.2782433 * bptreat
                - 0.0786239 * statin
                - 0.0628947 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.093204 * statin * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0710685 * (age - 55) / 10 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                + 0.0306363 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.0951455 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.3168231 * (age - 55) / 10 * dm
                - 0.1636391 * (age - 55) / 10 * smoking
                - 0.1265483 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (0.1142251 * math.log(adjust(uacr)) if not pd.isna(uacr) else -0.0055863)
            )

            logor_30yr_ASCVD = (
                -2.178888
                + 0.4629669 * ((age - 55) / 10)
                - 0.0902777 * (((age - 55) / 10) ** 2)
                + 0.1215214 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.1522069 * (mmol_conversion(hdl) - 1.3) / 0.3
                + 0.0092679 * (min(sbp, 110) - 110) / 20
                + 0.3113609 * (max(sbp, 110) - 130) / 20
                + 0.581256 * dm
                + 0.263167 * smoking
                + 0.0391726 * (min(egfr, 60) - 60) / (-15)
                + 0.0492959 * (max(egfr, 60) - 90) / (-15)
                + 0.1786178 * bptreat
                + 0.0131058 * statin
                - 0.0325135 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.0617093 * statin * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0489189 * (age - 55) / 10 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                + 0.0321079 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.1003185 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.2684574 * (age - 55) / 10 * dm
                - 0.1547301 * (age - 55) / 10 * smoking
                - 0.1130703 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (0.0903471 * math.log(adjust(uacr)) if not pd.isna(uacr) else -0.0145818)
            )

        if can_hf:
            logor_10yr_HF = (
                -4.841506
                + 0.9145975 * ((age - 55) / 10)
                - 0.4441346 * (min(sbp, 110) - 110) / 20
                + 0.3260323 * (max(sbp, 110) - 130) / 20
                + 0.9611365 * dm
                + 0.5755787 * smoking
                + 0.0008831 * (min(bmi, 30) - 25) / 5
                + 0.2988964 * (max(bmi, 30) - 30) / 5
                + 0.5915291 * (min(egfr, 60) - 60) / (-15)
                + 0.0556823 * (max(egfr, 60) - 90) / (-15)
                + 0.3314097 * bptreat
                - 0.1078596 * bptreat * (max(sbp, 110) - 130) / 20
                - 0.0875231 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.356859 * (age - 55) / 10 * dm
                - 0.1220248 * (age - 55) / 10 * smoking
                - 0.0053637 * (age - 55) / 10 * (max(bmi, 30) - 30) / 5
                - 0.1610389 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (0.2197281 * math.log(adjust(uacr)) if not pd.isna(uacr) else 0.0326667)
            )

            logor_30yr_HF = (
                -2.538952
                + 0.6319513 * ((age - 55) / 10)
                - 0.1009284 * (((age - 55) / 10) ** 2)
                - 0.3787175 * (min(sbp, 110) - 110) / 20
                + 0.2863393 * (max(sbp, 110) - 130) / 20
                + 0.7631221 * dm
                + 0.3355843 * smoking
                + 0.0677084 * (min(bmi, 30) - 25) / 5
                + 0.2517238 * (max(bmi, 30) - 30) / 5
                + 0.1940067 * (min(egfr, 60) - 60) / (-15)
                + 0.0664006 * (max(egfr, 60) - 90) / (-15)
                + 0.3171436 * bptreat
                - 0.0970661 * bptreat * (max(sbp, 110) - 130) / 20
                - 0.0896239 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.400743 * (age - 55) / 10 * dm
                - 0.2042041 * (age - 55) / 10 * smoking
                - 0.0054699 * (age - 55) / 10 * (max(bmi, 30) - 30) / 5
                - 0.13602 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (0.1486028 * math.log(adjust(uacr)) if not pd.isna(uacr) else 0.011608)
            )
    else:  # male
        if can_cvd_ascvd:
            logor_10yr_CVD = (
                -3.510705
                + 0.7768655 * ((age - 55) / 10)
                + 0.0659949 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0951111 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.420667 * (min(sbp, 110) - 110) / 20
                + 0.3120151 * (max(sbp, 110) - 130) / 20
                + 0.698521 * dm
                + 0.4314669 * smoking
                + 0.3841364 * (min(egfr, 60) - 60) / (-15)
                + 0.009384 * (max(egfr, 60) - 90) / (-15)
                + 0.2676494 * bptreat
                - 0.1390966 * statin
                - 0.0579315 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.1383719 * statin * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0488332 * (age - 55) / 10 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                + 0.0200406 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.102454 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.2236355 * (age - 55) / 10 * dm
                - 0.089485 * (age - 55) / 10 * smoking
                - 0.1321848 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (0.1887974 * math.log(adjust(uacr)) if not pd.isna(uacr) else 0.0916979)
            )

            logor_10yr_ASCVD = (
                -3.85146
                + 0.7141718 * ((age - 55) / 10)
                + 0.1602194 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.1139086 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.2719456 * (min(sbp, 110) - 110) / 20
                + 0.3058719 * (max(sbp, 110) - 130) / 20
                + 0.6600631 * dm
                + 0.3884022 * smoking
                + 0.2466316 * (min(egfr, 60) - 60) / (-15)
                + 0.0151852 * (max(egfr, 60) - 90) / (-15)
                + 0.186167 * bptreat
                - 0.0894395 * statin
                - 0.0411884 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.1058212 * statin * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.028089 * (age - 55) / 10 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                + 0.0240427 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.0912325 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.2004894 * (age - 55) / 10 * dm
                - 0.096936 * (age - 55) / 10 * smoking
                - 0.1022867 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (0.1510073 * math.log(adjust(uacr)) if not pd.isna(uacr) else 0.0556)
            )

            logor_30yr_CVD = (
                -1.398727
                + 0.464491 * ((age - 55) / 10)
                - 0.0998895 * (((age - 55) / 10) ** 2)
                + 0.0757606 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.1031778 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.1990714 * (min(sbp, 110) - 110) / 20
                + 0.2715816 * (max(sbp, 110) - 130) / 20
                + 0.4754637 * dm
                + 0.2069672 * smoking
                + 0.0331103 * (min(egfr, 60) - 60) / (-15)
                + 0.0540474 * (max(egfr, 60) - 90) / (-15)
                + 0.2189911 * bptreat
                - 0.0331044 * statin
                - 0.04534 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.1214535 * statin * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0483995 * (age - 55) / 10 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                + 0.0178997 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.1059324 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.2492861 * (age - 55) / 10 * dm
                - 0.1561543 * (age - 55) / 10 * smoking
                - 0.1012429 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (0.1007571 * math.log(adjust(uacr)) if not pd.isna(uacr) else 0.0572456)
            )

            logor_30yr_ASCVD = (
                -1.873449
                + 0.3995607 * ((age - 55) / 10)
                - 0.094557 * (((age - 55) / 10) ** 2)
                + 0.1686692 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.1202145 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.0555561 * (min(sbp, 110) - 110) / 20
                + 0.2633566 * (max(sbp, 110) - 130) / 20
                + 0.4362036 * dm
                + 0.1716233 * smoking
                - 0.0775282 * (min(egfr, 60) - 60) / (-15)
                + 0.0561236 * (max(egfr, 60) - 90) / (-15)
                + 0.1319331 * bptreat
                + 0.0102428 * statin
                - 0.0269294 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.0920557 * statin * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0297021 * (age - 55) / 10 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                + 0.0217935 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.0893347 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.2081467 * (age - 55) / 10 * dm
                - 0.1542716 * (age - 55) / 10 * smoking
                - 0.0597254 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (0.0684872 * math.log(adjust(uacr)) if not pd.isna(uacr) else 0.0193962)
            )

        if can_hf:
            logor_10yr_HF = (
                -4.556907
                + 0.9111795 * ((age - 55) / 10)
                - 0.6693649 * (min(sbp, 110) - 110) / 20
                + 0.3290082 * (max(sbp, 110) - 130) / 20
                + 0.8377655 * dm
                + 0.4978917 * smoking
                - 0.042749 * (min(bmi, 30) - 25) / 5
                + 0.3624165 * (max(bmi, 30) - 30) / 5
                + 0.5075796 * (min(egfr, 60) - 60) / (-15)
                + 0.0137716 * (max(egfr, 60) - 90) / (-15)
                + 0.2739963 * bptreat
                - 0.0645712 * bptreat * (max(sbp, 110) - 130) / 20
                - 0.1230039 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.3013297 * (age - 55) / 10 * dm
                - 0.1410318 * (age - 55) / 10 * smoking
                + 0.0021531 * (age - 55) / 10 * (max(bmi, 30) - 30) / 5
                - 0.1548018 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (0.2306299 * math.log(adjust(uacr)) if not pd.isna(uacr) else 0.1472194)
            )

            logor_30yr_HF = (
                -2.314872
                + 0.5750236 * ((age - 55) / 10)
                - 0.1062268 * (((age - 55) / 10) ** 2)
                - 0.4633994 * (min(sbp, 110) - 110) / 20
                + 0.2742874 * (max(sbp, 110) - 130) / 20
                + 0.612208 * dm
                + 0.2614987 * smoking
                + 0.0895459 * (min(bmi, 30) - 25) / 5
                + 0.2632424 * (max(bmi, 30) - 30) / 5
                + 0.1430472 * (min(egfr, 60) - 60) / (-15)
                + 0.0535184 * (max(egfr, 60) - 90) / (-15)
                + 0.2417468 * bptreat
                - 0.0498574 * bptreat * (max(sbp, 110) - 130) / 20
                - 0.1193827 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.316651 * (age - 55) / 10 * dm
                - 0.2046122 * (age - 55) / 10 * smoking
                - 0.0216878 * (age - 55) / 10 * (max(bmi, 30) - 30) / 5
                - 0.1165637 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (0.1366452 * math.log(adjust(uacr)) if not pd.isna(uacr) else 0.1078355)
            )

    return {
        "prevent_uacr_10yr_CVD": sigmoid_pct(logor_10yr_CVD) if not pd.isna(logor_10yr_CVD) else np.nan,
        "prevent_uacr_10yr_ASCVD": sigmoid_pct(logor_10yr_ASCVD) if not pd.isna(logor_10yr_ASCVD) else np.nan,
        "prevent_uacr_10yr_HF": sigmoid_pct(logor_10yr_HF) if not pd.isna(logor_10yr_HF) else np.nan,
        "prevent_uacr_30yr_CVD": sigmoid_pct(logor_30yr_CVD) if not pd.isna(logor_30yr_CVD) else np.nan,
        "prevent_uacr_30yr_ASCVD": sigmoid_pct(logor_30yr_ASCVD) if not pd.isna(logor_30yr_ASCVD) else np.nan,
        "prevent_uacr_30yr_HF": sigmoid_pct(logor_30yr_HF) if not pd.isna(logor_30yr_HF) else np.nan,
    }

