from __future__ import annotations

import numpy as np
import pandas as pd

from ._core import mmol_conversion, sigmoid_pct, validate_common_inputs


def prevent_base(sex, age, tc, hdl, sbp, dm, smoking, bmi, egfr, bptreat, statin) -> dict[str, float]:
    """
    R parity: AHAprevent::pred_risk_base (10yr + 30yr).
    """
    if not validate_common_inputs(age, sex, sbp, dm, smoking, egfr):
        return {
            "prevent_base_10yr_CVD": np.nan,
            "prevent_base_10yr_ASCVD": np.nan,
            "prevent_base_10yr_HF": np.nan,
            "prevent_base_30yr_CVD": np.nan,
            "prevent_base_30yr_ASCVD": np.nan,
            "prevent_base_30yr_HF": np.nan,
        }

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
                -3.307728
                + 0.7939329 * (age - 55) / 10
                + 0.0305239 * (mmol_conversion(tc - hdl) - 3.5)
                - 0.1606857 * (mmol_conversion(hdl) - 1.3) / (0.3)
                - 0.2394003 * (min(sbp, 110) - 110) / 20
                + 0.360078 * (max(sbp, 110) - 130) / 20
                + 0.8667604 * (dm)
                + 0.5360739 * (smoking)
                + 0.6045917 * (min(egfr, 60) - 60) / (-15)
                + 0.0433769 * (max(egfr, 60) - 90) / (-15)
                + 0.3151672 * (bptreat)
                - 0.1477655 * (statin)
                - 0.0663612 * (bptreat) * (max(sbp, 110) - 130) / 20
                + 0.1197879 * (statin) * (mmol_conversion(tc - hdl) - 3.5)
                - 0.0819715 * (age - 55) / 10 * (mmol_conversion(tc - hdl) - 3.5)
                + 0.0306769 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / (0.3)
                - 0.0946348 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.27057 * (age - 55) / 10 * (dm)
                - 0.078715 * (age - 55) / 10 * (smoking)
                - 0.1637806 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
            )

            logor_10yr_ASCVD = (
                -3.819975
                + 0.719883 * (age - 55) / 10
                + 0.1176967 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.151185 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.0835358 * (min(sbp, 110) - 110) / 20
                + 0.3592852 * (max(sbp, 110) - 130) / 20
                + 0.8348585 * (dm)
                + 0.4831078 * (smoking)
                + 0.4864619 * (min(egfr, 60) - 60) / (-15)
                + 0.0397779 * (max(egfr, 60) - 90) / (-15)
                + 0.2265309 * (bptreat)
                - 0.0592374 * (statin)
                - 0.0395762 * (bptreat) * (max(sbp, 110) - 130) / 20
                + 0.0844423 * (statin) * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0567839 * (age - 55) / 10 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                + 0.0325692 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.1035985 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.2417542 * (age - 55) / 10 * (dm)
                - 0.0791142 * (age - 55) / 10 * (smoking)
                - 0.1671492 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
            )

            logor_30yr_CVD = (
                -1.318827
                + 0.5503079 * (age - 55) / 10
                - 0.0928369 * (((age - 55) / 10) ** 2)
                + 0.0409794 * (mmol_conversion(tc - hdl) - 3.5)
                + (-0.1663306) * (mmol_conversion(hdl) - 1.3) / 0.3
                + (-0.1628654) * (min(sbp, 110) - 110) / 20
                + 0.3299505 * (max(sbp, 110) - 130) / 20
                + 0.6793894 * dm
                + 0.3196112 * smoking
                + 0.1857101 * (min(egfr, 60) - 60) / (-15)
                + 0.0553528 * (max(egfr, 60) - 90) / (-15)
                + 0.2894 * bptreat
                + (-0.075688) * statin
                + (-0.056367) * bptreat * (max(sbp, 110) - 130) / 20
                + 0.1071019 * statin * (mmol_conversion(tc - hdl) - 3.5)
                + (-0.0751438) * (age - 55) / 10 * (mmol_conversion(tc - hdl) - 3.5)
                + 0.0301786 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                + (-0.0998776) * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                + (-0.3206166) * (age - 55) / 10 * dm
                + (-0.1607862) * (age - 55) / 10 * smoking
                + (-0.1450788) * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
            )

            logor_30yr_ASCVD = (
                -1.974074
                + 0.4669202 * ((age - 55) / 10)
                - 0.0893118 * (((age - 55) / 10) ** 2)
                + 0.1256901 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.1542255 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.0018093 * (min(sbp, 110) - 110) / 20
                + 0.322949 * (max(sbp, 110) - 130) / 20
                + 0.6296707 * dm
                + 0.268292 * smoking
                + 0.100106 * (min(egfr, 60) - 60) / (-15)
                + 0.0499663 * (max(egfr, 60) - 90) / (-15)
                + 0.1875292 * bptreat
                + 0.0152476 * statin
                - 0.0276123 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.0736147 * statin * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0521962 * (age - 55) / 10 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                + 0.0316918 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.1046101 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.2727793 * (age - 55) / 10 * dm
                - 0.1530907 * (age - 55) / 10 * smoking
                - 0.1299149 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
            )

        if can_hf:
            logor_10yr_HF = (
                -4.310409
                + 0.8998235 * (age - 55) / 10
                - 0.4559771 * (min(sbp, 110) - 110) / 20
                + 0.3576505 * (max(sbp, 110) - 130) / 20
                + 1.038346 * (dm)
                + 0.583916 * (smoking)
                - 0.0072294 * (min(bmi, 30) - 25) / 5
                + 0.2997706 * (max(bmi, 30) - 30) / 5
                + 0.7451638 * (min(egfr, 60) - 60) / (-15)
                + 0.0557087 * (max(egfr, 60) - 90) / (-15)
                + 0.3534442 * (bptreat)
                - 0.0981511 * (bptreat) * (max(sbp, 110) - 130) / 20
                - 0.0946663 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.3581041 * (age - 55) / 10 * (dm)
                - 0.1159453 * (age - 55) / 10 * (smoking)
                - 0.003878 * (age - 55) / 10 * (max(bmi, 30) - 30) / 5
                - 0.1884289 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
            )

            logor_30yr_HF = (
                -2.205379
                + 0.6254374 * ((age - 55) / 10)
                - 0.0983038 * (((age - 55) / 10) ** 2)
                - 0.3919241 * (min(sbp, 110) - 110) / 20
                + 0.3142295 * (max(sbp, 110) - 130) / 20
                + 0.8330787 * dm
                + 0.3438651 * smoking
                + 0.0594874 * (min(bmi, 30) - 25) / 5
                + 0.2525536 * (max(bmi, 30) - 30) / 5
                + 0.2981642 * (min(egfr, 60) - 60) / (-15)
                + 0.0667159 * (max(egfr, 60) - 90) / (-15)
                + 0.333921 * bptreat
                - 0.0893177 * bptreat * (max(sbp, 110) - 130) / 20
                - 0.0974299 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.404855 * (age - 55) / 10 * dm
                - 0.1982991 * (age - 55) / 10 * smoking
                - 0.0035619 * (age - 55) / 10 * (max(bmi, 30) - 30) / 5
                - 0.1564215 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
            )
    else:  # male
        if can_cvd_ascvd:
            logor_10yr_CVD = (
                -3.031168
                + 0.7688528 * (age - 55) / 10
                + 0.0736174 * (mmol_conversion(tc - hdl) - 3.5)
                - 0.0954431 * (mmol_conversion(hdl) - 1.3) / (0.3)
                - 0.4347345 * (min(sbp, 110) - 110) / 20
                + 0.3362658 * (max(sbp, 110) - 130) / 20
                + 0.7692857 * (dm)
                + 0.4386871 * (smoking)
                + 0.5378979 * (min(egfr, 60) - 60) / (-15)
                + 0.0164827 * (max(egfr, 60) - 90) / (-15)
                + 0.288879 * (bptreat)
                - 0.1337349 * (statin)
                - 0.0475924 * (bptreat) * (max(sbp, 110) - 130) / 20
                + 0.150273 * (statin) * (mmol_conversion(tc - hdl) - 3.5)
                - 0.0517874 * (age - 55) / 10 * (mmol_conversion(tc - hdl) - 3.5)
                + 0.0191169 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / (0.3)
                - 0.1049477 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.2251948 * (age - 55) / 10 * (dm)
                - 0.0895067 * (age - 55) / 10 * (smoking)
                - 0.1543702 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
            )

            logor_10yr_ASCVD = (
                -3.500655
                + 0.7099847 * (age - 55) / 10
                + 0.1658663 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.1144285 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.2837212 * (min(sbp, 110) - 110) / 20
                + 0.3239977 * (max(sbp, 110) - 130) / 20
                + 0.7189597 * (dm)
                + 0.3956973 * (smoking)
                + 0.3690075 * (min(egfr, 60) - 60) / (-15)
                + 0.0203619 * (max(egfr, 60) - 90) / (-15)
                + 0.2036522 * (bptreat)
                - 0.0865581 * (statin)
                - 0.0322916 * (bptreat) * (max(sbp, 110) - 130) / 20
                + 0.114563 * (statin) * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0300005 * (age - 55) / 10 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                + 0.0232747 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.0927024 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.2018525 * (age - 55) / 10 * (dm)
                - 0.0970527 * (age - 55) / 10 * (smoking)
                - 0.1217081 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
            )

            logor_30yr_CVD = (
                -1.148204
                + 0.4627309 * (age - 55) / 10
                - 0.0984281 * (((age - 55) / 10) ** 2)
                + 0.0836088 * (mmol_conversion(tc - hdl) - 3.5)
                + (-0.1029824) * (mmol_conversion(hdl) - 1.3) / 0.3
                + (-0.2140352) * (min(sbp, 110) - 110) / 20
                + 0.2904325 * (max(sbp, 110) - 130) / 20
                + 0.5331276 * dm
                + 0.2141914 * smoking
                + 0.1155556 * (min(egfr, 60) - 60) / (-15)
                + 0.0603775 * (max(egfr, 60) - 90) / (-15)
                + 0.232714 * bptreat
                + (-0.0272112) * statin
                + (-0.0384488) * bptreat * (max(sbp, 110) - 130) / 20
                + 0.134192 * statin * (mmol_conversion(tc - hdl) - 3.5)
                - 0.0511759 * (age - 55) / 10 * (mmol_conversion(tc - hdl) - 3.5)
                + 0.0165865 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.1101437 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.2585943 * (age - 55) / 10 * dm
                - 0.1566406 * (age - 55) / 10 * smoking
                - 0.1166776 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
            )

            logor_30yr_ASCVD = (
                -1.736444
                + 0.3994099 * ((age - 55) / 10)
                - 0.0937484 * (((age - 55) / 10) ** 2)
                + 0.1744643 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.120203 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.0665117 * (min(sbp, 110) - 110) / 20
                + 0.2753037 * (max(sbp, 110) - 130) / 20
                + 0.4790257 * dm
                + 0.1782635 * smoking
                - 0.0218789 * (min(egfr, 60) - 60) / (-15)
                + 0.0602553 * (max(egfr, 60) - 90) / (-15)
                + 0.1421182 * bptreat
                + 0.0135996 * statin
                - 0.0218265 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.1013148 * statin * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0312619 * (age - 55) / 10 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                + 0.020673 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.0920935 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.2159947 * (age - 55) / 10 * dm
                - 0.1548811 * (age - 55) / 10 * smoking
                - 0.0712547 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
            )

        if can_hf:
            logor_10yr_HF = (
                -3.946391
                + 0.8972642 * (age - 55) / 10
                - 0.6811466 * (min(sbp, 110) - 110) / 20
                + 0.3634461 * (max(sbp, 110) - 130) / 20
                + 0.923776 * (dm)
                + 0.5023736 * (smoking)
                - 0.0485841 * (min(bmi, 30) - 25) / 5
                + 0.3726929 * (max(bmi, 30) - 30) / 5
                + 0.6926917 * (min(egfr, 60) - 60) / (-15)
                + 0.0251827 * (max(egfr, 60) - 90) / (-15)
                + 0.2980922 * (bptreat)
                - 0.0497731 * (bptreat) * (max(sbp, 110) - 130) / 20
                - 0.1289201 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.3040924 * (age - 55) / 10 * (dm)
                - 0.1401688 * (age - 55) / 10 * (smoking)
                + 0.0068126 * (age - 55) / 10 * (max(bmi, 30) - 30) / 5
                - 0.1797778 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
            )

            logor_30yr_HF = (
                -1.95751
                + 0.5681541 * ((age - 55) / 10)
                - 0.1048388 * (((age - 55) / 10) ** 2)
                - 0.4761564 * (min(sbp, 110) - 110) / 20
                + 0.30324 * (max(sbp, 110) - 130) / 20
                + 0.6840338 * dm
                + 0.2656273 * smoking
                + 0.0833107 * (min(bmi, 30) - 25) / 5
                + 0.26999 * (max(bmi, 30) - 30) / 5
                + 0.2541805 * (min(egfr, 60) - 60) / (-15)
                + 0.0638923 * (max(egfr, 60) - 90) / (-15)
                + 0.2583631 * bptreat
                - 0.0391938 * bptreat * (max(sbp, 110) - 130) / 20
                - 0.1269124 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.3273572 * (age - 55) / 10 * dm
                - 0.2043019 * (age - 55) / 10 * smoking
                - 0.0182831 * (age - 55) / 10 * (max(bmi, 30) - 30) / 5
                - 0.1342618 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
            )

    return {
        "prevent_base_10yr_CVD": sigmoid_pct(logor_10yr_CVD) if not pd.isna(logor_10yr_CVD) else np.nan,
        "prevent_base_10yr_ASCVD": sigmoid_pct(logor_10yr_ASCVD) if not pd.isna(logor_10yr_ASCVD) else np.nan,
        "prevent_base_10yr_HF": sigmoid_pct(logor_10yr_HF) if not pd.isna(logor_10yr_HF) else np.nan,
        "prevent_base_30yr_CVD": sigmoid_pct(logor_30yr_CVD) if not pd.isna(logor_30yr_CVD) else np.nan,
        "prevent_base_30yr_ASCVD": sigmoid_pct(logor_30yr_ASCVD) if not pd.isna(logor_30yr_ASCVD) else np.nan,
        "prevent_base_30yr_HF": sigmoid_pct(logor_30yr_HF) if not pd.isna(logor_30yr_HF) else np.nan,
    }

