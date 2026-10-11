from __future__ import annotations

import numpy as np
import pandas as pd

from ._core import invalid_hba1c, mmol_conversion, sigmoid_pct, validate_common_inputs


def _hba1c_nan() -> dict[str, float]:
    return {
        "prevent_hba1c_10yr_CVD": np.nan,
        "prevent_hba1c_10yr_ASCVD": np.nan,
        "prevent_hba1c_10yr_HF": np.nan,
        "prevent_hba1c_30yr_CVD": np.nan,
        "prevent_hba1c_30yr_ASCVD": np.nan,
        "prevent_hba1c_30yr_HF": np.nan,
    }


def prevent_hba1c(sex, age, tc, hdl, sbp, dm, smoking, bmi, egfr, bptreat, statin, hba1c) -> dict[str, float]:
    """
    R parity: AHAprevent::pred_risk_hba1c (10yr + 30yr).
    """
    if not validate_common_inputs(age, sex, sbp, dm, smoking, egfr):
        return _hba1c_nan()
    if invalid_hba1c(hba1c):
        return _hba1c_nan()

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
                -3.306162
                + 0.7858178 * ((age - 55) / 10)
                + 0.0194438 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.1521964 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.2296681 * (min(sbp, 110) - 110) / 20
                + 0.3465777 * (max(sbp, 110) - 130) / 20
                + 0.5366241 * dm
                + 0.5411682 * smoking
                + 0.5931898 * (min(egfr, 60) - 60) / (-15)
                + 0.0472458 * (max(egfr, 60) - 90) / (-15)
                + 0.3158567 * bptreat
                - 0.1535174 * statin
                - 0.0687752 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.1054746 * statin * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0761119 * (age - 55) / 10 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                + 0.0307469 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.0905966 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.2241857 * (age - 55) / 10 * dm
                - 0.080186 * (age - 55) / 10 * smoking
                - 0.1667286 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (
                    0.1338348 * (hba1c - 5.3) * dm + 0.1622409 * (hba1c - 5.3) * (1 - dm)
                    if not pd.isna(hba1c)
                    else -0.0142496
                )
            )

            logor_10yr_ASCVD = (
                -3.838746
                + 0.7111831 * ((age - 55) / 10)
                + 0.106797 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.1425745 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.0736824 * (min(sbp, 110) - 110) / 20
                + 0.3480844 * (max(sbp, 110) - 130) / 20
                + 0.5112951 * dm
                + 0.4880292 * smoking
                + 0.4754997 * (min(egfr, 60) - 60) / (-15)
                + 0.0438132 * (max(egfr, 60) - 90) / (-15)
                + 0.2259093 * bptreat
                - 0.0648872 * statin
                - 0.0437645 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.0697082 * statin * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0506382 * (age - 55) / 10 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                + 0.0327475 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.0996442 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.1924338 * (age - 55) / 10 * dm
                - 0.0803539 * (age - 55) / 10 * smoking
                - 0.1682586 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (
                    0.1339055 * (hba1c - 5.3) * dm + 0.1596461 * (hba1c - 5.3) * (1 - dm)
                    if not pd.isna(hba1c)
                    else 0.0015678
                )
            )

            logor_30yr_CVD = (
                -1.341059
                + 0.5343493 * ((age - 55) / 10)
                - 0.0952314 * (((age - 55) / 10) ** 2)
                + 0.0298124 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.1578451 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.1504488 * (min(sbp, 110) - 110) / 20
                + 0.3173368 * (max(sbp, 110) - 130) / 20
                + 0.4314738 * dm
                + 0.3209399 * smoking
                + 0.1771435 * (min(egfr, 60) - 60) / (-15)
                + 0.0582828 * (max(egfr, 60) - 90) / (-15)
                + 0.2888947 * bptreat
                - 0.0795886 * statin
                - 0.0600438 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.0920598 * statin * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0696108 * (age - 55) / 10 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                + 0.0308807 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.0954051 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.2763408 * (age - 55) / 10 * dm
                - 0.1623944 * (age - 55) / 10 * smoking
                - 0.1430514 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (
                    0.0940543 * (hba1c - 5.3) * dm + 0.1116486 * (hba1c - 5.3) * (1 - dm)
                    if not pd.isna(hba1c)
                    else -0.0024798
                )
            )

            logor_30yr_ASCVD = (
                -2.011533
                + 0.4555574 * ((age - 55) / 10)
                - 0.0903501 * (((age - 55) / 10) ** 2)
                + 0.1148321 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.1458754 * (mmol_conversion(hdl) - 1.3) / 0.3
                + 0.0089323 * (min(sbp, 110) - 110) / 20
                + 0.3139029 * (max(sbp, 110) - 130) / 20
                + 0.386281 * dm
                + 0.2714309 * smoking
                + 0.0930987 * (min(egfr, 60) - 60) / (-15)
                + 0.0532216 * (max(egfr, 60) - 90) / (-15)
                + 0.1862181 * bptreat
                + 0.0106964 * statin
                - 0.0329713 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.0583609 * statin * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0463273 * (age - 55) / 10 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                + 0.0324717 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.1004777 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.2266944 * (age - 55) / 10 * dm
                - 0.1541859 * (age - 55) / 10 * smoking
                - 0.1286005 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (
                    0.0875827 * (hba1c - 5.3) * dm + 0.1126417 * (hba1c - 5.3) * (1 - dm)
                    if not pd.isna(hba1c)
                    else 0.0124356
                )
            )

        if can_hf:
            logor_10yr_HF = (
                -4.288225
                + 0.8997391 * ((age - 55) / 10)
                - 0.4422749 * (min(sbp, 110) - 110) / 20
                + 0.3378691 * (max(sbp, 110) - 130) / 20
                + 0.681284 * dm
                + 0.5886005 * smoking
                - 0.0148657 * (min(bmi, 30) - 25) / 5
                + 0.2958374 * (max(bmi, 30) - 30) / 5
                + 0.73447 * (min(egfr, 60) - 60) / (-15)
                + 0.05926 * (max(egfr, 60) - 90) / (-15)
                + 0.3543475 * bptreat
                - 0.1002139 * bptreat * (max(sbp, 110) - 130) / 20
                - 0.0878765 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.303684 * (age - 55) / 10 * dm
                - 0.1178943 * (age - 55) / 10 * smoking
                - 0.008345 * (age - 55) / 10 * (max(bmi, 30) - 30) / 5
                - 0.1912183 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (
                    0.1856442 * (hba1c - 5.3) * dm + 0.1833083 * (hba1c - 5.3) * (1 - dm)
                    if not pd.isna(hba1c)
                    else -0.0143112
                )
            )

            logor_30yr_HF = (
                -2.193553
                + 0.6210856 * ((age - 55) / 10)
                - 0.1000972 * (((age - 55) / 10) ** 2)
                - 0.3773697 * (min(sbp, 110) - 110) / 20
                + 0.295316 * (max(sbp, 110) - 130) / 20
                + 0.5681692 * dm
                + 0.3449139 * smoking
                + 0.0540094 * (min(bmi, 30) - 25) / 5
                + 0.249767 * (max(bmi, 30) - 30) / 5
                + 0.2875781 * (min(egfr, 60) - 60) / (-15)
                + 0.0692013 * (max(egfr, 60) - 90) / (-15)
                + 0.3334936 * bptreat
                - 0.0922339 * bptreat * (max(sbp, 110) - 130) / 20
                - 0.0907885 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.3554646 * (age - 55) / 10 * dm
                - 0.2008846 * (age - 55) / 10 * smoking
                - 0.0079611 * (age - 55) / 10 * (max(bmi, 30) - 30) / 5
                - 0.156803 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (
                    0.1448336 * (hba1c - 5.3) * dm + 0.1277838 * (hba1c - 5.3) * (1 - dm)
                    if not pd.isna(hba1c)
                    else -0.0022589
                )
            )
    else:  # male
        if can_cvd_ascvd:
            logor_10yr_CVD = (
                -3.040901
                + 0.7699177 * ((age - 55) / 10)
                + 0.0605093 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0888525 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.417713 * (min(sbp, 110) - 110) / 20
                + 0.3288657 * (max(sbp, 110) - 130) / 20
                + 0.4759471 * dm
                + 0.4385663 * smoking
                + 0.5334616 * (min(egfr, 60) - 60) / (-15)
                + 0.0206431 * (max(egfr, 60) - 90) / (-15)
                + 0.2917524 * bptreat
                - 0.1383313 * statin
                - 0.0482622 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.1393796 * statin * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0463501 * (age - 55) / 10 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                + 0.0205926 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.1037717 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.1737697 * (age - 55) / 10 * dm
                - 0.0915839 * (age - 55) / 10 * smoking
                - 0.1637039 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (
                    0.13159 * (hba1c - 5.3) * dm + 0.1295185 * (hba1c - 5.3) * (1 - dm)
                    if not pd.isna(hba1c)
                    else -0.0128373
                )
            )

            logor_10yr_ASCVD = (
                -3.51835
                + 0.7064146 * ((age - 55) / 10)
                + 0.1532267 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.1082166 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.2675288 * (min(sbp, 110) - 110) / 20
                + 0.3173809 * (max(sbp, 110) - 130) / 20
                + 0.432604 * dm
                + 0.3958842 * smoking
                + 0.3665014 * (min(egfr, 60) - 60) / (-15)
                + 0.0250243 * (max(egfr, 60) - 90) / (-15)
                + 0.2061158 * bptreat
                - 0.0899988 * statin
                - 0.0334959 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.1034168 * statin * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0255406 * (age - 55) / 10 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                + 0.0247538 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.0917441 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.1499195 * (age - 55) / 10 * dm
                - 0.098089 * (age - 55) / 10 * smoking
                - 0.1305231 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (
                    0.1157161 * (hba1c - 5.3) * dm + 0.1288303 * (hba1c - 5.3) * (1 - dm)
                    if not pd.isna(hba1c)
                    else -0.0010001
                )
            )

            logor_30yr_CVD = (
                -1.180767
                + 0.4519873 * ((age - 55) / 10)
                - 0.101624 * (((age - 55) / 10) ** 2)
                + 0.0700456 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0968005 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.1923527 * (min(sbp, 110) - 110) / 20
                + 0.2827043 * (max(sbp, 110) - 130) / 20
                + 0.3417152 * dm
                + 0.2105272 * smoking
                + 0.1113291 * (min(egfr, 60) - 60) / (-15)
                + 0.0640135 * (max(egfr, 60) - 90) / (-15)
                + 0.2334248 * bptreat
                - 0.0299421 * statin
                - 0.0393204 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.1228854 * statin * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0463737 * (age - 55) / 10 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                + 0.0184599 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.1085744 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.2208049 * (age - 55) / 10 * dm
                - 0.1577978 * (age - 55) / 10 * smoking
                - 0.1179375 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (
                    0.0768169 * (hba1c - 5.3) * dm + 0.0777295 * (hba1c - 5.3) * (1 - dm)
                    if not pd.isna(hba1c)
                    else 0.0092204
                )
            )

            logor_30yr_ASCVD = (
                -1.777708
                + 0.3883267 * ((age - 55) / 10)
                - 0.0958114 * (((age - 55) / 10) ** 2)
                + 0.1613374 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.1144418 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.0474338 * (min(sbp, 110) - 110) / 20
                + 0.2691281 * (max(sbp, 110) - 130) / 20
                + 0.2859773 * dm
                + 0.1759553 * smoking
                - 0.0242898 * (min(egfr, 60) - 60) / (-15)
                + 0.0644523 * (max(egfr, 60) - 90) / (-15)
                + 0.142874 * bptreat
                + 0.0115062 * statin
                - 0.02333 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.0899664 * statin * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0275478 * (age - 55) / 10 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                + 0.022573 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.090802 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.1771894 * (age - 55) / 10 * dm
                - 0.1548847 * (age - 55) / 10 * smoking
                - 0.0732754 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (
                    0.0591089 * (hba1c - 5.3) * dm + 0.0821158 * (hba1c - 5.3) * (1 - dm)
                    if not pd.isna(hba1c)
                    else 0.0179755
                )
            )

        if can_hf:
            logor_10yr_HF = (
                -3.961954
                + 0.911787 * ((age - 55) / 10)
                - 0.6568071 * (min(sbp, 110) - 110) / 20
                + 0.3524645 * (max(sbp, 110) - 130) / 20
                + 0.5849752 * dm
                + 0.5014014 * smoking
                - 0.0512352 * (min(bmi, 30) - 25) / 5
                + 0.365294 * (max(bmi, 30) - 30) / 5
                + 0.6892219 * (min(egfr, 60) - 60) / (-15)
                + 0.0292377 * (max(egfr, 60) - 90) / (-15)
                + 0.3038296 * bptreat
                - 0.0515032 * bptreat * (max(sbp, 110) - 130) / 20
                - 0.1262343 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.2449514 * (age - 55) / 10 * dm
                - 0.1392217 * (age - 55) / 10 * smoking
                + 0.0009592 * (age - 55) / 10 * (max(bmi, 30) - 30) / 5
                - 0.1917105 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (
                    0.1652857 * (hba1c - 5.3) * dm + 0.1505859 * (hba1c - 5.3) * (1 - dm)
                    if not pd.isna(hba1c)
                    else -0.0113444
                )
            )

            logor_30yr_HF = (
                -1.974999
                + 0.5703729 * ((age - 55) / 10)
                - 0.1084544 * (((age - 55) / 10) ** 2)
                - 0.4471767 * (min(sbp, 110) - 110) / 20
                + 0.2910152 * (max(sbp, 110) - 130) / 20
                + 0.4507242 * dm
                + 0.259585 * smoking
                + 0.0850676 * (min(bmi, 30) - 25) / 5
                + 0.2637222 * (max(bmi, 30) - 30) / 5
                + 0.2454706 * (min(egfr, 60) - 60) / (-15)
                + 0.0675649 * (max(egfr, 60) - 90) / (-15)
                + 0.2611991 * bptreat
                - 0.0408908 * bptreat * (max(sbp, 110) - 130) / 20
                - 0.1241051 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.2849461 * (age - 55) / 10 * dm
                - 0.2032308 * (age - 55) / 10 * smoking
                - 0.0239714 * (age - 55) / 10 * (max(bmi, 30) - 30) / 5
                - 0.138301 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (
                    0.1101184 * (hba1c - 5.3) * dm + 0.0949198 * (hba1c - 5.3) * (1 - dm)
                    if not pd.isna(hba1c)
                    else 0.0084192
                )
            )

    return {
        "prevent_hba1c_10yr_CVD": sigmoid_pct(logor_10yr_CVD) if not pd.isna(logor_10yr_CVD) else np.nan,
        "prevent_hba1c_10yr_ASCVD": sigmoid_pct(logor_10yr_ASCVD) if not pd.isna(logor_10yr_ASCVD) else np.nan,
        "prevent_hba1c_10yr_HF": sigmoid_pct(logor_10yr_HF) if not pd.isna(logor_10yr_HF) else np.nan,
        "prevent_hba1c_30yr_CVD": sigmoid_pct(logor_30yr_CVD) if not pd.isna(logor_30yr_CVD) else np.nan,
        "prevent_hba1c_30yr_ASCVD": sigmoid_pct(logor_30yr_ASCVD) if not pd.isna(logor_30yr_ASCVD) else np.nan,
        "prevent_hba1c_30yr_HF": sigmoid_pct(logor_30yr_HF) if not pd.isna(logor_30yr_HF) else np.nan,
    }

