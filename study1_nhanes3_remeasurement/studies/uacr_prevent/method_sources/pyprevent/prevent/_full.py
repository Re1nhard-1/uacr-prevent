from __future__ import annotations

import math

import numpy as np
import pandas as pd

from ._core import (
    adjust,
    invalid_hba1c,
    invalid_sdi_decile,
    invalid_uacr,
    mmol_conversion,
    sdicat,
    sigmoid_pct,
    validate_common_inputs,
)


def _full_nan() -> dict[str, float]:
    return {
        "prevent_full_10yr_CVD": np.nan,
        "prevent_full_10yr_ASCVD": np.nan,
        "prevent_full_10yr_HF": np.nan,
        "prevent_full_30yr_CVD": np.nan,
        "prevent_full_30yr_ASCVD": np.nan,
        "prevent_full_30yr_HF": np.nan,
    }


def prevent_full(sex, age, tc, hdl, sbp, dm, smoking, bmi, egfr, bptreat, statin, uacr, hba1c, sdi) -> dict[str, float]:
    """R parity: AHAprevent::pred_risk_full (10yr + 30yr)."""
    if not validate_common_inputs(age, sex, sbp, dm, smoking, egfr):
        return _full_nan()
    if invalid_uacr(uacr) or invalid_hba1c(hba1c) or invalid_sdi_decile(sdi):
        return _full_nan()

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
                -3.860385
                + 0.7716794 * ((age - 55) / 10)
                + 0.0062109 * (mmol_conversion(tc - hdl) - 3.5)
                - 0.1547756 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.1933123 * (min(sbp, 110) - 110) / 20
                + 0.3071217 * (max(sbp, 110) - 130) / 20
                + 0.496753 * dm
                + 0.466605 * smoking
                + 0.4780697 * (min(egfr, 60) - 60) / (-15)
                + 0.0529077 * (max(egfr, 60) - 90) / (-15)
                + 0.3034892 * bptreat
                - 0.1556524 * statin
                - 0.0667026 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.1061825 * statin * (mmol_conversion(tc - hdl) - 3.5)
                - 0.0742271 * ((age - 55) / 10) * (mmol_conversion(tc - hdl) - 3.5)
                + 0.0288245 * ((age - 55) / 10) * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.0875188 * ((age - 55) / 10) * (max(sbp, 110) - 130) / 20
                - 0.2267102 * ((age - 55) / 10) * dm
                - 0.0676125 * ((age - 55) / 10) * smoking
                - 0.1493231 * ((age - 55) / 10) * (min(egfr, 60) - 60) / (-15)
                + (
                    0.1361989 * (2 - sdicat(sdi)) * sdicat(sdi)
                    + 0.2261596 * (sdicat(sdi) - 1) * (0.5 * sdicat(sdi))
                    if not pd.isna(sdi)
                    else 0.1804508
                )
                + (0.1645922 * math.log(adjust(uacr)) if not pd.isna(uacr) else 0.0198413)
                + (
                    0.1298513 * (hba1c - 5.3) * dm + 0.1412555 * (hba1c - 5.3) * (1 - dm)
                    if not pd.isna(hba1c)
                    else -0.0031658
                )
            )

            logor_10yr_ASCVD = (
                -4.291503
                + 0.7023067 * ((age - 55) / 10)
                + 0.0898765 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.1407316 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.0256648 * (min(sbp, 110) - 110) / 20
                + 0.314511 * (max(sbp, 110) - 130) / 20
                + 0.4799217 * dm
                + 0.4062049 * smoking
                + 0.3847744 * (min(egfr, 60) - 60) / (-15)
                + 0.0495174 * (max(egfr, 60) - 90) / (-15)
                + 0.2133861 * bptreat
                - 0.0678552 * statin
                - 0.0451416 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.0788187 * statin * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0535985 * ((age - 55) / 10) * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                + 0.0291762 * ((age - 55) / 10) * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.0961839 * ((age - 55) / 10) * (max(sbp, 110) - 130) / 20
                - 0.2001466 * ((age - 55) / 10) * dm
                - 0.0586472 * ((age - 55) / 10) * smoking
                - 0.1537791 * ((age - 55) / 10) * (min(egfr, 60) - 60) / (-15)
                + (
                    0.1413965 * ((2 - sdicat(sdi)) * sdicat(sdi))
                    + 0.228136 * ((sdicat(sdi) - 1) * (0.5 * sdicat(sdi)))
                    if not pd.isna(sdi)
                    else 0.1588908
                )
                + (0.1371824 * math.log(adjust(uacr)) if not pd.isna(uacr) else 0.0061613)
                + (
                    0.123192 * (hba1c - 5.3) * dm + 0.1410572 * (hba1c - 5.3) * (1 - dm)
                    if not pd.isna(hba1c)
                    else 0.005866
                )
            )

            logor_30yr_CVD = (
                -1.748475
                + 0.5073749 * ((age - 55) / 10)
                - 0.0981751 * (((age - 55) / 10) ** 2)
                + 0.0162303 * (mmol_conversion(tc - hdl) - 3.5)
                - 0.1617147 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.1111241 * (min(sbp, 110) - 110) / 20
                + 0.282946 * (max(sbp, 110) - 130) / 20
                + 0.4004069 * dm
                + 0.2918701 * smoking
                + 0.1017102 * (min(egfr, 60) - 60) / (-15)
                + 0.0622643 * (max(egfr, 60) - 90) / (-15)
                + 0.2872416 * bptreat
                - 0.0768135 * statin
                - 0.0557282 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.0917585 * statin * (mmol_conversion(tc - hdl) - 3.5)
                - 0.0679131 * (age - 55) / 10 * (mmol_conversion(tc - hdl) - 3.5)
                + 0.029076 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.0907755 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.2702118 * (age - 55) / 10 * dm
                - 0.1373216 * (age - 55) / 10 * smoking
                - 0.1255864 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (
                    0.1067741 * (2 - sdicat(sdi)) * sdicat(sdi)
                    + 0.1853138 * (sdicat(sdi) - 1) * (0.5 * sdicat(sdi))
                    if not pd.isna(sdi)
                    else 0.1567115
                )
                + (0.1028065 * math.log(adjust(uacr)) if not pd.isna(uacr) else -0.0006181)
                + (
                    0.0925285 * (hba1c - 5.3) * dm + 0.0975598 * (hba1c - 5.3) * (1 - dm)
                    if not pd.isna(hba1c)
                    else 0.0101713
                )
            )

            logor_30yr_ASCVD = (
                -2.314066
                + 0.4386739 * ((age - 55) / 10)
                - 0.0921956 * (((age - 55) / 10) ** 2)
                + 0.0977728 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.1453525 * (mmol_conversion(hdl) - 1.3) / 0.3
                + 0.0590925 * (min(sbp, 110) - 110) / 20
                + 0.2862862 * (max(sbp, 110) - 130) / 20
                + 0.3669136 * dm
                + 0.2354695 * smoking
                + 0.0354338 * (min(egfr, 60) - 60) / (-15)
                + 0.0573093 * (max(egfr, 60) - 90) / (-15)
                + 0.1840085 * bptreat
                + 0.0117504 * statin
                - 0.0331945 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.0664311 * statin * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0492826 * (age - 55) / 10 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                + 0.0288888 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.0964709 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.2279648 * (age - 55) / 10 * dm
                - 0.120405 * (age - 55) / 10 * smoking
                - 0.1157635 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (
                    0.1107632 * ((2 - sdicat(sdi)) * sdicat(sdi))
                    + 0.1840367 * ((sdicat(sdi) - 1) * (0.5 * sdicat(sdi)))
                    if not pd.isna(sdi)
                    else 0.1308962
                )
                + (0.0810739 * math.log(adjust(uacr)) if not pd.isna(uacr) else -0.0147785)
                + (
                    0.0794709 * (hba1c - 5.3) * dm + 0.1002615 * (hba1c - 5.3) * (1 - dm)
                    if not pd.isna(hba1c)
                    else 0.017301
                )
            )

        if can_hf:
            logor_10yr_HF = (
                -4.896524
                + 0.884209 * ((age - 55) / 10)
                - 0.421474 * (min(sbp, 110) - 110) / 20
                + 0.3002919 * (max(sbp, 110) - 130) / 20
                + 0.6170359 * dm
                + 0.5380269 * smoking
                - 0.0191335 * (min(bmi, 30) - 25) / 5
                + 0.2764302 * (max(bmi, 30) - 30) / 5
                + 0.5975847 * (min(egfr, 60) - 60) / (-15)
                + 0.0654197 * (max(egfr, 60) - 90) / (-15)
                + 0.3313614 * bptreat
                - 0.1002304 * bptreat * (max(sbp, 110) - 130) / 20
                - 0.0845363 * ((age - 55) / 10) * (max(sbp, 110) - 130) / 20
                - 0.2989062 * ((age - 55) / 10) * dm
                - 0.1111354 * ((age - 55) / 10) * smoking
                + 0.0008104 * ((age - 55) / 10) * (max(bmi, 30) - 30) / 5
                - 0.1666635 * ((age - 55) / 10) * (min(egfr, 60) - 60) / (-15)
                + (
                    0.1213034 * (2 - sdicat(sdi)) * sdicat(sdi)
                    + 0.2314147 * (sdicat(sdi) - 1) * (0.5 * sdicat(sdi))
                    if not pd.isna(sdi)
                    else 0.1819138
                )
                + (0.1948135 * math.log(adjust(uacr)) if not pd.isna(uacr) else 0.0395368)
                + (
                    0.176668 * (hba1c - 5.3) * dm + 0.1614911 * (hba1c - 5.3) * (1 - dm)
                    if not pd.isna(hba1c)
                    else -0.0010583
                )
            )

            logor_30yr_HF = (
                -2.642208
                + 0.5927507 * ((age - 55) / 10)
                - 0.1028754 * (((age - 55) / 10) ** 2)
                - 0.3593781 * (min(sbp, 110) - 110) / 20
                + 0.2628556 * (max(sbp, 110) - 130) / 20
                + 0.5113472 * dm
                + 0.347344 * smoking
                + 0.0564656 * (min(bmi, 30) - 25) / 5
                + 0.2363857 * (max(bmi, 30) - 30) / 5
                + 0.1971295 * (min(egfr, 60) - 60) / (-15)
                + 0.0735227 * (max(egfr, 60) - 90) / (-15)
                + 0.3219386 * bptreat
                - 0.0880321 * bptreat * (max(sbp, 110) - 130) / 20
                - 0.0863132 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.3425359 * (age - 55) / 10 * dm
                - 0.181405 * (age - 55) / 10 * smoking
                + 0.0031285 * (age - 55) / 10 * (max(bmi, 30) - 30) / 5
                - 0.1356989 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (
                    0.0847634 * ((2 - sdicat(sdi)) * sdicat(sdi))
                    + 0.18397 * ((sdicat(sdi) - 1) * (0.5 * sdicat(sdi)))
                    if not pd.isna(sdi)
                    else 0.1485802
                )
                + (0.1273306 * math.log(adjust(uacr)) if not pd.isna(uacr) else 0.0167008)
                + (
                    0.1378342 * (hba1c - 5.3) * dm + 0.1138832 * (hba1c - 5.3) * (1 - dm)
                    if not pd.isna(hba1c)
                    else 0.0138979
                )
            )
    else:  # male
        if can_cvd_ascvd:
            logor_10yr_CVD = (
                -3.631387
                + 0.7847578 * ((age - 55) / 10)
                + 0.0534485 * (mmol_conversion(tc - hdl) - 3.5)
                - 0.0911282 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.4921973 * (min(sbp, 110) - 110) / 20
                + 0.2972415 * (max(sbp, 110) - 130) / 20
                + 0.4527054 * dm
                + 0.3726641 * smoking
                + 0.3886854 * (min(egfr, 60) - 60) / (-15)
                + 0.0081661 * (max(egfr, 60) - 90) / (-15)
                + 0.2508052 * bptreat
                - 0.1538484 * statin
                - 0.0474695 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.1415382 * statin * (mmol_conversion(tc - hdl) - 3.5)
                - 0.0436455 * ((age - 55) / 10) * (mmol_conversion(tc - hdl) - 3.5)
                + 0.0199549 * ((age - 55) / 10) * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.1022686 * ((age - 55) / 10) * (max(sbp, 110) - 130) / 20
                - 0.1762507 * ((age - 55) / 10) * dm
                - 0.0715873 * ((age - 55) / 10) * smoking
                - 0.1428668 * ((age - 55) / 10) * (min(egfr, 60) - 60) / (-15)
                + (
                    0.0802431 * (2 - sdicat(sdi)) * sdicat(sdi)
                    + 0.275073 * (sdicat(sdi) - 1) * (0.5 * sdicat(sdi))
                    if not pd.isna(sdi)
                    else 0.144759
                )
                + (0.1772853 * math.log(adjust(uacr)) if not pd.isna(uacr) else 0.1095674)
                + (
                    0.1165698 * (hba1c - 5.3) * dm + 0.1048297 * (hba1c - 5.3) * (1 - dm)
                    if not pd.isna(hba1c)
                    else -0.0230072
                )
            )

            logor_10yr_ASCVD = (
                -3.969788
                + 0.7128741 * ((age - 55) / 10)
                + 0.1465201 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.1125794 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.3387216 * (min(sbp, 110) - 110) / 20
                + 0.2980252 * (max(sbp, 110) - 130) / 20
                + 0.399583 * dm
                + 0.3379111 * smoking
                + 0.2582604 * (min(egfr, 60) - 60) / (-15)
                + 0.0147769 * (max(egfr, 60) - 90) / (-15)
                + 0.1686621 * bptreat
                - 0.1073619 * statin
                - 0.0381038 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.1034169 * statin * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0228755 * ((age - 55) / 10) * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                + 0.0267453 * ((age - 55) / 10) * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.0897449 * ((age - 55) / 10) * (max(sbp, 110) - 130) / 20
                - 0.1497464 * ((age - 55) / 10) * dm
                - 0.077206 * ((age - 55) / 10) * smoking
                - 0.1198368 * ((age - 55) / 10) * (min(egfr, 60) - 60) / (-15)
                + (
                    0.0651121 * ((2 - sdicat(sdi)) * sdicat(sdi))
                    + 0.2676683 * ((sdicat(sdi) - 1) * (0.5 * sdicat(sdi)))
                    if not pd.isna(sdi)
                    else 0.1388492
                )
                + (0.1375837 * math.log(adjust(uacr)) if not pd.isna(uacr) else 0.0652944)
                + (
                    0.101282 * (hba1c - 5.3) * dm + 0.1092726 * (hba1c - 5.3) * (1 - dm)
                    if not pd.isna(hba1c)
                    else -0.0112852
                )
            )

            logor_30yr_CVD = (
                -1.504558
                + 0.4427595 * ((age - 55) / 10)
                - 0.1064108 * (((age - 55) / 10) ** 2)
                + 0.0629381 * (mmol_conversion(tc - hdl) - 3.5)
                - 0.1015427 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.2542326 * (min(sbp, 110) - 110) / 20
                + 0.2549679 * (max(sbp, 110) - 130) / 20
                + 0.333835 * dm
                + 0.1873833 * smoking
                + 0.0246102 * (min(egfr, 60) - 60) / (-15)
                + 0.0552014 * (max(egfr, 60) - 90) / (-15)
                + 0.1979729 * bptreat
                - 0.0407714 * statin
                - 0.0365522 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.1232822 * statin * (mmol_conversion(tc - hdl) - 3.5)
                - 0.0441334 * (age - 55) / 10 * (mmol_conversion(tc - hdl) - 3.5)
                + 0.0177865 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.1046657 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.2116113 * (age - 55) / 10 * dm
                - 0.1277905 * (age - 55) / 10 * smoking
                - 0.0955922 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (
                    0.0256704 * (2 - sdicat(sdi)) * sdicat(sdi)
                    + 0.1887637 * (sdicat(sdi) - 1) * (0.5 * sdicat(sdi))
                    if not pd.isna(sdi)
                    else 0.089241
                )
                + (0.0894596 * math.log(adjust(uacr)) if not pd.isna(uacr) else 0.0710124)
                + (
                    0.0676202 * (hba1c - 5.3) * dm + 0.063409 * (hba1c - 5.3) * (1 - dm)
                    if not pd.isna(hba1c)
                    else 0.0038783
                )
            )

            logor_30yr_ASCVD = (
                -1.985368
                + 0.3743566 * ((age - 55) / 10)
                - 0.0995499 * (((age - 55) / 10) ** 2)
                + 0.1544808 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.1215297 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.1083968 * (min(sbp, 110) - 110) / 20
                + 0.2555179 * (max(sbp, 110) - 130) / 20
                + 0.2696998 * dm
                + 0.1628432 * smoking
                - 0.077507 * (min(egfr, 60) - 60) / (-15)
                + 0.0583407 * (max(egfr, 60) - 90) / (-15)
                + 0.1120322 * bptreat
                - 0.0025063 * statin
                - 0.0256116 * bptreat * (max(sbp, 110) - 130) / 20
                + 0.0886745 * statin * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                - 0.0254507 * (age - 55) / 10 * ((mmol_conversion(tc) - mmol_conversion(hdl)) - 3.5)
                + 0.0244639 * (age - 55) / 10 * (mmol_conversion(hdl) - 1.3) / 0.3
                - 0.0869146 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.165745 * (age - 55) / 10 * dm
                - 0.1244714 * (age - 55) / 10 * smoking
                - 0.0624552 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (
                    0.015675 * (2 - sdicat(sdi)) * sdicat(sdi)
                    + 0.1864231 * (sdicat(sdi) - 1) * (0.5 * sdicat(sdi))
                    if not pd.isna(sdi)
                    else 0.0845697
                )
                + (0.0560171 * math.log(adjust(uacr)) if not pd.isna(uacr) else 0.0252244)
                + (
                    0.0501422 * (hba1c - 5.3) * dm + 0.0722905 * (hba1c - 5.3) * (1 - dm)
                    if not pd.isna(hba1c)
                    else 0.0114945
                )
            )

        if can_hf:
            logor_10yr_HF = (
                -4.663513
                + 0.9095703 * ((age - 55) / 10)
                - 0.6765184 * (min(sbp, 110) - 110) / 20
                + 0.3111651 * (max(sbp, 110) - 130) / 20
                + 0.5535052 * dm
                + 0.4326811 * smoking
                - 0.0854286 * (min(bmi, 30) - 25) / 5
                + 0.3551736 * (max(bmi, 30) - 30) / 5
                + 0.5102245 * (min(egfr, 60) - 60) / (-15)
                + 0.015472 * (max(egfr, 60) - 90) / (-15)
                + 0.2570964 * bptreat
                - 0.0591177 * bptreat * (max(sbp, 110) - 130) / 20
                - 0.1219056 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.2437577 * (age - 55) / 10 * dm
                - 0.105363 * (age - 55) / 10 * smoking
                + 0.0037907 * (age - 55) / 10 * (max(bmi, 30) - 30) / 5
                - 0.1660207 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (
                    0.1106372 * ((2 - sdicat(sdi)) * sdicat(sdi))
                    + 0.3371204 * ((sdicat(sdi) - 1) * (0.5 * sdicat(sdi)))
                    if not pd.isna(sdi)
                    else 0.1694628
                )
                + (0.2164607 * math.log(adjust(uacr)) if not pd.isna(uacr) else 0.1702805)
                + (
                    0.148297 * (hba1c - 5.3) * dm + 0.1234088 * (hba1c - 5.3) * (1 - dm)
                    if not pd.isna(hba1c)
                    else -0.0234637
                )
            )

            logor_30yr_HF = (
                -2.425439
                + 0.5478829 * ((age - 55) / 10)
                - 0.1111928 * (((age - 55) / 10) ** 2)
                - 0.4547346 * (min(sbp, 110) - 110) / 20
                + 0.2527602 * (max(sbp, 110) - 130) / 20
                + 0.4385384 * dm
                + 0.2397952 * smoking
                + 0.0640931 * (min(bmi, 30) - 25) / 5
                + 0.2643081 * (max(bmi, 30) - 30) / 5
                + 0.1354588 * (min(egfr, 60) - 60) / (-15)
                + 0.0570689 * (max(egfr, 60) - 90) / (-15)
                + 0.220666 * bptreat
                - 0.0436769 * bptreat * (max(sbp, 110) - 130) / 20
                - 0.1168376 * (age - 55) / 10 * (max(sbp, 110) - 130) / 20
                - 0.2730055 * (age - 55) / 10 * dm
                - 0.1573691 * (age - 55) / 10 * smoking
                - 0.0174998 * (age - 55) / 10 * (max(bmi, 30) - 30) / 5
                - 0.1128676 * (age - 55) / 10 * (min(egfr, 60) - 60) / (-15)
                + (
                    0.057746 * ((2 - sdicat(sdi)) * sdicat(sdi))
                    + 0.2446441 * ((sdicat(sdi) - 1) * (0.5 * sdicat(sdi)))
                    if not pd.isna(sdi)
                    else 0.1076782
                )
                + (0.1233486 * math.log(adjust(uacr)) if not pd.isna(uacr) else 0.1274796)
                + (
                    0.0985062 * (hba1c - 5.3) * dm + 0.0804844 * (hba1c - 5.3) * (1 - dm)
                    if not pd.isna(hba1c)
                    else 0.0022806
                )
            )

    return {
        "prevent_full_10yr_CVD": sigmoid_pct(logor_10yr_CVD) if not pd.isna(logor_10yr_CVD) else np.nan,
        "prevent_full_10yr_ASCVD": sigmoid_pct(logor_10yr_ASCVD) if not pd.isna(logor_10yr_ASCVD) else np.nan,
        "prevent_full_10yr_HF": sigmoid_pct(logor_10yr_HF) if not pd.isna(logor_10yr_HF) else np.nan,
        "prevent_full_30yr_CVD": sigmoid_pct(logor_30yr_CVD) if not pd.isna(logor_30yr_CVD) else np.nan,
        "prevent_full_30yr_ASCVD": sigmoid_pct(logor_30yr_ASCVD) if not pd.isna(logor_30yr_ASCVD) else np.nan,
        "prevent_full_30yr_HF": sigmoid_pct(logor_30yr_HF) if not pd.isna(logor_30yr_HF) else np.nan,
    }

