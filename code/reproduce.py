from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.special import expit, logit
from scipy.stats import t
ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / 'data'
sys.path.insert(0, str(ROOT / 'third_party'))
from prevent._uacr import prevent_uacr
LABELS = ['<3%', '3–<5%', '5–<10%', '>=10%']
THRESHOLDS = np.array([3.0, 5.0, 10.0])

def survey_mean(all_data, values, domain, weight='WTMEC2YR', proportion=False):
    z = all_data[['SDMVSTRA', 'SDMVPSU', weight]].copy()
    values = pd.Series(values, index=all_data.index, dtype=float)
    domain = pd.Series(domain, index=all_data.index).fillna(False).astype(bool)
    valid_design = z[weight].gt(0) & z[['SDMVSTRA', 'SDMVPSU']].notna().all(axis=1)
    use = domain & valid_design & values.notna()
    assert use.sum() > 0
    weights = z.loc[use, weight].to_numpy()
    y = values.loc[use].to_numpy()
    denominator = weights.sum()
    point = float(np.dot(weights, y) / denominator)
    z['linearized'] = 0.0
    z.loc[use, 'linearized'] = weights * (y - point) / denominator
    clusters = z.loc[valid_design].groupby(['SDMVSTRA', 'SDMVPSU']).linearized.sum()
    variance = 0.0
    paired_variance = 0.0
    all_paired = True
    for _, vals in clusters.groupby(level=0):
        m = len(vals)
        assert m > 1, 'Lonely PSU requires an explicit variance policy'
        variance += m / (m - 1) * float(np.sum((vals.to_numpy() - vals.mean()) ** 2))
        if m == 2:
            paired_variance += float((vals.iloc[0] - vals.iloc[1]) ** 2)
        else:
            all_paired = False
    if all_paired:
        assert np.isclose(variance, paired_variance, atol=1e-15, rtol=1e-12)
    represented = z.loc[use, ['SDMVSTRA', 'SDMVPSU']].drop_duplicates()
    df = len(represented) - represented.SDMVSTRA.nunique()
    assert df > 0
    se = float(np.sqrt(variance))
    critical = float(t.ppf(0.975, df))
    if proportion and 0 < point < 1:
        half = critical * se / (point * (1 - point))
        lower, upper = expit([logit(point) - half, logit(point) + half])
        ci_method = 'Taylor-logit, t critical value'
    elif proportion:
        lower = upper = None
        ci_method = 'boundary proportion: interval not estimated; do not infer certainty'
    else:
        lower, upper = (point - critical * se, point + critical * se)
        ci_method = 'Taylor-linear, t critical value'
    return dict(n=int(use.sum()), events=int(y.sum()) if proportion else None, estimate=point, se=se, lower95=None if lower is None else float(lower), upper95=None if upper is None else float(upper), df=int(df), strata=int(represented.SDMVSTRA.nunique()), psus=len(represented), sum_weights=float(denominator), ci_method=ci_method)

def weighted_quantile(values, weights, probabilities=(0.1, 0.5, 0.9)):
    order = np.argsort(values)
    v, w = (np.asarray(values)[order], np.asarray(weights)[order])
    cw = np.cumsum(w) / w.sum()
    return {str(q): float(v[min(np.searchsorted(cw, q), len(v) - 1)]) for q in probabilities}

def core_summary(d, s, weight='WTMEC2YR'):
    dom = d.index.isin(s.index)

    def prop(values):
        result = survey_mean(d, values.reindex(d.index), dom, weight, True)
        event_units = s.loc[values.astype(bool), ['SDMVSTRA', 'SDMVPSU']].drop_duplicates()
        result['psus_with_events'] = len(event_units)
        return result

    def mean(values):
        return survey_mean(d, values.reindex(d.index), dom, weight)
    a, b = (s.category_first, s.category_second)
    result = {'any_category_change': prop(a.ne(b)), 'lower_category': prop(b.lt(a)), 'higher_category': prop(b.gt(a)), 'mean_delta_pp': mean(s.delta_pp), 'mean_absolute_delta_pp': mean(s.delta_pp.abs()), 'absolute_delta_ge1pp': prop(s.delta_pp.abs().ge(1)), 'absolute_delta_quantiles_pp': weighted_quantile(s.delta_pp.abs(), s[weight]), 'thresholds': {}, 'weight_variable': weight, 'category_labels': LABELS, 'maximum_categories_moved': int((b - a).abs().max()), 'male_n': int(s.sex.eq(0).sum()), 'female_n': int(s.sex.eq(1).sum())}
    for thr in THRESHOLDS:
        result['thresholds'][str(int(thr))] = prop(s.risk_first.ge(thr).ne(s.risk_second.ge(thr)))
    counts = pd.crosstab(a, b).reindex(index=range(4), columns=range(4), fill_value=0)
    pct = np.zeros((4, 4))
    for i in range(4):
        for j in range(4):
            pct[i, j] = 100 * s.loc[a.eq(i) & b.eq(j), weight].sum() / s[weight].sum()
    assert counts.to_numpy().sum() == len(s) and np.isclose(pct.sum(), 100)
    result['transition_counts'] = counts.to_numpy().tolist()
    result['transition_percent'] = pct.tolist()
    return result

def response_features(d):
    age = (d.age - 55) / 10
    return pd.DataFrame({'intercept': 1.0, 'age10': age, 'age10_squared': age ** 2, 'female': d.sex, 'smoking': d.smoking, 'diabetes': d.dm, 'bp_treatment': d.bptreat, 'statin': d.statin, 'egfr15': (d.egfr - 90) / 15, 'sbp20': (d.sbp - 130) / 20, 'nonhdl40': (d.tc - d.hdl - 140) / 40, 'hdl10': (d.hdl - 50) / 10, 'log_random_uacr': np.log(d.URDACT)}, index=d.index)

def fit_response(x, r, weights):
    fit = sm.GLM(r, x, family=sm.families.Binomial(), freq_weights=weights / weights.mean()).fit(maxiter=100, tol=1e-10)
    assert fit.converged, 'Response model did not converge'
    prob = np.asarray(fit.predict(x))
    assert np.isfinite(prob).all() and ((prob > 0) & (prob < 1)).all()
    return (fit, prob)

def effective_n(w):
    return float(w.sum() ** 2 / np.dot(w, w))

def response_analysis(d, target, paired):
    z = d.loc[target].copy()
    r = z.index.isin(paired.index).astype(float)
    assert int(r.sum()) == len(paired)
    x = response_features(z)
    w = z.WTMEC2YR.to_numpy()
    y = paired.category_random.ne(paired.category_morning).reindex(z.index).fillna(False).to_numpy(float)
    fit, prob = fit_response(x, r, w)
    wr = w * r / prob
    estimate = float(np.dot(wr, y) / wr.sum())
    unadjusted = float(np.dot(w * r, y) / (w * r).sum())
    const = w * r / np.average(r, weights=w)
    assert np.isclose(np.dot(const, y) / const.sum(), unadjusted)
    checks = []
    for name in x.columns[1:]:
        v = x[name].to_numpy()
        tm = np.average(v, weights=w)
        scale = np.sqrt(np.average((v - tm) ** 2, weights=w))
        checks.append(dict(variable=name, target_mean=float(tm), responders_mean=float(np.average(v, weights=w * r)), ipw_mean=float(np.average(v, weights=wr)), standardized_difference_before=float((np.average(v, weights=w * r) - tm) / scale), standardized_difference_after=float((np.average(v, weights=wr) - tm) / scale)))
    design = d.loc[d.WTMEC2YR.gt(0), ['SDMVSTRA', 'SDMVPSU']].drop_duplicates()
    replicates = []
    for h, units in design.groupby('SDMVSTRA'):
        psus = units.SDMVPSU.to_numpy()
        m = len(psus)
        assert m >= 2
        in_h = z.SDMVSTRA.eq(h).to_numpy()
        for psu in psus:
            mult = np.ones(len(z))
            mult[in_h] = m / (m - 1)
            mult[in_h & z.SDMVPSU.eq(psu).to_numpy()] = 0.0
            rw = w * mult
            _, rp = fit_response(x, r, rw)
            adjusted = rw * r / rp
            theta = float(np.dot(adjusted, y) / adjusted.sum())
            raw = float(np.dot(rw * r, y) / (rw * r).sum())
            replicates.append(dict(stratum=int(h), deleted_psu=int(psu), factor=(m - 1) / m, response_weighted=theta, unadjusted=raw, converged=True))
    df = len(design) - design.SDMVSTRA.nunique()
    var = sum((v['factor'] * (v['response_weighted'] - estimate) ** 2 for v in replicates))
    var_raw = sum((v['factor'] * (v['unadjusted'] - unadjusted) ** 2 for v in replicates))
    se = float(np.sqrt(var))
    interval = expit(logit(estimate) + np.array([-1, 1]) * t.ppf(0.975, df) * se / (estimate * (1 - estimate)))
    return_rate = float(np.average(r, weights=w))
    observed_changes_target = float(np.dot(w * r, y) / w.sum())
    minimum = float(prob.min())
    maximum_multiplier = float((1 / prob[r == 1]).max())
    result = dict(target_n=len(z), responders_n=int(r.sum()), nonresponders_n=int((1 - r).sum()), weighted_return_rate=return_rate, unadjusted=unadjusted, ipw_estimate=estimate, ipw_se_jkn=se, ipw_lower95=float(interval[0]), ipw_upper95=float(interval[1]), df=int(df), replicates=len(replicates), all_replicates_converged=True, unadjusted_se_jkn=float(np.sqrt(var_raw)), response_probability_quantiles={str(q): float(np.quantile(prob, q)) for q in [0, 0.01, 0.1, 0.5, 0.9, 0.99, 1]}, maximum_inverse_response_multiplier=maximum_multiplier, weight_stability_gate_passed=minimum >= 0.1 and maximum_multiplier <= 10, kish_effective_n_before=effective_n(w * r), kish_effective_n_after=effective_n(wr), extreme_missing_outcome_bounds=[observed_changes_target, observed_changes_target + 1 - return_rate], bound_interpretation='Identification extremes, not confidence intervals; excludes missing baseline predictors', response_coefficients={k: float(v) for k, v in zip(x.columns, fit.params)}, balance=checks, max_absolute_standardized_difference_before=max((abs(v['standardized_difference_before']) for v in checks)), max_absolute_standardized_difference_after=max((abs(v['standardized_difference_after']) for v in checks)), limitation='IPW assumes response exchangeability conditional on the chosen measured covariates and a sufficient response model; agreement does not establish absence of unmeasured selection')
    prediction = pd.DataFrame({'SEQN': z.SEQN, 'response': r, 'response_probability': prob, 'mec_weight': w, 'ipw_weight_responders': wr, 'any_change_if_observed': np.where(r == 1, y, np.nan)})
    prediction.to_csv(OUT / 'response_model_predictions.csv', index=False)
    pd.DataFrame(replicates).to_csv(OUT / 'jackknife_replicates.csv', index=False)
    pd.DataFrame(checks).to_csv(OUT / 'response_balance.csv', index=False)
    return result

def interpretation_summary(d, domain, first, second, weight, prefix, scenarios=False):
    domain = pd.Series(domain, index=d.index).astype(bool)
    changed = d.category_first.ne(d.category_second) & domain
    down = d.category_second.lt(d.category_first) & domain
    up = d.category_second.gt(d.category_first) & domain
    high1, high2 = (d[first].ge(30), d[second].ge(30))
    flag_change = high1.ne(high2)
    both_low = ~high1 & ~high2
    assert (down | up).equals(changed)
    assert (high2.ne(high1) & domain).equals(flag_change & domain)
    assert (~high2 & ~high1 & domain).equals(both_low & domain)
    assert (d.category_second.ne(d.category_first) & domain).equals(changed)

    def prop(label, values, within=None):
        use = domain if within is None else domain & within
        if not use.any():
            return {'n': 0, 'events': 0, 'estimate': None, 'status': 'undefined empty denominator'}
        v = pd.Series(values, index=d.index).astype(float)
        result = survey_mean(d, v, use, weight, True)
        constant = v.loc[use].nunique() == 1
        if constant:
            result.update(estimate=float(v.loc[use].iloc[0]), se=0.0, lower95=None, upper95=None, ci_method='constant indicator; no interval estimated, do not infer certainty')
        exported = pd.read_csv(FROZEN / ('primary_participant_results.csv' if prefix == 'primary' else 'external_bridged_participant_results.csv'), float_precision='round_trip').set_index('SEQN')
        ids = d.loc[use, 'SEQN']
        numerator = float(np.dot(exported.loc[ids, weight], v.loc[use]))
        denominator = float(exported.loc[ids, weight].sum())
        assert abs(numerator / denominator - result['estimate']) < 1e-13
        result['numerator_weight'] = numerator
        result['denominator_weight'] = denominator
        result['psus_with_numerator'] = len(d.loc[use & v.eq(1), ['SDMVSTRA', 'SDMVPSU']].drop_duplicates())
        return result
    result = {'n': int(domain.sum()), 'change_n': int(changed.sum()), 'weight': weight, 'overall_disagreement': prop('overall', changed), 'joint_flag_risk': []}
    for flag in [False, True]:
        for risk in [False, True]:
            values = flag_change.eq(flag) & changed.eq(risk)
            result['joint_flag_risk'].append({'uacr_flag_changed': flag, 'ascvd_category_changed': risk, 'population_fraction': prop(f'joint_{flag}_{risk}', values)})
    assert sum((x['population_fraction']['events'] for x in result['joint_flag_risk'])) == int(domain.sum())
    assert abs(sum((x['population_fraction']['estimate'] for x in result['joint_flag_risk'])) - 1) < 1e-13
    result['same_uacr_flag_share_of_changes'] = prop('same_flag_share', ~flag_change, changed)
    result['both_uacr_below30'] = {'population_fraction': prop('both_low_population', both_low), 'disagreement_within': prop('both_low_disagreement', changed, both_low), 'share_of_all_changes': prop('both_low_share_changes', both_low, changed), 'down_n': int((domain & both_low & down).sum()), 'up_n': int((domain & both_low & up).sum())}
    if scenarios:
        result['selection_scenarios'] = {}
        for anchor, uacr in [('first', first), ('second', second)]:
            risk = d['risk_' + anchor]
            near = pd.concat([(risk - thr).abs().le(0.5) for thr in [3.0, 5.0, 10.0]], axis=1).any(axis=1)
            selectors = {'everyone': domain, 'uacr_ge30': d[uacr].ge(30), 'risk_within_half_pp': near}
            rows = []
            for name, selected in selectors.items():
                key = anchor + '/' + name
                row = {'scenario': name, 'selected_n': int((selected & domain).sum()), 'selected_change_n': int((selected & changed).sum()), 'missed_change_n': int((~selected & changed).sum()), 'workload': prop(key + '/workload', selected), 'coverage': prop(key + '/coverage', selected, changed), 'yield': prop(key + '/yield', changed, selected), 'missed_fraction_of_changes': prop(key + '/missed', ~selected, changed), 'selected_disagreement_population': prop(key + '/selected_changes', selected & changed)}
                assert row['selected_change_n'] + row['missed_change_n'] == int(changed.sum())
                assert abs(row['coverage']['estimate'] + row['missed_fraction_of_changes']['estimate'] - 1) < 1e-13
                a = row['workload']['estimate'] * row['yield']['estimate']
                b = row['coverage']['estimate'] * result['overall_disagreement']['estimate']
                assert abs(a - b) < 1e-13
                rows.append(row)
                d[anchor + '_' + name] = (selected & domain).astype(int)
            result['selection_scenarios'][anchor] = rows
    fields = ['SEQN', first, second, 'risk_first', 'risk_second', 'category_first', 'category_second', 'SDMVSTRA', 'SDMVPSU', weight]
    if scenarios:
        fields += [a + '_' + n for a in ['first', 'second'] for n in ['everyone', 'uacr_ge30', 'risk_within_half_pp']]
    export = d.loc[domain, fields].copy()
    export['uacr_flag_changed'] = flag_change.loc[domain].astype(int)
    export['both_uacr_below30'] = both_low.loc[domain].astype(int)
    export['ascvd_category_changed'] = changed.loc[domain].astype(int)
    export.to_csv(OUT / (prefix + '_interpretability_participants.csv'), index=False)
    return result

def csv(name):
    return pd.read_csv(FROZEN / name, float_precision='round_trip')

def design_frame(s, kind):
    units = csv(kind + '_design.csv')
    return (pd.concat([s, units], ignore_index=True, sort=False), np.arange(len(s)))

def rescore(s, first, second):
    errors = []
    for label, urine in [('first', first), ('second', second)]:
        risks = []
        for _, row in s.iterrows():
            x = {k: float(row[k]) for k in ['sex', 'age', 'tc', 'hdl', 'sbp', 'dm', 'smoking', 'bmi', 'egfr', 'bptreat', 'statin']}
            x['uacr'] = float(row[urine])
            risks.append(prevent_uacr(**x)['prevent_uacr_10yr_ASCVD'])
        error = np.abs(np.asarray(risks) - s['risk_' + label].to_numpy())
        assert np.isfinite(error).all() and error.max() < 1e-10
        assert np.array_equal(np.searchsorted(THRESHOLDS, risks, side='right'), s['category_' + label])
        errors.append(float(error.max()))
    return max(errors)

def main():
    global OUT
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    OUT = args.output.resolve()
    OUT.mkdir(parents=True, exist_ok=False)
    p = csv('primary_participant_results.csv')
    e = csv('external_bridged_participant_results.csv')
    native = csv('external_native.csv')
    alt = csv('alternative_bp.csv')
    errors = {'primary': rescore(p, 'URDACT', 'URDACT2'), 'external_bridged': rescore(e, 'uacr_first_bridged', 'uacr_repeat_bridged'), 'external_native': rescore(native, 'uacr_first', 'uacr_repeat'), 'alternative_bp': rescore(alt, 'URDACT', 'URDACT2')}
    d, idx = design_frame(p, 'primary')
    primary = d.loc[idx]
    z, ei = design_frame(e, 'external')
    external = z.loc[ei]
    totals = {'primary': core_summary(d, primary), 'external_bridged': core_summary(z, external, 'unit_weight')}
    sens = {}
    for name, selection in [('exclude_indeterminate_pregnancy', ~primary.pregnancy_indeterminate.astype(bool)), ('exclude_albumin_below_lod', primary.URXUMA.ge(0.3) & primary.URXUMA2.ge(0.3)), ('matched_age30_69', primary.age.le(69)), ('fasting_treatment_context', primary.WTSAF2YR.gt(0) & primary.LBDLDL.between(70, 189) & primary.dm.eq(0) & primary.statin.eq(0))]:
        sens[name] = core_summary(d, primary.loc[selection], 'WTSAF2YR' if name == 'fasting_treatment_context' else 'WTMEC2YR')
    ad, ai = design_frame(alt, 'primary')
    sens['alternative_sbp'] = core_summary(ad, ad.loc[ai])
    nd, ni = design_frame(native, 'external')
    totals.update(primary_sensitivities=sens, external_sensitivities={'native_assay': core_summary(nd, nd.loc[ni], 'unit_weight'), 'exclude_indeterminate_pregnancy': core_summary(z, external.loc[~external.pregnancy_indeterminate.astype(bool)], 'unit_weight')})
    target = csv('response_target.csv')
    target = target.merge(p[['SEQN', 'category_first', 'category_second']], on='SEQN', how='left', validate='one_to_one')
    td, ti = design_frame(target, 'primary')
    tm = td.index.isin(ti)
    observed = td.loc[tm & td.responded.eq(1)].rename(columns={'category_first': 'category_random', 'category_second': 'category_morning'})
    response = response_analysis(td, tm, observed)
    totals['primary_response_sensitivity'] = response
    ep = interpretation_summary(d, d.index.isin(idx), 'URDACT', 'URDACT2', 'WTMEC2YR', 'primary', True)
    ee = interpretation_summary(z, z.index.isin(ei), 'uacr_first_bridged', 'uacr_repeat_bridged', 'unit_weight', 'external')
    use = d.index.isin(idx) & d.WTSAF2YR.gt(0) & d.LBDLDL.between(70, 189) & d.dm.eq(0) & d.statin.eq(0)
    clinical = {'n': int(use.sum()), 'any_category_change': survey_mean(d, d.category_first.ne(d.category_second), use, 'WTSAF2YR', True), 'thresholds': {}}
    for boundary in [3, 5, 10]:
        a = d.risk_first.ge(boundary)
        b = d.risk_second.ge(boundary)
        vals = {'up': ~a & b, 'down': a & ~b, 'gross': a.ne(b), 'above_first': a, 'above_second': b, 'net': b.astype(int) - a.astype(int)}
        clinical['thresholds'][str(boundary)] = {k: survey_mean(d, v, use, 'WTSAF2YR', k != 'net') for k, v in vals.items()}
    for name, value in [('core_recomputed', totals), ('interpretation_recomputed', {'primary': ep, 'external': ee}), ('clinical_context_recomputed', clinical)]:
        (OUT / (name + '.json')).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    status = {'status': 'COMPLETE', 'rescoring_max_error_pp': errors, 'response_refits': response['replicates'], 'all_response_refits_converged': response['all_replicates_converged'], 'scope': 'Re-execution from processed participant inputs, not raw-data cleaning or independent clinical validation'}
    (OUT / 'verification.json').write_text(json.dumps(status, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(status, indent=2))
if __name__ == '__main__':
    main()
