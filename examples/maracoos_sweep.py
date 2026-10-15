"""Sweep of OI settings and QC thresholds for the MARACOOS 6 km totals, scored by withheld-site prediction.

Companion to maracoos_study.py, run on the same archived week:

    python examples/maracoos_sweep.py --data ~/Documents/hfradar_parity/2026-10-06_week

Each configuration combines the production radials of every site except one and predicts the withheld site's
radial velocities (projection of the total on the radial direction), as in maracoos_study.py. Predictions are
stored before the uncertainty and site-share QC, so every QC threshold can be scored afterwards from the same
solutions. Configurations:

    A. OI length scale Sx = Sy (8, 10, 12, 15 km) x variance ratio mdlvar / errvar (12.7, 6.4, 3.2, 1.6)
    B. t-1/t/t+1 radials, unweighted or weighted with an OI time scale T (1, 3, 6, 12 h), for Sx = 10, 12, 15 km,
       and for Sx = 10 km with errvar 132
    C. Each radial's own error variance (ETMP squared) instead of a single errvar, with and without the window

The OI velocity and normalized uncertainty depend on mdlvar and errvar only through their ratio, so mdlvar is fixed
at production's 420 (cm/s)^2 and errvar varies. Production is Sx = 15 km, errvar 66, no window, Uerr/Verr <= 0.6;
Nwankwo et al. (2026) used 625 / 100 (ratio 6.25), Sx = 10 km and the unweighted window.

Settings are selected on the first part of the week (--split) and checked on the rest. Results (tables, figures,
sweep_report.md) are written to --output; finished solutions are cached there and reused.
"""
import argparse
import copy
import json
import logging
import time
import warnings
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from joblib import Parallel, delayed, parallel_config
from scipy.spatial import cKDTree

import realtime_totals as rt
from hfradarpy.totals import Total, combineRadials
from maracoos_study import MAX_DISTANCE_KM, xyz
from maracoos_study_analysis import paired_interval

EXAMPLES = Path(__file__).resolve().parent
MDLVAR = 420.0
PRODUCTION = 'sx15_e66'
MAX_SPEED = 300.0
UNCERTAINTY = np.round(np.arange(0.30, 1.001, 0.05), 2)


def configurations():
    configs = {}
    for sx in (8, 10, 12, 15):
        for errvar in (33, 66, 132, 264):
            configs[f'sx{sx}_e{errvar}'] = dict(sx=float(sx), errvar=float(errvar), window=0, time_scale=None,
                                                errcol=None, family='length scale x variance ratio')
    for sx in (10, 12, 15):
        configs[f'sx{sx}_e66_window'] = dict(sx=float(sx), errvar=66.0, window=1, time_scale=None, errcol=None,
                                             family='window')
        for T in (1, 3, 6, 12):
            configs[f'sx{sx}_e66_window_T{T}'] = dict(sx=float(sx), errvar=66.0, window=1, time_scale=float(T),
                                                      errcol=None, family='window')
    for T in (None, 6):
        configs[f'sx10_e132_window{"_T%d" % T if T else ""}'] = dict(sx=10.0, errvar=132.0, window=1,
                                                                     time_scale=float(T) if T else None,
                                                                     errcol=None, family='window')
    for sx in (10, 15):
        configs[f'sx{sx}_etmp'] = dict(sx=float(sx), errvar=66.0, window=0, time_scale=None, errcol='ETMP',
                                       family='radial error')
    configs['sx10_etmp_window'] = dict(sx=10.0, errvar=66.0, window=1, time_scale=None, errcol='ETMP',
                                       family='radial error')
    return configs


def setup(data):
    data = Path(data).expanduser().resolve()
    doc = rt.load_config('realtime_06km_maracoos_v03', data / 'prod/totals_configs_export.json')
    overrides = json.loads((data / 'overrides.json').read_text())
    driver = copy.deepcopy(rt.DRIVER_DEFAULTS)
    rt.deep_update(driver, overrides.pop('driver', {}))
    rt.deep_update(doc, overrides)
    conf = doc['conf']
    grid = rt.load_grid(conf['Totals']['GridFile'])
    # Production QC305: grid points inside the land polygon fail; computed with the production code path
    T = Total(grid=grid)
    T.initialize_qc()
    T.qc_valid_location_polygon(*rt.load_polygon(conf['QC']['QC305_MaskFile']), flag_inside=True, name='land')
    land = (T.data['Q305'] == 4).to_numpy()
    return conf, driver, grid, land


def load_inputs(conf, hours):
    rconf = conf['Radials']
    sites, jobs = rt.pattern_jobs(rconf)
    types = dict(jobs)['BestChoice']
    mask = rt.load_polygon(rconf.get('MaskDir'))
    inputs = {}
    for t in pd.date_range(hours[0] - pd.Timedelta('1h'), hours[-1] + pd.Timedelta('1h'), freq='h'):
        files = rt.radial_filenames(rconf['BaseDir'], sites, types, t, rconf.get('MonthFlag', False),
                                    rconf.get('TypeFlag', False))
        radials, _, _ = rt.load_radials(files)
        rt.clean_radials(radials, remove_flagged=conf['QC'].get('Process', False), mask=mask)
        for r in radials.values():
            keep = [c for c in ('LOND', 'LATD', 'VELO', 'HEAD', 'ETMP') if c in r.data]
            r.data = r.data[keep].copy()
            for attr in ('_tables', 'diagnostics_radial', 'diagnostics_hardware', 'range_information', 'crad_data'):
                r.__dict__.pop(attr, None)
        inputs[t] = radials
    return inputs


def loso_task(t, site, name, config, radials, grid, land, conf, cache):
    """Predict the withheld site's radials at hour t; store them with the QC quantities of each prediction."""
    warnings.simplefilter('ignore', category=FutureWarning)
    logging.disable(logging.WARNING)
    path = cache / name / f'{t:%Y%m%dT%H%M}_{site}.npz'
    if path.exists():
        return
    target = radials[t][site].data
    dist, point = cKDTree(xyz(grid.x.to_numpy(), grid.y.to_numpy())).query(
        xyz(target.LOND.to_numpy(), target.LATD.to_numpy()))
    eligible = dist <= MAX_DISTANCE_KM
    selected, inverse = np.unique(point[eligible], return_inverse=True)
    solve_grid = grid.iloc[selected].reset_index(drop=True)

    rows, sites = [], []
    for k in range(-config['window'], config['window'] + 1):
        for s, r in radials.get(t + pd.Timedelta(hours=k), {}).items():
            if s != site:
                rows.append(r)
                sites.append(s)
    rdf = pd.DataFrame({'Radial': pd.Series(rows, index=sites, dtype=object)})
    tc = conf['Totals']
    n = len(target)
    out = {k: np.full(n, np.nan) for k in ('pred', 'uerr', 'msfr', 'ror')}
    if rdf.index.nunique() >= tc.get('MinNumSites', 2) and len(solve_grid):
        kwargs = dict(gRes=6000, tStp=t.to_pydatetime(), minContrSites=tc.get('MinNumSites', 2),
                      minContrRads=tc.get('MinNumRads', 3), useLandMask=False, dropPoints=False, engine='vectorized',
                      method='oi', mdlvar=MDLVAR, errvar=config['errvar'], sx=config['sx'], sy=config['sx'],
                      oi_option='exponential', sRad=2 * config['sx'] * 1000, oi_time_scale=config['time_scale'],
                      oi_errvar_column=config['errcol'])
        if config['window']:
            kwargs['tempthreshold'] = pd.Timedelta(hours=config['window'])
        total, _ = combineRadials(rdf, solve_grid, **kwargs)
        d = total.data
        u, v = d.VELU.to_numpy(float)[inverse], d.VELV.to_numpy(float)[inverse]
        # Production QC303 (speed) and QC305 (land); uncertainty and site share are applied when scoring
        ok = np.isfinite(u) & np.isfinite(v) & (np.hypot(u, v) <= MAX_SPEED) & ~land[selected][inverse]
        head = np.deg2rad(target.HEAD.to_numpy(float)[eligible])
        out['pred'][eligible] = np.where(ok, u * np.sin(head) + v * np.cos(head), np.nan)
        out['uerr'][eligible] = np.maximum(d.UQAL.to_numpy(float), d.VQAL.to_numpy(float))[inverse]
        out['msfr'][eligible] = d.MSFR.to_numpy(float)[inverse]
        out['ror'][eligible] = d.ROR.to_numpy(float)[inverse]
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.npz.tmp')
    with tmp.open('wb') as f:
        np.savez_compressed(f, observed=target.VELO.to_numpy(float), eligible=eligible, **out)
    tmp.replace(path)


def compute(configs, inputs, hours, grid, land, conf, cache, jobs):
    order = sorted(configs, key=lambda c: (-configs[c]['window'], -configs[c]['sx']))  # slowest first
    tasks = [(t, site, name) for name in order for t in hours for site in sorted(inputs[t])]
    started = time.monotonic()
    with parallel_config(backend='loky', inner_max_num_threads=1):
        done = Parallel(n_jobs=jobs, return_as='generator_unordered', pre_dispatch=2 * jobs)(
            delayed(loso_task)(t, site, name, configs[name],
                               {t + pd.Timedelta(hours=k): inputs.get(t + pd.Timedelta(hours=k), {})
                                for k in (-1, 0, 1)}, grid, land, conf, cache)
            for t, site, name in tasks)
        for i, _ in enumerate(done, 1):
            if i % 200 == 0 or i == len(tasks):
                print(f'{i}/{len(tasks)} tasks, {(time.monotonic() - started) / 60:.1f} min', flush=True)


def load(configs, hours, inputs, cache):
    frames = []
    for name in configs:
        for t in hours:
            for site in sorted(inputs[t]):
                with np.load(cache / name / f'{t:%Y%m%dT%H%M}_{site}.npz') as a:
                    n = len(a['observed'])
                    frames.append(pd.DataFrame({'config': name, 'time': t, 'site': site, 'bin': np.arange(n),
                                                **{k: a[k] for k in ('observed', 'eligible', 'pred', 'uerr', 'msfr',
                                                                     'ror')}}))
    d = pd.concat(frames, ignore_index=True)
    return d[d['eligible'] & np.isfinite(d['observed'])]


def accepted(d, uerr=0.6, share=None, ror=None):
    keep = np.isfinite(d['pred']) & (d['uerr'] <= uerr)
    if share is not None:
        keep &= d['msfr'] <= share
    if ror is not None:
        keep &= d['ror'] <= ror
    return keep


def curves(d, n_eligible):
    """Coverage and RMSE of every configuration at every uncertainty threshold."""
    rows = []
    for name, g in d.groupby('config'):
        err = g['pred'] - g['observed']
        for thr in UNCERTAINTY:
            ok = accepted(g, thr)
            rows.append(dict(config=name, uncertainty=thr, coverage_pct=100 * ok.sum() / n_eligible,
                             rmse_cm_s=np.sqrt(np.mean(err[ok] ** 2)) if ok.any() else np.nan))
    return pd.DataFrame(rows)


def rmse_at_coverage(curve, levels):
    """RMSE of each configuration interpolated to fixed coverage levels (equal-coverage comparison)."""
    rows = []
    for name, g in curve.groupby('config'):
        g = g.dropna().sort_values('coverage_pct')
        row = {'config': name}
        for level in levels:
            inside = g['coverage_pct'].min() <= level <= g['coverage_pct'].max()
            row[f'rmse_at_{level}pct'] = np.interp(level, g['coverage_pct'], g['rmse_cm_s']) if inside else np.nan
        rows.append(row)
    return pd.DataFrame(rows).set_index('config')


def paired_vs_production(d, configs, uerr=0.6):
    """Each configuration at a fixed uncertainty threshold, compared with production on the same observations."""
    key = ['time', 'site', 'bin']
    prod = d[d['config'] == PRODUCTION].set_index(key)
    prod_ok = accepted(prod, 0.6)
    rows = []
    n_eligible = len(prod)
    for name in configs:
        g = d[d['config'] == name].set_index(key).reindex(prod.index)
        ok = accepted(g, uerr)
        both = ok & prod_ok
        e = (g['pred'] - g['observed'])[both]
        r = (prod['pred'] - prod['observed'])[both]
        low, high, _ = paired_interval(e.to_numpy(), r.to_numpy(), e.index.get_level_values('time'))
        rows.append(dict(config=name, uncertainty=uerr, coverage_pct=100 * ok.sum() / n_eligible,
                         paired_bins=int(both.sum()), rmse_cm_s=np.sqrt(np.mean(e ** 2)),
                         production_rmse_cm_s=np.sqrt(np.mean(r ** 2)),
                         delta_rmse_cm_s=np.sqrt(np.mean(e ** 2)) - np.sqrt(np.mean(r ** 2)),
                         ci_low_24h=low, ci_high_24h=high))
    return pd.DataFrame(rows).set_index('config')


def selection_rules(d, name, n_eligible):
    """Does a site-share or ROR limit select better vectors than the uncertainty limit alone, at equal coverage?"""
    g = d[d['config'] == name]
    err2 = (g['pred'] - g['observed']) ** 2
    rows = []
    for rule, kwargs in [('uncertainty only', {}), ('share <= 0.75', dict(share=0.75)), ('share <= 0.9', dict(share=0.9)),
                         ('ROR <= 4', dict(ror=4)), ('ROR <= 8', dict(ror=8))]:
        for thr in UNCERTAINTY:
            ok = accepted(g, thr, **kwargs)
            rows.append(dict(rule=rule, uncertainty=thr, coverage_pct=100 * ok.sum() / n_eligible,
                             rmse_cm_s=np.sqrt(err2[ok].mean()) if ok.any() else np.nan))
    return pd.DataFrame(rows)


def check_production(d, study_dir):
    """The production configuration at Uerr/Verr <= 0.6 must reproduce maracoos_study.py's OI production predictions."""
    latest = study_dir / 'latest.json'
    if not latest.exists():
        return 'maracoos_study.py output not found; production check skipped'
    run = Path(json.loads(latest.read_text())['run'])
    study = pd.read_parquet(run / 'withheld_predictions.parquet')['OI production']
    prod = d[d['config'] == PRODUCTION].set_index(['time', 'site', 'bin'])
    ours = prod['pred'].where(accepted(prod, 0.6))
    common = ours.index.intersection(study.index)
    a, b = ours.loc[common], study.loc[common]
    same_mask = bool((a.isna() == b.isna()).all())
    diff = float(np.nanmax(np.abs(a - b))) if a.notna().any() else np.nan
    return (f'Production configuration vs maracoos_study.py OI production, {len(common):,} withheld bins: '
            f'same predicted bins {same_mask}, largest difference {diff:.2e} cm/s')


def figures(curve, configs, out):
    fig, axes = plt.subplots(1, 3, figsize=(17, 5), constrained_layout=True, sharey=True)
    families = ['length scale x variance ratio', 'window', 'radial error']
    for ax, family in zip(axes, families):
        names = [n for n, c in configs.items() if c['family'] == family or n == PRODUCTION]
        cmap = plt.get_cmap('tab20')
        for k, name in enumerate(names):
            g = curve[curve['config'] == name]
            style = dict(color='black', lw=2.4, zorder=5) if name == PRODUCTION else dict(color=cmap(k % 20), lw=1.4)
            ax.plot(g['coverage_pct'], g['rmse_cm_s'], label=name, **style)
            at = g[np.isclose(g['uncertainty'], 0.6)]
            ax.plot(at['coverage_pct'], at['rmse_cm_s'], 'o', ms=4, color=style['color'])
        ax.set_title(family)
        ax.set_xlabel('Withheld radials predicted (%)')
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, frameon=False, ncol=2)
    axes[0].set_ylabel('RMS error of predicted radials (cm/s)')
    fig.suptitle('Accuracy versus coverage as the Uerr/Verr limit varies from 0.3 to 1.0 (dots: 0.6)')
    fig.savefig(out / 'figures' / 'accuracy_coverage_sweep.png', dpi=130)
    plt.close(fig)


def markdown(df, digits=3):
    df = df.reset_index() if df.index.name else df
    head = '| ' + ' | '.join(map(str, df.columns)) + ' |\n|' + '---|' * len(df.columns) + '\n'
    fmt = lambda v: f'{v:.{digits}f}' if isinstance(v, (float, np.floating)) else str(v)
    return head + '\n'.join('| ' + ' | '.join(fmt(v) for v in row) + ' |' for row in df.itertuples(index=False))


def analyze(d, configs, hours, split, out, study_dir):
    (out / 'tables').mkdir(exist_ok=True)
    (out / 'figures').mkdir(exist_ok=True)
    n_eligible = len(d[d['config'] == PRODUCTION])
    check = check_production(d, study_dir)

    curve = curves(d, n_eligible)
    curve.to_csv(out / 'tables' / 'accuracy_coverage.csv', index=False)
    prod_cov = curve[(curve['config'] == PRODUCTION) & np.isclose(curve['uncertainty'], 0.6)]['coverage_pct'].item()
    levels = [30, 40, int(round(prod_cov)), 50]
    equal = rmse_at_coverage(curve, levels).loc[list(configs)]
    equal.to_csv(out / 'tables' / 'rmse_at_equal_coverage.csv')
    paired = pd.concat([paired_vs_production(d, configs, u) for u in (0.6, 0.7)])
    paired.to_csv(out / 'tables' / 'paired_vs_production.csv')
    figures(curve, configs, out)

    # Selection on the first part of the week, check on the rest
    first, second = d[d['time'] < split], d[d['time'] >= split]
    n1, n2 = len(first[first['config'] == PRODUCTION]), len(second[second['config'] == PRODUCTION])
    level = int(round(prod_cov))
    e1 = rmse_at_coverage(curves(first, n1), [level])[f'rmse_at_{level}pct']
    e2 = rmse_at_coverage(curves(second, n2), [level])[f'rmse_at_{level}pct']
    best = e1.idxmin()
    holdout = pd.DataFrame({'selection period RMSE': e1, 'check period RMSE': e2}).loc[list(configs)]
    holdout.to_csv(out / 'tables' / 'split_check.csv')

    rules = selection_rules(d, best, n_eligible)
    rules.to_csv(out / 'tables' / 'selection_rules.csv', index=False)
    rule_equal = rmse_at_coverage(rules.rename(columns={'rule': 'config'}), [30, 40, level])

    report = f"""# MARACOOS OI settings and QC threshold sweep

Week {hours[0]:%Y-%m-%d %H:%M} to {hours[-1]:%Y-%m-%d %H:%M}, withheld-site test every {int((hours[1] - hours[0]) / pd.Timedelta('1h'))} hours,
{n_eligible:,} withheld radial observations. mdlvar = {MDLVAR:g} (cm/s)^2; the OI solution depends on mdlvar and errvar
only through their ratio. Production is `{PRODUCTION}` (Sx = Sy = 15 km, errvar 66) with Uerr/Verr <= 0.6.

{check}

## Accuracy versus coverage

Every configuration is scored at Uerr/Verr limits from 0.3 to 1.0. A curve below another is more accurate at the same
coverage. Production at 0.6 predicts {prod_cov:.1f} % of the withheld radials.

![Accuracy versus coverage](figures/accuracy_coverage_sweep.png)

RMS error (cm/s) at equal coverage, interpolated along each curve:

{markdown(equal.round(3))}

## Paired with production

Each configuration at a fixed limit, against production (0.6) on the observations both predict; 24-hour block
bootstrap intervals (exploratory, seven blocks).

{markdown(paired.round(3))}

## Selection check

Configurations ranked by RMSE at {level} % coverage on {hours[0]:%d %b} to {split - pd.Timedelta('1h'):%d %b %H:%M} and checked on the rest
of the week. Best in the selection period: `{best}`.

{markdown(holdout.round(3))}

## Site-share and ROR limits

For `{best}`: RMS error at equal coverage when the uncertainty limit is combined with a site-share or ROR limit.
A lower value than "uncertainty only" means the extra rule removes worse vectors than the uncertainty limit does.

{markdown(rule_equal.round(3))}
"""
    (out / 'sweep_report.md').write_text(report)
    print(report)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--data', type=Path, default=Path('~/Documents/hfradar_parity/2026-10-06_week'))
    parser.add_argument('--start', default='2026-10-06T00:00')
    parser.add_argument('--end', default='2026-10-12T23:00')
    parser.add_argument('--loso-step', type=int, default=6, help='Hours between withheld-site tests (default 6)')
    parser.add_argument('--split', default='2026-10-10T00:00', help='Select on hours before this, check on the rest')
    parser.add_argument('--jobs', type=int, default=4)
    parser.add_argument('--output', type=Path, default=EXAMPLES / 'output/maracoos_sweep')
    parser.add_argument('--study-output', type=Path, default=EXAMPLES / 'output/maracoos_totals_study')
    args = parser.parse_args()
    warnings.simplefilter('ignore', category=FutureWarning)
    logging.disable(logging.WARNING)

    configs = configurations()
    all_hours = pd.date_range(args.start, args.end, freq='h')
    hours = all_hours[::args.loso_step]
    conf, driver, grid, land = setup(args.data)
    inputs = load_inputs(conf, all_hours)
    cache = args.output / 'solutions'
    compute(configs, inputs, hours, grid, land, conf, cache, args.jobs)
    d = load(configs, hours, inputs, cache)
    analyze(d, configs, hours, pd.Timestamp(args.split), args.output, args.study_output)


if __name__ == '__main__':
    main()
