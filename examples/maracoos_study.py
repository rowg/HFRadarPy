"""Reproducible, offline MARACOOS OI comparison; no production files are modified.

Run from any directory with the hfradarpy environment::

    python examples/maracoos_study.py --data ~/Documents/hfradar_parity/2026-10-06_week

The dated input archive must contain prod/totals_configs_export.json and overrides.json.
Each expensive task is checkpointed atomically under an input/code fingerprint. The
study uses existing UWLS as a reference and does not change any combination algorithm.
"""
import argparse
import copy
import hashlib
import inspect
import json
import logging
import os
from pathlib import Path
import time
import warnings

import numpy as np
import pandas as pd
import xarray as xr
from joblib import Parallel, delayed, parallel_config
from scipy.spatial import cKDTree

import realtime_totals as rt
from hfradarpy.gapfill import fill_gaps
from hfradarpy.totals import combineRadials

EXAMPLES = Path(__file__).resolve().parent
BASELINE = 'OI production'
NW = dict(mdlvar=625., errvar=100., sx=10., sy=10., weighting_function='exponential')
VARIABLES = ['VELU', 'VELV', 'PRIM', 'UQAL', 'VQAL', 'MSFR', 'NRAD']
MAX_DISTANCE_KM = 4.3
SEED = 20261014


def configurations(prod):
    return {
        'uwls_prod': dict(method='uwls', window=0, outputs={'UWLS production': (None, None)}),
        'oi_prod': dict(method='oi', oi=prod, window=0, outputs={BASELINE: (None, None)}),
        'oi_prod_window': dict(method='oi', oi=prod, window=1,
                               outputs={'OI production + window': (None, None)}),
        'oi_nw': dict(method='oi', oi=NW, window=0,
                      outputs={'OI Nwankwo parameters': (None, None)}),
        'oi_nw_window': dict(method='oi', oi=NW, window=1, outputs={
            'OI Nwankwo parameters + window': (None, None),
            'OI Nwankwo + uncertainty QC': (0.7, None),
            'OI Nwankwo + site QC': (None, 0.75),
            'Nwankwo (full)': (0.7, 0.75),
        }),
    }


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def atomic_json(path, value):
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, default=str) + '\n')
    tmp.replace(path)


def verify_environment(run):
    """Prevent reuse across dependency changes even when the numerical source matches."""
    current = {m: __import__(m).__version__ for m in ('numpy', 'pandas', 'scipy', 'xarray', 'geopandas')}
    path = run / 'environment.json'
    if path.exists() and json.loads(path.read_text()) != current:
        raise RuntimeError('Cached dependency versions differ; use a new --output directory to recompute')
    if not path.exists():
        atomic_json(path, current)


def xyz(lon, lat):
    lon, lat = np.deg2rad(lon), np.deg2rad(lat)
    return 6371. * np.column_stack((np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)))


def setup(data, start, end, output, step):
    data = Path(data).expanduser().resolve()
    doc = rt.load_config('realtime_06km_maracoos_v03', data / 'prod/totals_configs_export.json')
    overrides = json.loads((data / 'overrides.json').read_text())
    driver = copy.deepcopy(rt.DRIVER_DEFAULTS)
    rt.deep_update(driver, overrides.pop('driver', {}))
    rt.deep_update(doc, overrides)
    conf = doc['conf']
    rconf, tconf = conf['Radials'], conf['Totals']
    grid = rt.load_grid(tconf['GridFile'])
    sites, jobs = rt.pattern_jobs(rconf)
    types = dict(jobs)['BestChoice']
    prod = {k: tconf['OI'][k] for k in ('mdlvar', 'errvar', 'sx', 'sy')}
    prod['weighting_function'] = tconf['OI'].get('weighting_function') or 'exponential'
    combos = configurations(prod)
    hours = pd.date_range(start, end, freq='h')
    if not len(hours):
        raise ValueError('The study requires at least one hour')
    filenames = {t: rt.radial_filenames(rconf['BaseDir'], sites, types, t,
                                       rconf.get('MonthFlag', False), rconf.get('TypeFlag', False))
                 for t in pd.date_range(hours[0] - pd.Timedelta('1h'), hours[-1] + pd.Timedelta('1h'), freq='h')}
    source_paths = [EXAMPLES / 'realtime_totals.py']
    source_paths += sorted((EXAMPLES.parent / 'hfradarpy').glob('*.py'))
    source_hashes = {str(p.relative_to(EXAMPLES.parent)): digest(p) for p in source_paths}
    # Analysis/report edits do not invalidate the expensive numerical task cache.
    compute_code = '\n'.join(inspect.getsource(f) for f in
                             (configurations, xyz, load_inputs, combine, quality_control, compute_task, setup))
    settings_files = [data / 'prod/totals_configs_export.json', data / 'overrides.json', Path(tconf['GridFile'])]
    settings_files += [Path(p) for p in (rconf.get('MaskDir'), conf['QC'].get('QC305_MaskFile')) if p]
    inventory = {str(p): digest(p) for paths in filenames.values() for p in paths if p.exists()}
    manifest = dict(start=str(hours[0]), end=str(hours[-1]), loso_step_hours=step, configurations=combos,
                    max_distance_km=MAX_DISTANCE_KM, sources=source_hashes,
                    study_compute_sha256=hashlib.sha256(compute_code.encode()).hexdigest(),
                    settings={str(p): digest(p) for p in settings_files}, radial_sha256=inventory,
                    seed=SEED, conf=conf, driver=driver,
                    notes=['No independent in situ data supplied.',
                           'Only the supplied week is evaluated; no parameter tuning is performed.',
                           'Centered windows require one additional hour of input availability.'])
    signature = hashlib.sha256(json.dumps(manifest, sort_keys=True, default=str).encode()).hexdigest()
    root = Path(output).resolve()
    run = root / signature[:16]
    for sub in ('totals', 'loso', 'tables', 'figures'):
        (run / sub).mkdir(parents=True, exist_ok=True)
    manifest['fingerprint'] = signature
    manifest['versions'] = {m: __import__(m).__version__ for m in ('numpy', 'pandas', 'scipy', 'xarray', 'geopandas')}
    atomic_json(run / 'manifest.json', manifest)
    atomic_json(root / 'latest.json', {'run': str(run), 'fingerprint': signature})
    ctx = dict(conf=conf, driver=driver, grid=grid, combos=combos, hours=hours, filenames=filenames,
               polygon=rt.load_polygon(conf['QC'].get('QC305_MaskFile')),
               radial_mask=rt.load_polygon(rconf.get('MaskDir')), run=run)
    print(f'Run {signature[:16]}: {len(hours)} hours, {len(grid)} grid points, {len(inventory)} input files', flush=True)
    return ctx


def load_inputs(ctx):
    inputs = {}
    for t, files in ctx['filenames'].items():
        radials, _, _ = rt.load_radials(files)
        rt.clean_radials(radials, remove_flagged=ctx['conf']['QC'].get('Process', False), mask=ctx['radial_mask'])
        for r in radials.values():
            r.data = r.data[['LOND', 'LATD', 'VELO', 'HEAD']].copy()
            for attr in ('_tables', 'diagnostics_radial', 'diagnostics_hardware', 'range_information', 'crad_data'):
                r.__dict__.pop(attr, None)
        inputs[t] = radials
    counts = pd.DataFrame([dict(time=t, site=site, bins=len(r.data))
                           for t, rs in inputs.items() for site, r in rs.items()])
    counts.to_csv(ctx['run'] / 'tables/input_counts.csv', index=False)
    print(f'Loaded {len(counts)} radial files', flush=True)
    return inputs


def combine(t, radials, combo, grid, conf, driver, exclude=None):
    rows, sites = [], []
    for k in range(-combo['window'], combo['window'] + 1):
        for site, radial in radials.get(t + pd.Timedelta(hours=k), {}).items():
            if site != exclude:
                rows.append(radial)
                sites.append(site)
    rdf = pd.DataFrame({'Radial': pd.Series(rows, index=sites, dtype=object)})
    tc = conf['Totals']
    if rdf.index.nunique() < tc.get('MinNumSites', 2) or grid.empty:
        return None
    kwargs = dict(gRes=6000, tStp=t.to_pydatetime(), minContrSites=tc.get('MinNumSites', 2),
                  minContrRads=tc.get('MinNumRads', 3), useLandMask=False, dropPoints=False, engine='vectorized')
    if combo['window']:
        kwargs['tempthreshold'] = pd.Timedelta(hours=combo['window'])
    if combo['method'] == 'uwls':
        total, _ = combineRadials(rdf, grid, sRad=tc['UWLS']['spatthresh'] * 1000, method='uwls',
                                 distance=driver['uwls_search_distance'], **kwargs)
    else:
        oi = combo['oi']
        total, _ = combineRadials(rdf, grid, sRad=2 * max(oi['sx'], oi['sy']) * 1000, method='oi',
                                 mdlvar=float(oi['mdlvar']), errvar=float(oi['errvar']), sx=float(oi['sx']),
                                 sy=float(oi['sy']), oi_option=oi['weighting_function'], **kwargs)
    return total if total.data.VELU.notna().any() else None


def quality_control(total, combo, variant, conf, driver, polygon):
    c = copy.deepcopy(conf)
    uncertainty, share = variant
    if uncertainty is not None:
        for test in ('QC306_OI_Uerr', 'QC307_OI_Verr'):
            c['QC']['GDOP_OI'][test]['MaxError'] = uncertainty
    rt.quality_control(total, c, combo['method'], polygon, driver)
    if share is not None:
        total.qc_site_contribution(max_share=share)
        total.qc_qartod_primary_flag(recalculate=True)
    return total


def compute_task(t, key, combo, radials, grid, conf, driver, polygon, run, site=None):
    warnings.simplefilter('ignore', category=FutureWarning)
    logging.disable(logging.WARNING)
    stamp = t.strftime('%Y%m%dT%H%M')
    path = run / ('loso' if site else 'totals') / f'{key}_{stamp}{"_" + site if site else ""}.npz'
    if path.exists():
        return str(path)
    data = {}
    if site:
        target = radials[t][site].data
        dist, point = cKDTree(xyz(grid.x.to_numpy(), grid.y.to_numpy())).query(
            xyz(target.LOND.to_numpy(), target.LATD.to_numpy()))
        eligible = dist <= MAX_DISTANCE_KM
        selected, inverse = np.unique(point[eligible], return_inverse=True)
        solve_grid = grid.iloc[selected].reset_index(drop=True)
        data.update(observed=target.VELO.to_numpy(float), eligible=eligible, distance_km=dist,
                    lon=target.LOND.to_numpy(float), lat=target.LATD.to_numpy(float))
    else:
        solve_grid = grid
    total = combine(t, radials, combo, solve_grid, conf, driver, exclude=site)
    for name, variant in combo['outputs'].items():
        tq = quality_control(copy.deepcopy(total), combo, variant, conf, driver, polygon) if total is not None else None
        if site:
            pred = np.full(len(target), np.nan)
            if tq is not None:
                d = tq.data
                u, v = d.VELU.to_numpy(float)[inverse], d.VELV.to_numpy(float)[inverse]
                ok = (d.PRIM.to_numpy()[inverse] == 1) & np.isfinite(u) & np.isfinite(v)
                head = np.deg2rad(target.HEAD.to_numpy(float)[eligible])
                pred[eligible] = np.where(ok, u * np.sin(head) + v * np.cos(head), np.nan)
            data[name] = pred
        else:
            for variable in VARIABLES:
                data[f'{name}::{variable}'] = (tq.data[variable].to_numpy(float) if tq is not None and variable in tq.data
                                               else np.full(len(grid), np.nan))
    tmp = path.with_suffix('.npz.tmp')
    with tmp.open('wb') as f:
        np.savez_compressed(f, **data)
    tmp.replace(path)
    return str(path)


def compute(ctx, inputs, jobs, step):
    started = time.monotonic()
    order = ['oi_prod_window', 'oi_nw_window', 'oi_prod', 'oi_nw', 'uwls_prod']
    # One preflight case per algorithm/QC family, before scheduling the full week.
    t = ctx['hours'][0]
    for key in order:
        compute_task(t, key, ctx['combos'][key], {h: inputs[h] for h in inputs if abs(h-t) <= pd.Timedelta('1h')},
                     ctx['grid'], ctx['conf'], ctx['driver'], ctx['polygon'], ctx['run'])
    print('Preflight passed for all methods and QC variants', flush=True)
    for stage in ('totals', 'loso'):
        tasks = [(key, t, site) for key in order
                 for t in (ctx['hours'] if stage == 'totals' else ctx['hours'][::step])
                 for site in ([None] if stage == 'totals' else sorted(inputs[t]))]
        with parallel_config(backend='loky', inner_max_num_threads=1):
            completed = Parallel(n_jobs=jobs, return_as='generator_unordered', pre_dispatch=jobs)(
                delayed(compute_task)(t, key, ctx['combos'][key],
                    {t + pd.Timedelta(hours=k): inputs.get(t + pd.Timedelta(hours=k), {})
                     for k in range(-ctx['combos'][key]['window'], ctx['combos'][key]['window'] + 1)},
                    ctx['grid'], ctx['conf'], ctx['driver'], ctx['polygon'], ctx['run'], site)
                for key, t, site in tasks)
            for i, _ in enumerate(completed, 1):
                if i == 1 or i % 25 == 0 or i == len(tasks):
                    print(f'{stage}: {i}/{len(tasks)} tasks; {(time.monotonic()-started)/60:.1f} min elapsed', flush=True)
    atomic_json(ctx['run'] / 'compute_complete.json', dict(elapsed_seconds=time.monotonic()-started))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=Path('~/Documents/hfradar_parity/2026-10-06_week'))
    parser.add_argument('--start', default='2026-10-06T00:00')
    parser.add_argument('--end', default='2026-10-12T23:00')
    parser.add_argument('--output', type=Path, default=EXAMPLES / 'output/maracoos_totals_study')
    parser.add_argument('--jobs', type=int, default=4)
    parser.add_argument('--loso-step', type=int, default=6)
    parser.add_argument('--stage', choices=('compute', 'analysis', 'all'), default='all')
    parser.add_argument('--insitu-file', type=Path, help='Optional independent CSV: time, lon, lat, u, v')
    parser.add_argument('--insitu-units', choices=('cm/s', 'm/s'), default='cm/s')
    args = parser.parse_args()
    if args.loso_step < 1 or args.jobs < 1:
        parser.error('--jobs and --loso-step must be positive')
    warnings.simplefilter('ignore', category=FutureWarning)
    logging.disable(logging.WARNING)
    ctx = setup(args.data, args.start, args.end, args.output, args.loso_step)
    verify_environment(ctx['run'])
    ctx['insitu_file'], ctx['insitu_units'] = args.insitu_file, args.insitu_units
    if args.stage in ('compute', 'all'):
        inputs = load_inputs(ctx)
        compute(ctx, inputs, args.jobs, args.loso_step)
    if args.stage in ('analysis', 'all'):
        from maracoos_study_analysis import analyze
        analyze(ctx)


if __name__ == '__main__':
    main()
