#!/usr/bin/env python
"""
@author Mike Smith
@email michaesm@marine.rutgers.edu
@purpose Benchmark total vector combination (unweighted least squares and optimal interpolation) on real CODAR
         radial files, comparing the original loop engine with the vectorized engine of combineRadials.

Radial files are expected in one directory per site:
    <radials>/<SITE>/RDLm_<SITE>_YYYY_MM_DD_HHMM.ruv

Each hour is combined from whichever sites have a radial file for that hour. A site that is missing an hour is
simply left out of that hour's combination; combineRadials still requires at least two contributing sites and
three radials at each grid point.

Example usage:
    python radials_to_totals_benchmark.py ~/Downloads/radials_clean --start 2026-10-06 --hours 168
    python radials_to_totals_benchmark.py ~/Downloads/radials_clean --sites NANT BLCK AMAG --hours 24
    python radials_to_totals_benchmark.py ~/Downloads/radials_clean --engines vectorized --workers 4
    python radials_to_totals_benchmark.py ~/Downloads/radials_clean --pattern RDLi   # ideal pattern radials
    python radials_to_totals_benchmark.py ~/Downloads/radials_clean --time-window 1   # radials from t-1, t, t+1

The loop engine (the original per grid point implementation) only runs for the first --loop-hours hours by
default. Its time for the full period is extrapolated, and its output is checked against the vectorized engine.
"""
import argparse
import logging
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from hfradarpy.common import timestamp_from_lluv_filename
from hfradarpy.radials import Radial
from synthetic_totals_benchmark import TOTAL_COLUMNS, compare_engines, load_grid, run_combination

logger = logging.getLogger(__name__)
log_level = 'INFO'
log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
logging.basicConfig(stream=sys.stdout, format=log_format, level=log_level)

# pandas emits a lot of FutureWarnings from inside combineRadials; they are not useful for a benchmark
warnings.simplefilter('ignore', category=FutureWarning)


def find_radial_files(radial_dir, sites, pattern, times):
    """
    Find the radial file of every site for every hour.

    :return: DataFrame (index: hour, columns: sites) of file paths, None where a site has no file for that hour
    """
    files = pd.DataFrame(None, index=times, columns=sites, dtype=object)
    for site in sites:
        for f in sorted((Path(radial_dir) / site).glob(f'{pattern}_{site}_*.ruv')):
            ts = pd.Timestamp(timestamp_from_lluv_filename(f.name))
            if ts in files.index:
                files.loc[ts, site] = f
    return files


def load_radial(fname):
    """Load a radial file. Returns None if it is corrupt or has no data."""
    try:
        r = Radial(fname)
    except Exception:
        logger.exception(f'Could not load {fname}')
        return None
    if r._iscorrupt or r.data.empty:
        logger.warning(f'Skipping empty or corrupt radial file {fname}')
        return None
    return r


def load_radial_sets(files, workers=1, times=None, window=0):
    """
    Load every radial file and build one rDF per hour (index: site code, column: Radial) for combineRadials.
    Sites without a usable file for an hour are left out of that hour.

    With window > 0, the rDF of hour t holds the radials of every site from t - window to t + window hours (one row
    per radial, so a site code repeats), and files must also cover those neighboring hours. times are the hours to
    build an rDF for (default: every hour in files).
    """
    paths = [f for f in files.to_numpy().ravel() if f is not None and not pd.isna(f)]
    if workers > 1:
        with ProcessPoolExecutor(workers) as pool:
            loaded = dict(zip(paths, pool.map(load_radial, paths, chunksize=16)))
    else:
        loaded = {f: load_radial(f) for f in paths}

    radial_sets = []
    for ts in (files.index if times is None else times):
        rows = files.loc[ts - pd.Timedelta(hours=window):ts + pd.Timedelta(hours=window)]
        sites, radials = [], []
        for _, row in rows.iterrows():
            for site, f in row.items():
                if f is not None and not pd.isna(f) and loaded.get(f) is not None:
                    sites.append(site)
                    radials.append(loaded[f])
        radial_sets.append(pd.DataFrame({'Radial': pd.Series(radials, index=sites, dtype=object)}))
    return radial_sets


def main(radial_dir, grid_file, save_dir, sites, pattern, start, hours, methods, engines, search_radius, grid_res,
         mdlvar, errvar, sx, sy, loop_hours=2, workers=1, oi_parallel=True, save=True, time_window=0):
    save_dir = Path(save_dir)
    times = pd.date_range(start, periods=hours, freq='h')
    if time_window and 'loop' in engines:
        logger.warning('The loop engine cannot combine several radials per site; running the vectorized engine only')
        engines = [e for e in engines if e != 'loop'] or ['vectorized']
    # Radials from time_window hours before the first and after the last hour are needed too
    file_times = pd.date_range(times[0] - pd.Timedelta(hours=time_window), times[-1] + pd.Timedelta(hours=time_window),
                               freq='h')

    if not sites:
        sites = sorted(p.name for p in Path(radial_dir).iterdir() if p.is_dir() and any(p.glob(f'{pattern}_*.ruv')))

    grid = load_grid(grid_file)
    grid_lon, grid_lat = grid.x.to_numpy(), grid.y.to_numpy()

    files = find_radial_files(radial_dir, sites, pattern, file_times)
    available = files.loc[times].notna()
    logger.info(f'Found {int(available.to_numpy().sum())} {pattern} files for {len(sites)} sites over {hours} hours')
    for site in sites:
        logger.info(f'  {site}: {int(available[site].sum())}/{hours} hours')

    t0 = time.perf_counter()
    radial_sets = load_radial_sets(files, workers, times=times, window=time_window)
    load_seconds = time.perf_counter() - t0
    n_sites = np.array([rdf.index.nunique() for rdf in radial_sets])
    n_bins = np.array([sum(len(r.data) for r in rdf['Radial']) for rdf in radial_sets])
    logger.info(f'Loaded radials in {load_seconds:.1f} s: {n_sites.mean():.1f} sites and {n_bins.mean():.0f} '
                f'radial bins per hour on average ({n_sites.min()}-{n_sites.max()} sites)')

    method_kwargs = {
        'uwls': dict(method='uwls'),
        'wls': dict(method='wls'),
        'oi': dict(method='oi', mdlvar=mdlvar, errvar=errvar, sx=sx, sy=sy, oi_option='exponential',
                   runParallel=oi_parallel),
    }
    if time_window:
        for kwargs in method_kwargs.values():
            kwargs['tempthreshold'] = pd.Timedelta(hours=time_window)

    # The loop engine is slow, so by default it only runs the first few hours and its total is extrapolated
    engine_hours = {e: hours for e in engines}
    if 'loop' in engines and loop_hours > 0:
        engine_hours['loop'] = min(hours, loop_hours)

    summary = []
    comparisons = {}
    for method in methods:
        results = {}
        for engine in sorted(engines, key=lambda e: e != 'vectorized'):  # vectorized first
            n = engine_hours[engine]
            # Hours are only spread over processes for the vectorized engine; the loop engine's OI already
            # runs in parallel with joblib, and a pool around it would oversubscribe the CPU.
            nproc = workers if engine == 'vectorized' else 1
            logger.info(f'===== {method.upper()} / {engine} ({n} hours, {nproc} process(es)) =====')
            res = run_combination(grid, grid_file, times[:n], radial_sets[:n], method, engine, search_radius,
                                  grid_res, method_kwargs, workers=nproc)
            results[engine] = res

            row = {
                'method': method,
                'engine': engine,
                'hours_run': n,
                'processes': nproc,
                'mean_s/hour': np.nanmean(res['durations']),
                'max_s/hour': np.nanmax(res['durations']),
                'wall_s/hour': res['wall'] / n,
                f'est_{hours}h_min': res['wall'] / n * hours / 60,
                'vectors/hour': np.isfinite(res['u']).sum() / n,
            }
            if engine != 'vectorized' and 'vectorized' in results:
                # Agreement with the vectorized engine, column by column, over the hours both engines ran
                comparison = compare_engines(res, results['vectorized'])
                comparisons[method] = comparison
                row['outputs_identical'] = bool(comparison['identical'].all())
            summary.append(row)

            if save:
                save_dir.mkdir(parents=True, exist_ok=True)
                ds = xr.Dataset(
                    {col: (('time', 'point'), res['total'][:, :, k]) for k, col in enumerate(TOTAL_COLUMNS)},
                    coords={'time': times[:n], 'lon': ('point', grid_lon), 'lat': ('point', grid_lat)},
                    attrs={'method': method, 'engine': engine, 'pattern': pattern, 'sites': ' '.join(sites),
                           'grid_file': str(grid_file)},
                )
                ds['compute_seconds'] = (('time',), res['durations'])
                ds['n_sites'] = (('time',), n_sites[:n])
                ds['VELU'].attrs['units'] = ds['VELV'].attrs['units'] = ds['VELO'].attrs['units'] = 'cm s-1'
                fname = save_dir / f'totals_{method}_{engine}_{times[0]:%Y%m%dT%H}_{times[n - 1]:%Y%m%dT%H}.nc'
                ds.to_netcdf(fname)
                logger.info(f'Saved {fname}')

    summary = pd.DataFrame(summary).set_index(['method', 'engine'])
    if {'loop', 'vectorized'} <= set(engines):
        wall = summary['wall_s/hour']
        summary['speedup'] = np.nan
        for method in methods:
            summary.loc[(method, 'vectorized'), 'speedup'] = wall[(method, 'loop')] / wall[(method, 'vectorized')]
    with pd.option_context('display.float_format', '{:.3f}'.format, 'display.width', 250,
                           'display.max_columns', 50):
        print(f'\nRadial loading: {load_seconds:.1f} s for {int(available.to_numpy().sum())} files')
        print('\n' + summary.to_string())
    with pd.option_context('display.float_format', '{:.2e}'.format, 'display.width', 250):
        for method, comparison in comparisons.items():
            print(f'\n{method}: loop vs vectorized engine, first {engine_hours["loop"]} hour(s), all output columns')
            print(comparison.to_string())
        if comparisons and not all(c['identical'].all() for c in comparisons.values()):
            logger.error('The loop and vectorized engines produced different outputs!')
    return summary


if __name__ == '__main__':
    examples = Path(__file__).resolve().parent

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('radials', type=Path, help='Directory with one sub-directory of radial files per site')
    parser.add_argument('--sites', nargs='+', help='Site codes to use (default: every site directory with files)')
    parser.add_argument('--pattern', default='RDLm', choices=['RDLm', 'RDLi'],
                        help='Radial file type: measured (RDLm) or ideal (RDLi) antenna pattern')
    parser.add_argument('--grid', default=examples / 'data' / 'grid_files' / 'maracoos_grid_6km_extended.txt',
                        type=Path, help='Two column (lon lat) grid file')
    parser.add_argument('--output', default=examples / 'output' / 'totals' / 'real', type=Path,
                        help='Directory for the output NetCDF files')
    parser.add_argument('--start', required=True, help='First timestamp (UTC), e.g. 2026-10-06T00:00')
    parser.add_argument('--hours', default=24, type=int, help='Number of hourly timesteps (default: 1 day)')
    parser.add_argument('--methods', nargs='+', default=['uwls', 'oi'], choices=['uwls', 'oi', 'wls'],
                        help='Combination methods to benchmark')
    parser.add_argument('--engines', nargs='+', default=['vectorized', 'loop'], choices=['vectorized', 'loop'],
                        help='combineRadials engines to benchmark')
    parser.add_argument('--loop-hours', default=2, type=int,
                        help='Only run the slow loop engine for this many hours and extrapolate its total time '
                             '(0 = run every hour)')
    parser.add_argument('--time-window', default=0, type=int,
                        help='Also combine each site\'s radials from this many hours before and after each hour '
                             '(Nwankwo et al. 2026 used 1). Vectorized engine only')
    parser.add_argument('--workers', default=1, type=int,
                        help='Spread radial loading and the vectorized hours over this many processes')
    parser.add_argument('--search-radius', default=10.0, type=float,
                        help='Search radius for (u)wls in km (OI uses 2 * decorrelation length)')
    parser.add_argument('--grid-res', default=6.0, type=float, help='Grid resolution in km (metadata only)')
    parser.add_argument('--mdlvar', default=420, type=int, help='OI a priori model variance (cm/s)^2')
    parser.add_argument('--errvar', default=66, type=int, help='OI observation error variance (cm/s)^2')
    parser.add_argument('--sx', default=10.0, type=float, help='OI zonal decorrelation length scale (km)')
    parser.add_argument('--sy', default=10.0, type=float, help='OI meridional decorrelation length scale (km)')
    parser.add_argument('--oi-serial', action='store_true', help='Disable joblib parallelism for OI (loop engine)')
    parser.add_argument('--no-save', action='store_true', help='Do not write NetCDF output')
    args = parser.parse_args()

    main(args.radials, args.grid, args.output, args.sites, args.pattern, args.start, args.hours, args.methods,
         args.engines, search_radius=args.search_radius * 1000, grid_res=args.grid_res * 1000,
         mdlvar=args.mdlvar, errvar=args.errvar, sx=args.sx, sy=args.sy, loop_hours=args.loop_hours,
         workers=args.workers, oi_parallel=not args.oi_serial, save=not args.no_save, time_window=args.time_window)
