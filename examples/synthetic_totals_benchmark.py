#!/usr/bin/env python
"""
@author Mike Smith
@email michaesm@marine.rutgers.edu
@purpose Benchmark total vector combination (unweighted least squares and optimal interpolation) on the
         MARACOOS 6 km extended grid using one week of synthetic, hourly radial data.

The synthetic radials are generated in memory from a known analytic surface current field so that the
recovered totals can be compared against the truth. Nothing is read from disk except the grid file.

The current field (cm/s) is the sum of:
    * a steady south-westward shelf flow
    * an M2 tidal current (12.42 h period) whose phase varies across the domain
    * a diurnal (24 h) wind-driven oscillation
    * an anticyclonic Gaussian eddy that drifts south-west through the domain

Each site sees every 5 degree bearing out to a time-varying maximum range, restricted to bins that fall over
the ocean grid. Random bins are dropped, sites occasionally go offline, and Gaussian noise is added to the
radial velocities.

Example usage:
    python synthetic_totals_benchmark.py                       # one day, uwls + oi, vectorized + loop engines
    python synthetic_totals_benchmark.py --hours 168           # one week
    python synthetic_totals_benchmark.py --engines vectorized  # skip the slow loop engine
    python synthetic_totals_benchmark.py --loop-hours 0        # run the loop engine for every hour (slow!)
    python synthetic_totals_benchmark.py --workers 4           # spread the vectorized hours over 4 processes
    python synthetic_totals_benchmark.py --methods uwls wls oi --hours 6
    python synthetic_totals_benchmark.py --oi-serial           # loop engine OI without joblib parallelism

The loop engine (the original per grid point implementation) only runs for the first --loop-hours hours by
default. Its time for the full period is extrapolated, and its output is checked against the vectorized engine.
"""
import argparse
import logging
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from collections import OrderedDict
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from pyproj import Geod
from scipy.spatial import cKDTree
from shapely.geometry import Point

from hfradarpy.radials import Radial
from hfradarpy.totals import combineRadials

logger = logging.getLogger(__name__)
log_level = 'INFO'
log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
logging.basicConfig(stream=sys.stdout, format=log_format, level=log_level)

# pandas emits a lot of FutureWarnings from inside combineRadials; they are not useful for a benchmark
warnings.simplefilter('ignore', category=FutureWarning)

GEOD = Geod(ellps='WGS84')

# Approximate locations of long-range (5 MHz) sites along the U.S. East Coast. Positions are illustrative
# only and are used solely to give the synthetic radials a realistic geometry.
SITES = OrderedDict([
    ('NANT', (41.2420, -69.9880)),
    ('MVCO', (41.3250, -70.5650)),
    ('BLCK', (41.1720, -71.5560)),
    ('AMAG', (40.9950, -72.1330)),
    ('MRCH', (40.7860, -72.7520)),
    ('HEMP', (40.5860, -73.5530)),
    ('HOOK', (40.4390, -73.9940)),
    ('LOVE', (39.7220, -74.1150)),
    ('BRIG', (39.4050, -74.3730)),
    ('WILD', (38.9430, -74.8650)),
    ('ASSA', (37.9800, -75.3630)),
    ('CEDR', (37.6600, -75.5990)),
    ('LISL', (36.6920, -75.9190)),
    ('DUCK', (36.1820, -75.7510)),
    ('HATY', (35.2560, -75.5180)),
    ('CORE', (34.7470, -76.4120)),
])

# Reference point used to convert lon/lat to km for the analytic current field
REF_LON, REF_LAT = -72.0, 39.0


def lonlat_to_km(lon, lat):
    """Equirectangular projection (km) about REF_LON/REF_LAT. Accurate enough for an analytic field."""
    x = (np.asarray(lon) - REF_LON) * 111.32 * np.cos(np.deg2rad(REF_LAT))
    y = (np.asarray(lat) - REF_LAT) * 110.57
    return x, y


def true_currents(lon, lat, hours):
    """
    Analytic surface current field.

    :param lon: longitudes (degrees)
    :param lat: latitudes (degrees)
    :param hours: hours since the start of the simulation
    :return: u, v (cm/s)
    """
    x, y = lonlat_to_km(lon, lat)

    # Steady south-westward shelf flow
    u = np.full(x.shape, -6.0)
    v = np.full(x.shape, -5.0)

    # M2 tidal current with a phase that propagates across the domain
    m2 = 2 * np.pi * hours / 12.42 - 0.004 * x - 0.002 * y
    u = u + 12.0 * np.cos(m2)
    v = v + 8.0 * np.sin(m2)

    # Diurnal wind-driven oscillation
    k1 = 2 * np.pi * hours / 24.0
    u = u + 5.0 * np.sin(k1)
    v = v + 3.0 * np.cos(k1)

    # Anticyclonic (clockwise) Gaussian eddy drifting south-west at ~10 km/day
    vmax, length = 50.0, 50.0  # cm/s, km
    xc = 150.0 - 8.0 * hours / 24.0
    yc = 20.0 - 6.0 * hours / 24.0
    psi_amp = vmax * length * np.exp(0.5)
    g = np.exp(-((x - xc) ** 2 + (y - yc) ** 2) / (2 * length ** 2))
    u = u + psi_amp * (y - yc) / length ** 2 * g
    v = v - psi_amp * (x - xc) / length ** 2 * g

    return u, v


def load_grid(grid_file):
    """Load a two column (lon lat) grid file into a GeoSeries."""
    lonlat = np.loadtxt(grid_file)
    return gpd.GeoSeries([Point(xy) for xy in lonlat], crs='EPSG:4326')


def site_geometry(lat0, lon0, ocean_tree, range_res=5.8, max_range=220.0, bearing_res=5.0, ocean_tol=5.0):
    """
    Pre-compute the static radial bin geometry for a site, keeping only the bins that fall over the ocean grid.

    :return: DataFrame with LOND, LATD, BEAR, HEAD, RNGE for each bin
    """
    ranges = np.arange(1, int(max_range / range_res) + 1) * range_res
    bearings = np.arange(0, 360, bearing_res)
    rr, bb = np.meshgrid(ranges, bearings)
    rr, bb = rr.ravel(), bb.ravel()

    lon, lat, back_az = GEOD.fwd(np.full(rr.size, lon0), np.full(rr.size, lat0), bb, rr * 1000)

    # Keep bins that are within ocean_tol km of an ocean grid point
    dist, _ = ocean_tree.query(np.column_stack(lonlat_to_km(lon, lat)))
    keep = dist <= ocean_tol

    return pd.DataFrame({
        'LOND': lon[keep],
        'LATD': lat[keep],
        'BEAR': bb[keep],
        'HEAD': np.mod(back_az[keep], 360),  # direction from the bin toward the radar (deg true)
        'RNGE': rr[keep],
    })


def make_radial(site, lat0, lon0, geometry, timestamp, hours, rng, range_res=5.8, noise=5.0, dropout=0.10):
    """
    Build an in-memory Radial object for one site and one hour.

    VELO follows the CODAR convention: positive toward the radar.
    """
    # Maximum range breathes with a 24 h cycle plus some randomness (mimics diurnal ionospheric noise)
    max_range = 180.0 + 30.0 * np.sin(2 * np.pi * hours / 24.0 + list(SITES).index(site)) + rng.normal(0, 10)
    df = geometry[geometry['RNGE'] <= max_range]
    df = df[rng.random(len(df)) > dropout].copy()

    u, v = true_currents(df['LOND'].to_numpy(), df['LATD'].to_numpy(), hours)
    head = np.deg2rad(df['HEAD'].to_numpy())
    velo = u * np.sin(head) + v * np.cos(head) + rng.normal(0, noise, len(df))

    df['VELO'] = velo
    df['VELU'] = velo * np.sin(head)
    df['VELV'] = velo * np.cos(head)
    df['ETMP'] = np.abs(rng.normal(noise, 1.0, len(df)))
    df = df.reset_index(drop=True)

    r = Radial.__new__(Radial)
    r.data = df
    r.time = timestamp
    r.is_wera = False
    r._iscorrupt = False
    r.file_name = f'RDLm_{site}_{timestamp:%Y_%m_%d_%H%M}.ruv'
    r.file_path = ''
    r.metadata = OrderedDict([
        ('Site', site),
        ('TimeStamp', f'{timestamp:%Y %m %d  %H %M %S}'),
        ('TimeZone', '"UTC" +0.000 0 "Atlantic/Reykjavik"'),
        ('TimeCoverage', '180.000 Minutes'),
        ('Origin', f'{lat0:.7f}  {lon0:.7f}'),
        ('GreatCircle', '"WGS84" 6378137.000  298.257223562997'),
        ('RangeResolutionKMeters', f'{range_res:.6f}'),
        ('AntennaBearing', '0.0 True'),
        ('PatternType', 'Measured'),
    ])
    return r


def make_radial_dataframe(geometries, timestamp, hours, rng, outage=0.05):
    """Create the rDF DataFrame (index: site code, column: Radial) expected by combineRadials."""
    radials = {}
    for site, (lat0, lon0) in SITES.items():
        if rng.random() < outage:
            continue  # site is offline this hour
        radials[site] = make_radial(site, lat0, lon0, geometries[site], timestamp, hours, rng)
    return pd.DataFrame({'Radial': pd.Series(radials, dtype=object)})


TOTAL_COLUMNS = ['VELU', 'VELV', 'VELO', 'HEAD', 'UQAL', 'VQAL', 'CQAL', 'GDOP', 'NRAD']


def combine_hour(grid, ts, rdf, method, engine, search_radius, grid_res, method_kwargs):
    """Combine one hour of radials. Returns the compute time, the output columns and any warning."""
    t0 = time.perf_counter()
    tot, warn = combineRadials(rdf, grid, search_radius, grid_res, ts,
                               minContrSites=2, minContrRads=3,
                               useLandMask=False, dropPoints=False, engine=engine,
                               **method_kwargs[method])
    duration = time.perf_counter() - t0
    return duration, (None if warn else tot.data[TOTAL_COLUMNS].to_numpy(dtype=float)), warn


# Each worker process loads the grid once instead of receiving it with every hour
_WORKER_GRID = None


def _init_worker(grid_file):
    global _WORKER_GRID
    _WORKER_GRID = load_grid(grid_file)


def _combine_hour_in_worker(args):
    return combine_hour(_WORKER_GRID, *args)


def run_combination(grid, grid_file, times, radial_sets, method, engine, search_radius, grid_res, method_kwargs,
                    workers=1):
    """
    Combine every hour of radials with one method/engine. Returns per-hour timings and gridded results.

    With workers > 1 the hours are spread over a pool of processes; 'wall' is then the elapsed time for all
    hours, which is what matters for throughput, while 'durations' are the per-hour times inside the workers.
    """
    hours = len(radial_sets)
    out = {'durations': np.full(hours, np.nan),
           'total': np.full((hours, len(grid), len(TOTAL_COLUMNS)), np.nan)}  # every combineRadials output column

    def record(i, ts, rdf, result):
        duration, total, warn = result
        out['durations'][i] = duration
        if warn:
            logger.warning(f'{ts}: {warn}')
        else:
            out['total'][i] = total
        logger.info(f'{method}/{engine} {ts:%Y-%m-%d %H:%M} | {rdf.index.nunique()} sites | '
                    f'{np.isfinite(out["total"][i, :, 0]).sum():5d} vectors | {duration:8.3f} s')

    t0 = time.perf_counter()
    if workers > 1:
        tasks = [(ts, rdf, method, engine, search_radius, grid_res, method_kwargs)
                 for ts, rdf in zip(times, radial_sets)]
        with ProcessPoolExecutor(workers, initializer=_init_worker, initargs=(grid_file,)) as pool:
            for i, result in enumerate(pool.map(_combine_hour_in_worker, tasks)):
                record(i, times[i], radial_sets[i], result)
    else:
        for i, (ts, rdf) in enumerate(zip(times, radial_sets)):
            record(i, ts, rdf, combine_hour(grid, ts, rdf, method, engine, search_radius, grid_res, method_kwargs))
    out['wall'] = time.perf_counter() - t0

    out['u'], out['v'], out['gdop'] = (out['total'][:, :, TOTAL_COLUMNS.index(c)] for c in ('VELU', 'VELV', 'GDOP'))
    return out


def compare_engines(res, ref, rtol=1e-6, atol=1e-6):
    """
    Compare every output column of two runs over the hours they have in common.

    The engines sum the radial contributions in a different order, so values agree to floating point
    round-off rather than bit for bit; rtol/atol decide what counts as the same.
    """
    n = len(res['durations'])
    a, b = res['total'], ref['total'][:n]
    rows = []
    for k, col in enumerate(TOTAL_COLUMNS):
        x, y = a[:, :, k], b[:, :, k]
        same_mask = np.array_equal(np.isnan(x), np.isnan(y))
        both = ~np.isnan(x) & ~np.isnan(y)
        diff = np.abs(x[both] - y[both])
        if col == 'HEAD':
            diff = np.minimum(diff, 360 - diff)  # 0 and 360 degrees are the same direction
        rel = diff / np.maximum(np.abs(y[both]), np.finfo(float).tiny)
        tol_ok = bool(np.all(diff <= atol + rtol * np.abs(y[both])))
        rows.append({'column': col, 'values_compared': both.sum(), 'same_nan_mask': same_mask,
                     'max_abs_diff': diff.max() if diff.size else 0.0,
                     'max_rel_diff': rel.max() if rel.size else 0.0,
                     'identical': same_mask and tol_ok})
    return pd.DataFrame(rows).set_index('column')


def score(method, res, grid_lon, grid_lat, max_gdop):
    """Compare a run against the analytic truth."""
    hours = len(res['durations'])
    u_true, v_true = true_currents(grid_lon[None, :], grid_lat[None, :], np.arange(hours)[:, None])
    valid = np.isfinite(res['u'])
    if method != 'oi':
        # Least squares solutions with poor site geometry blow up; score them the way they would be
        # used operationally, after the QARTOD GDOP test. (OI stores a normalized error in GDOP instead.)
        valid &= res['gdop'] <= max_gdop
    err = np.hypot(res['u'] - u_true, res['v'] - v_true)[valid]
    return {
        'scored/hour': valid.sum() / hours,
        'median_err_cm/s': np.median(err) if err.size else np.nan,
        'rmse_u_cm/s': np.sqrt(np.mean((res['u'][valid] - u_true[valid]) ** 2)) if valid.any() else np.nan,
        'rmse_v_cm/s': np.sqrt(np.mean((res['v'][valid] - v_true[valid]) ** 2)) if valid.any() else np.nan,
    }, u_true, v_true


def main(grid_file, save_dir, start, hours, methods, engines, search_radius, grid_res, mdlvar, errvar, sx, sy,
         loop_hours=2, workers=1, max_gdop=2.0, oi_parallel=True, seed=42, save=True):
    save_dir = Path(save_dir)
    rng = np.random.default_rng(seed)

    # Grid
    t0 = time.perf_counter()
    grid = load_grid(grid_file)
    grid_lon, grid_lat = grid.x.to_numpy(), grid.y.to_numpy()
    ocean_tree = cKDTree(np.column_stack(lonlat_to_km(grid_lon, grid_lat)))
    geometries = {site: site_geometry(lat0, lon0, ocean_tree) for site, (lat0, lon0) in SITES.items()}
    logger.info(f'Loaded {len(grid)} grid points and built geometry for {len(SITES)} sites '
                f'in {time.perf_counter() - t0:.2f} s')
    for site, geom in geometries.items():
        logger.debug(f'{site}: {len(geom)} ocean bins')

    times = pd.date_range(start, periods=hours, freq='h')

    # Pre-generate all of the synthetic radials so that the combination timing is not polluted by it
    t0 = time.perf_counter()
    radial_sets = [make_radial_dataframe(geometries, ts, i, rng) for i, ts in enumerate(times)]
    n_bins = np.mean([sum(len(r.data) for r in rdf['Radial']) for rdf in radial_sets])
    logger.info(f'Generated {hours} hours of synthetic radials ({n_bins:.0f} radial bins/hour on average) '
                f'in {time.perf_counter() - t0:.2f} s')

    method_kwargs = {
        'uwls': dict(method='uwls'),
        'wls': dict(method='wls'),
        'oi': dict(method='oi', mdlvar=mdlvar, errvar=errvar, sx=sx, sy=sy, oi_option='exponential',
                   runParallel=oi_parallel),
    }

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
            stats, u_true, v_true = score(method, res, grid_lon, grid_lat, max_gdop)

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
                **stats,
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
                    {
                        'u': (('time', 'point'), res['u'],
                              {'units': 'cm s-1', 'long_name': f'{method} eastward velocity'}),
                        'v': (('time', 'point'), res['v'],
                              {'units': 'cm s-1', 'long_name': f'{method} northward velocity'}),
                        'u_true': (('time', 'point'), u_true, {'units': 'cm s-1'}),
                        'v_true': (('time', 'point'), v_true, {'units': 'cm s-1'}),
                        'gdop': (('time', 'point'), res['gdop'],
                                 {'comment': 'GDOP for least squares; normalized OI total error for oi'}),
                        'compute_seconds': (('time',), res['durations']),
                    },
                    coords={'time': times[:n], 'lon': ('point', grid_lon), 'lat': ('point', grid_lat)},
                    attrs={'method': method, 'engine': engine, 'grid_file': str(grid_file),
                           'comment': 'Synthetic benchmark data'},
                )
                fname = save_dir / f'synthetic_totals_{method}_{engine}_{times[0]:%Y%m%dT%H}_{times[n - 1]:%Y%m%dT%H}.nc'
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
    parser.add_argument('--grid', default=examples / 'data' / 'grid_files' / 'maracoos_grid_6km_extended.txt',
                        type=Path, help='Two column (lon lat) grid file')
    parser.add_argument('--output', default=examples / 'output' / 'totals' / 'synthetic', type=Path,
                        help='Directory for the output NetCDF files')
    parser.add_argument('--start', default='2025-01-01T00:00', help='First timestamp (UTC)')
    parser.add_argument('--hours', default=24, type=int, help='Number of hourly timesteps (default: 1 day)')
    parser.add_argument('--methods', nargs='+', default=['uwls', 'oi'], choices=['uwls', 'oi', 'wls'],
                        help='Combination methods to benchmark')
    parser.add_argument('--engines', nargs='+', default=['vectorized', 'loop'], choices=['vectorized', 'loop'],
                        help='combineRadials engines to benchmark')
    parser.add_argument('--loop-hours', default=2, type=int,
                        help='Only run the slow loop engine for this many hours and extrapolate its total time '
                             '(0 = run every hour)')
    parser.add_argument('--workers', default=1, type=int,
                        help='Spread the hours over this many processes (vectorized engine only). '
                             'About 4 is the useful limit on a 10-core laptop')
    parser.add_argument('--search-radius', default=10.0, type=float,
                        help='Search radius for (u)wls in km (OI uses 2 * decorrelation length)')
    parser.add_argument('--grid-res', default=6.0, type=float, help='Grid resolution in km (metadata only)')
    parser.add_argument('--mdlvar', default=420, type=int, help='OI a priori model variance (cm/s)^2')
    parser.add_argument('--errvar', default=66, type=int, help='OI observation error variance (cm/s)^2')
    parser.add_argument('--sx', default=10.0, type=float, help='OI zonal decorrelation length scale (km)')
    parser.add_argument('--sy', default=10.0, type=float, help='OI meridional decorrelation length scale (km)')
    parser.add_argument('--max-gdop', default=2.0, type=float,
                        help='GDOP threshold applied to least squares totals when scoring against the truth')
    parser.add_argument('--oi-serial', action='store_true', help='Disable joblib parallelism for OI (loop engine)')
    parser.add_argument('--seed', default=42, type=int, help='Random seed for the synthetic data')
    parser.add_argument('--no-save', action='store_true', help='Do not write NetCDF output')
    args = parser.parse_args()

    main(args.grid, args.output, args.start, args.hours, args.methods, args.engines,
         search_radius=args.search_radius * 1000, grid_res=args.grid_res * 1000,
         mdlvar=args.mdlvar, errvar=args.errvar, sx=args.sx, sy=args.sy, loop_hours=args.loop_hours,
         workers=args.workers,
         max_gdop=args.max_gdop, oi_parallel=not args.oi_serial, seed=args.seed, save=not args.no_save)
