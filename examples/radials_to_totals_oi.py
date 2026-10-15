#!/usr/bin/env python
"""
@author Mike Smith
@email michaesm@marine.rutgers.edu
@purpose Compute optimal interpolation (OI) total vectors following Nwankwo et al. (2026), "Review of the Texas
         shelf currents using high frequency radar surface currents data", Ocean Dynamics 76:57.

Steps (each one can be turned off):
    1. OI with radials from t-1, t and t+1 hours (--time-window 1), Sx = Sy = 10 km, signal variance 625 and
       error variance 100 (cm/s)^2
    2. QC of the totals: maximum speed, normalized uncertainty of u and v (0.7), and site contribution (no site
       contributing more than 75 % of a vector's radials, or Fang et al. (2015) ROR < 4)
    3. Daily means (--daily)
    4. Gap filling by nearest-neighbor averaging and optionally Gauss-Markov interpolation (Nwankwo et al. 2023)

Radial files are expected in one directory per site: <radials>/<SITE>/RDLm_<SITE>_YYYY_MM_DD_HHMM.ruv
Each hour uses whichever sites have radials. The output is one NetCDF file of gridded totals (time, lat, lon).

Example usage:
    python radials_to_totals_oi.py ~/Downloads/radials_clean --start 2026-10-06 --hours 168
    python radials_to_totals_oi.py ~/Downloads/radials_clean --start 2026-10-06 --hours 168 --daily --gap-fill
    python radials_to_totals_oi.py ~/Downloads/radials_clean --start 2026-10-06 --hours 24 --time-window 0 \\
        --max-uncertainty 0 --max-site-share 0          # plain OI, no QC
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

from hfradarpy.gapfill import fill_gaps
from hfradarpy.totals import combineRadials
from radials_to_totals_benchmark import find_radial_files, load_radial_sets
from synthetic_totals_benchmark import load_grid

logger = logging.getLogger(__name__)
log_level = 'INFO'
log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
logging.basicConfig(stream=sys.stdout, format=log_format, level=log_level)

# pandas emits a lot of FutureWarnings from inside combineRadials
warnings.simplefilter('ignore', category=FutureWarning)

# Variables written to the output, per hour and grid point
OUTPUT = ['VELU', 'VELV', 'UQAL', 'VQAL', 'CQAL', 'NRAD', 'NSIT', 'ROR', 'MSFR']
# QC flags: maximum speed, u and v uncertainty, site contribution, primary flag (NaN when a test is not run)
QC_FLAGS = {'Q303': 'maximum speed', 'Q306': 'u uncertainty', 'Q307': 'v uncertainty', 'SCON': 'site contribution',
            'PRIM': 'primary flag'}
COLUMNS = OUTPUT + list(QC_FLAGS)

_WORKER_GRID = None


def _init_worker(grid_file):
    global _WORKER_GRID
    warnings.simplefilter('ignore', category=FutureWarning)
    _WORKER_GRID = load_grid(grid_file)


def combine_hour(grid, ts, rdf, combine_kwargs, qc):
    """
    OI total vectors for one hour, with QC flags.

    :return: (array (points, len(COLUMNS)), warning string)
    """
    out = np.full((len(grid), len(COLUMNS)), np.nan)
    tot, warn = combineRadials(rdf, grid, 10000, 6000, ts, method='oi', engine='vectorized',
                               useLandMask=False, dropPoints=False, **combine_kwargs)
    if warn:
        return out, warn

    tot.initialize_qc()
    if qc['max_speed']:
        tot.qc_qartod_maximum_velocity(max_speed=qc['max_speed'], high_speed=qc['max_speed'])
    if qc['max_uncertainty']:
        tot.qc_qartod_u_uncertainty(max_Uerr=qc['max_uncertainty'])
        tot.qc_qartod_v_uncertainty(max_Verr=qc['max_uncertainty'])
    if qc['max_site_share'] or qc['max_ror']:
        tot.qc_site_contribution(max_ror=qc['max_ror'] or None, max_share=qc['max_site_share'] or None)
    if tot.metadata['QCTest']:
        tot.qc_qartod_primary_flag()
    else:
        tot.data['PRIM'] = np.where(tot.data['VELU'].notna(), 1, np.nan)

    return tot.data.reindex(columns=COLUMNS).to_numpy(dtype=float), ''


def _combine_hour_in_worker(args):
    return combine_hour(_WORKER_GRID, *args)


def to_gridded(values, times, grid_lon, grid_lat):
    """(time, point, variable) array on a regular grid -> xarray.Dataset with dims (time, lat, lon)."""
    lons = np.unique(np.round(grid_lon, 6))
    lats = np.unique(np.round(grid_lat, 6))
    ix = np.searchsorted(lons, np.round(grid_lon, 6))
    iy = np.searchsorted(lats, np.round(grid_lat, 6))
    ds = xr.Dataset(coords=dict(time=times, lat=lats, lon=lons))
    for k, name in enumerate(COLUMNS):
        field = np.full((len(times), lats.size, lons.size), np.nan)
        field[:, iy, ix] = values[:, :, k]
        ds[name] = (('time', 'lat', 'lon'), field)
        if name in QC_FLAGS:
            ds[name].attrs.update(long_name=f'QC flag: {QC_FLAGS[name]}', flag_values='1 2 3 4',
                                  flag_meanings='pass not_evaluated suspect fail')
    return ds


def main(radial_dir, grid_file, output, sites, pattern, start, hours, time_window, mdlvar, errvar, sx, sy,
         search_factor, min_sites, min_radials, qc, keep_flagged, daily, gap_fill, gauss_markov, workers):
    times = pd.date_range(start, periods=hours, freq='h')
    if not sites:
        sites = sorted(p.name for p in Path(radial_dir).iterdir() if p.is_dir() and any(p.glob(f'{pattern}_*.ruv')))

    grid = load_grid(grid_file)
    file_times = pd.date_range(times[0] - pd.Timedelta(hours=time_window),
                               times[-1] + pd.Timedelta(hours=time_window), freq='h')
    files = find_radial_files(radial_dir, sites, pattern, file_times)
    t0 = time.perf_counter()
    radial_sets = load_radial_sets(files, workers, times=times, window=time_window)
    logger.info(f'Loaded {int(files.notna().to_numpy().sum())} {pattern} files for {len(sites)} sites '
                f'in {time.perf_counter() - t0:.1f} s')

    combine_kwargs = dict(mdlvar=mdlvar, errvar=errvar, sx=sx, sy=sy, minContrSites=min_sites,
                          minContrRads=min_radials, oi_search_factor=search_factor)
    if time_window:
        combine_kwargs['tempthreshold'] = pd.Timedelta(hours=time_window)

    t0 = time.perf_counter()
    tasks = [(ts, rdf, combine_kwargs, qc) for ts, rdf in zip(times, radial_sets)]
    if workers > 1:
        with ProcessPoolExecutor(workers, initializer=_init_worker, initargs=(grid_file,)) as pool:
            results = list(pool.map(_combine_hour_in_worker, tasks))
    else:
        results = [combine_hour(grid, *task) for task in tasks]
    for ts, (_, warn) in zip(times, results):
        if warn:
            logger.warning(f'{ts}: {warn}')
    values = np.stack([r[0] for r in results])
    logger.info(f'Combined {hours} hours in {time.perf_counter() - t0:.1f} s')

    ds = to_gridded(values, times, grid.x.to_numpy(), grid.y.to_numpy())
    has_vector = np.isfinite(ds['VELU'])
    computed = int(has_vector.sum())
    flagged = int(((ds['PRIM'] == 4) & has_vector).sum())
    logger.info(f'{computed / hours:.0f} vectors per hour; {flagged / max(computed, 1):.1%} failed QC')
    for name, label in QC_FLAGS.items():
        if name != 'PRIM' and np.isfinite(ds[name]).any():
            failed = int(((ds[name] == 4) & has_vector).sum())
            logger.info(f'  {label}: {failed / max(computed, 1):.1%} failed')
    if not keep_flagged:
        bad = ds['PRIM'] == 4
        ds['VELU'] = ds['VELU'].where(~bad)
        ds['VELV'] = ds['VELV'].where(~bad)

    if daily:
        ds = ds[['VELU', 'VELV']].resample(time='1D').mean()
    if gap_fill:
        before = int(np.isfinite(ds['VELU']).sum())
        ds = fill_gaps(ds, gauss_markov=gauss_markov)
        flags = ds['fill_flag'].values
        logger.info(f'Gap filling: {before} observed values, {(flags == 1).sum()} filled by nearest-neighbor '
                    f'averaging, {(flags == 2).sum()} by Gauss-Markov interpolation')

    ds['VELU'].attrs.update(units='cm s-1', long_name='Eastward sea water velocity')
    ds['VELV'].attrs.update(units='cm s-1', long_name='Northward sea water velocity')
    ds.attrs.update(method='optimal interpolation (hfradarpy combineRadials, engine=vectorized)',
                    reference='Nwankwo et al. (2026), Ocean Dynamics 76:57', sites=' '.join(sites),
                    radial_type=pattern, time_window_hours=time_window, mdlvar=mdlvar, errvar=errvar, sx_km=sx,
                    sy_km=sy, oi_search_factor=search_factor, qc=str(qc), qc_flagged_removed=int(not keep_flagged))
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(output)
    logger.info(f'Saved {output}')
    return ds


if __name__ == '__main__':
    examples = Path(__file__).resolve().parent

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('radials', type=Path, help='Directory with one sub-directory of radial files per site')
    parser.add_argument('--sites', nargs='+', help='Site codes to use (default: every site directory with files)')
    parser.add_argument('--pattern', default='RDLm', choices=['RDLm', 'RDLi'], help='Radial file type')
    parser.add_argument('--grid', default=examples / 'data' / 'grid_files' / 'maracoos_grid_6km_extended.txt',
                        type=Path, help='Two column (lon lat) grid file')
    parser.add_argument('--output', default=examples / 'output' / 'totals' / 'oi_totals.nc', type=Path,
                        help='Output NetCDF file')
    parser.add_argument('--start', required=True, help='First hour (UTC), e.g. 2026-10-06T00:00')
    parser.add_argument('--hours', default=24, type=int, help='Number of hours (default 24)')
    parser.add_argument('--time-window', default=1, type=int,
                        help='Also use each site\'s radials this many hours before and after (default 1; 0 = off)')
    parser.add_argument('--mdlvar', default=625.0, type=float, help='Signal variance (cm/s)^2 (default 625)')
    parser.add_argument('--errvar', default=100.0, type=float, help='Error variance (cm/s)^2 (default 100)')
    parser.add_argument('--sx', default=10.0, type=float, help='East-west decorrelation length in km (default 10)')
    parser.add_argument('--sy', default=10.0, type=float, help='North-south decorrelation length in km (default 10)')
    parser.add_argument('--search-factor', default=2.0, type=float,
                        help='Use radials within this many decorrelation lengths (default 2)')
    parser.add_argument('--min-sites', default=2, type=int, help='Minimum contributing sites (paper also tried 3)')
    parser.add_argument('--min-radials', default=3, type=int, help='Minimum contributing radials')
    parser.add_argument('--max-speed', default=150.0, type=float, help='QC: maximum speed in cm/s (0 = off)')
    parser.add_argument('--max-uncertainty', default=0.7, type=float,
                        help='QC: maximum normalized uncertainty of u and v (default 0.7; 0 = off)')
    parser.add_argument('--max-site-share', default=0.75, type=float,
                        help='QC: maximum share of radials from one site (default 0.75; 0 = off)')
    parser.add_argument('--max-ror', default=0.0, type=float,
                        help='QC: maximum ratio of overlapping radials, Fang et al. 2015 (e.g. 4; default 0 = off)')
    parser.add_argument('--keep-flagged', action='store_true', help='Keep vectors that fail QC in VELU/VELV')
    parser.add_argument('--daily', action='store_true', help='Average to daily fields before gap filling')
    parser.add_argument('--gap-fill', action='store_true', help='Fill gaps by nearest-neighbor averaging')
    parser.add_argument('--gauss-markov', action='store_true', help='Then fill remaining gaps by Gauss-Markov')
    parser.add_argument('--workers', default=1, type=int, help='Processes for loading radials and combining hours')
    args = parser.parse_args()

    qc = dict(max_speed=args.max_speed, max_uncertainty=args.max_uncertainty, max_site_share=args.max_site_share,
              max_ror=args.max_ror)
    main(args.radials, args.grid, args.output, args.sites, args.pattern, args.start, args.hours, args.time_window,
         args.mdlvar, args.errvar, args.sx, args.sy, args.search_factor, args.min_sites, args.min_radials, qc,
         args.keep_flagged, args.daily, args.gap_fill or args.gauss_markov, args.gauss_markov, args.workers)
