#!/usr/bin/env python
"""
@author Mike Smith
@email michaesm@marine.rutgers.edu
@purpose Generic, configuration driven realtime total vector processing with hfradarpy. Python replacement for the
         MATLAB HFRProgs/Rutgers totals toolbox drivers (realtime_maracoos_06km_combine.m,
         realtime_yucatan_long_range.m, create_hourly_totals.m and codar_driver_totals.m) and for the .mat to NetCDF
         conversion that followed them (realtime_qc_totals_to_netcdf.sh).

The processing configuration is the same 'totals_configs' document the MATLAB drivers use (and hfradmin edits), read
from a JSON file with --config-file (a mongoexport or hfradmin export; pick a document with --config when the file holds
several), or from MongoDB by name with --mongo --config NAME. An --overrides JSON file is merged on top of it for
settings that differ on this machine (radial, grid, mask and output paths) and for the driver options below.

MongoDB is optional: without --mongo the driver never connects to it and pymongo does not need to be installed.

For every hour of the processing window (newest first) and every pattern type enabled in Radials.ProcessByPattern:
    1. Build the radial file name of each site (filenames_standard_filesystem.m). RDL types are radial (.ruv) files and
       ELT types are elliptical (.euv) files from bistatic sites; both are combined the same way.
    2. Load the files that exist and are not empty (loadRDLFile_all.m)
    3. When QC.Process is true, drop radial vectors flagged as bad (PRIM == 4) by the radial QC (cleanQCedRadials.m)
    4. When Radials.MaskDir is set, drop radial vectors inside that (land) polygon (maskRadials.m)
    5. Combine them into totals with the methods enabled in Totals.UWLS.Process and Totals.OI.Process
       (makeTotals.m / tuvLS.m and makeTotalsOI.m / tuvOI.m)
    6. Quality control the totals as codar_driver_totals.m does:
         UWLS: QC303 maximum speed, QC305 valid location, QC302 GDOP (HFRProgs definition)
         OI:   QC303 maximum speed, QC305 valid location, QC306 Uerr, QC307 Verr,
               and QC308 total error / QC309 UV covariance when they are defined in QC.GDOP_OI
         PRIM: worst flag of the tests above
    7. Save each method as a CF/NCEI NetCDF file in the production format (u, v, u_err, ..., qc303_maxspeed, ...)
    8. With --mongo, save the hour's status in the 'totals_generated' field of the configuration document (hfradmin)

Driver options (the "driver" section of the overrides file). The two MATLAB drivers hard code these differently, so
they are not part of the MongoDB configuration:
    qc305_polygon:         'land' (default, MARACOOS): grid points inside QC305_MaskFile fail QC305
                           'water' (Yucatan): grid points outside QC305_MaskFile fail QC305
    qc305_remove_vectors:  true (default, MARACOOS): velocities and errors of grid points failing QC305 are removed
                           false (Yucatan): they are only flagged
    uwls_search_distance:  'flat' (default): UWLS radials are selected with the flat earth distance of HFRProgs
                           makeTotals (lonlat2dist.m), so the same radials are used as in MATLAB
                           'geodesic': the hfradarpy default, distances along the WGS84 ellipsoid
    netcdf_missing_integers: 'zero' (default): number_of_radials and the QC flags are 0 where they have no value, as in
                           the production files (realtime_qc_totals_to_netcdf.py casts NaN to integers)
                           'fill': they are set to their _FillValue (-999 and 25)
    netcdf_attributes:     global NetCDF attributes merged over the document's 'netcdf_attributes'

Reprocessing: an hour is only (re)processed when one of its total files does not exist, when the set or the content
of its radial files changed (a SHA-1 checksum of each radial file is recorded in the total files; modification times
are not used because rsync -a preserves them), or when the processing settings changed (a signature of the settings
is recorded in the total files). Use --force to reprocess every hour.

MongoDB bookkeeping (--mongo): after each hour is processed, its record is saved in
totals_generated.<pattern>.<YYYYmmddTHHMMSS> with the same fields as create_hourly_totals.m. Records older than
HoursAgoToCheck + 24 hours are removed. Each total file keeps its record (also without --mongo), so a record missing
from MongoDB (e.g. MongoDB was down, or the hour was processed without --mongo) is restored on the next run with
--mongo. If MongoDB cannot be reached, totals are still processed. The connection string is read from the
HFRADAR_MONGO_URI environment variable (default mongodb://127.0.0.1:27017).

Differences from the MATLAB drivers:
    * MATLAB decides whether to reprocess from the radial count stored in MongoDB; here the decision uses the total
      files themselves (see above). MATLAB also removes a single old record per run, which leaves records behind when
      runs are missed.
    * The Yucatan driver forced the radial month flag for the 'Measured' pattern; here Radials.MonthFlag is used.
    * NetCDF: UWLS GDOP errors are written as unitless values, and qc302_gdop is included in UWLS files. The Yucatan
      converter rounded the unitless GDOP values to integers with a cm2/s2 scale factor, dropped vectors whose rounded
      u_err/v_err/uv_covariance exceeded 1, and left qc302 out of the file although the primary flag uses it.

Example usage:
    # No MongoDB: configuration from an export, no status records
    python realtime_totals.py --config-file codar.totals_configs.json --overrides configs/maracoos_06km_local.json \\
        --start 2026-10-06T00:00 --end 2026-10-07T00:00
    # Operational: configuration from MongoDB and status records for hfradmin
    python realtime_totals.py --mongo --config realtime_06km_maracoos_v03 --overrides configs/maracoos_06km_local.json
    python realtime_totals.py --mongo --config yucatan_realtime_testing --overrides configs/yucatan_local.json
"""
import argparse
import copy
import datetime as dt
import hashlib
import json
import logging
import os
import sys
import time
import warnings
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from shapely.geometry import Point

from hfradarpy.calc import inpolygon
from hfradarpy.io.nc import required_global_attributes
from hfradarpy.radials import Radial
from hfradarpy.totals import combineRadials

logger = logging.getLogger(__name__)
log_level = 'INFO'
log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
logging.basicConfig(stream=sys.stdout, format=log_format, level=log_level)

# pandas emits a lot of FutureWarnings from inside combineRadials
warnings.simplefilter('ignore', category=FutureWarning)

MONGO_URI = os.environ.get('HFRADAR_MONGO_URI', 'mongodb://127.0.0.1:27017')
MONGO_DATABASE = 'codar'
MONGO_COLLECTION = 'totals_configs'
MONGO_KEY_FORMAT = '%Y%m%dT%H%M%S'  # datestr(t, 30) in MATLAB; parsed by hfradmin

DRIVER_DEFAULTS = {'qc305_polygon': 'land', 'qc305_remove_vectors': True, 'uwls_search_distance': 'flat',
                   'netcdf_missing_integers': 'zero', 'netcdf_attributes': {}}

FILE_SUFFIXES = {'RDL': '.ruv', 'ELT': '.euv', 'QCV': '.qcv'}  # filenames_standard_filesystem.m
METHODS = {'uwls': 'UWLS', 'oi': 'OI'}

# NetCDF encoding of the production files
FLAG_FILL = np.int8(25)
FLOAT_FILL = np.float32(-999)
NRAD_FILL = np.int16(-999)
COORD_FILL = np.float32(999)
FLAG_VALUES = np.array([1, 2, 3, 4, 9], dtype=np.int8)
QC_INFO = 'Quality control reference: IOOS QARTOD HF Radar ver 1.0 May 2016\n' \
          'QCFlagDefinitions: 1 = pass, 2 = not_evaluated, 3 = suspect, 4 = fail, 9 = missing_data\n'


# ----------------------------------------------------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------------------------------------------------
def mongo_collection():
    """The MongoDB totals configuration collection (mongo_database.m). pymongo is only needed when MongoDB is used."""
    try:
        from pymongo import MongoClient
    except ImportError:
        raise ImportError('pymongo is required to use MongoDB (--mongo): conda install pymongo') from None
    return MongoClient(MONGO_URI, serverSelectionTimeoutMS=5000)[MONGO_DATABASE][MONGO_COLLECTION]


def load_config(name=None, config_file=None):
    """
    Load a totals configuration document from a JSON file (mongoexport array, one document per line, or a single
    document as returned by hfradmin) or, without config_file, from MongoDB by name.

    :return: configuration document (without _id and totals_generated)
    """
    if config_file:
        text = Path(config_file).read_text().strip()
        try:
            docs = json.loads(text)
        except json.JSONDecodeError:  # one document per line
            docs = [json.loads(line) for line in text.splitlines() if line.strip()]
        docs = docs if isinstance(docs, list) else [docs]
        if name:
            docs = [d for d in docs if d.get('name') == name]
        if len(docs) != 1:
            raise ValueError(f'{config_file} contains {len(docs)} matching configuration documents; use --config to pick one')
        doc = docs[0]
    else:
        doc = mongo_collection().find_one({'name': name}, {'totals_generated': 0})
        if doc is None:
            raise ValueError(f'No "{name}" document in {MONGO_DATABASE}.{MONGO_COLLECTION}')
    doc.pop('_id', None)
    doc.pop('totals_generated', None)
    return doc


def deep_update(base, overrides):
    """Recursively merge the overrides dictionary into base (in place)."""
    for key, value in overrides.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            deep_update(base[key], value)
        else:
            base[key] = value
    return base


def site_list(rconf):
    """
    Site codes and best choice pattern types. Sites are either a list of {code, type} dictionaries (hfradmin) or a
    list of codes with a separate Types list (older documents).
    """
    sites = rconf['Sites']
    if sites and isinstance(sites[0], dict):
        return [s['code'] for s in sites], [s.get('type') or 'RDLm' for s in sites]
    return list(sites), list(rconf['Types'])


def pattern_jobs(rconf):
    """
    Pattern types to process: list of (pattern_field, file type of each site).
    pattern_field is the key of the totals_generated record and the name of the output sub-directory.
    """
    codes, best = site_list(rconf)
    jobs = []
    for pattern, enabled in (rconf.get('ProcessByPattern') or {}).items():
        if not enabled:
            continue
        if pattern == 'Ideal':
            jobs.append((pattern, ['RDLi'] * len(codes)))
        elif pattern in ('Measured', 'MeasuredOnly'):
            jobs.append((pattern, ['RDLm'] * len(codes)))
        elif pattern == 'BestChoice':
            jobs.append((pattern, best))
        else:
            logger.warning(f'Unknown ProcessByPattern entry "{pattern}" ignored')
    return codes, jobs


def enabled_methods(conf):
    return [m for m, key in METHODS.items() if conf['Totals'].get(key, {}).get('Process')]


def processing_signature(conf, driver):
    """
    Short hash of the settings that change the totals. Recorded in each total file so that changing them triggers
    reprocessing. Directories, file names and the realtime window are left out on purpose.
    """
    tconf = conf['Totals']
    settings = {
        'ProcessingVersion': conf.get('ProcessingVersion'),
        'Sites': site_list(conf['Radials']),
        'RadialMask': Path(conf['Radials'].get('MaskDir') or '').name,
        'Grid': Path(tconf['GridFile']).name,
        'MinNumSites': tconf.get('MinNumSites', 2),
        'MinNumRads': tconf.get('MinNumRads', 3),
        'UWLS': {k: tconf['UWLS'].get(k) for k in ('Process', 'spatthresh')} if 'UWLS' in tconf else None,
        'OI': {k: tconf['OI'].get(k) for k in ('Process', 'mdlvar', 'errvar', 'sx', 'sy', 'weighting_function')}
        if 'OI' in tconf else None,
        'QC': {k: Path(v).name if k.endswith('File') else v for k, v in conf['QC'].items()},
        'driver': {k: v for k, v in driver.items() if k != 'netcdf_attributes'},
    }
    return hashlib.sha1(json.dumps(settings, sort_keys=True, default=str).encode()).hexdigest()[:12]


# ----------------------------------------------------------------------------------------------------------------------
# Files
# ----------------------------------------------------------------------------------------------------------------------
def radial_filenames(base_dir, sites, types, t, month_flag=False, type_flag=False):
    """
    Build the file name of every site for one hour (filenames_standard_filesystem.m):
        base_dir/SITE[/TYPE][/YYYY_MM]/TYPE_SITE_YYYY_MM_DD_HHMM.ruv (.euv for ELT types)
    """
    files = []
    for site, ftype in zip(sites, types):
        directory = Path(base_dir) / site
        if type_flag:
            directory = directory / ftype
        if month_flag:
            directory = directory / f'{t:%Y_%m}'
        suffix = FILE_SUFFIXES.get(ftype[:3].upper(), '.ruv')
        files.append(directory / f'{ftype}_{site}_{t:%Y_%m_%d_%H%M}{suffix}')
    return files


def output_filename(conf, method, pattern_field, t):
    """
    Total file name, laid out like the MATLAB drivers and the NetCDF conversion:
        BaseDir/nc/GRIDTYPE/[pattern]/[YYYY_MM]/PREFIXYYYY_MM_DD_HHMM.nc
    GRIDTYPE is conf.GridType (conf.SystemType for older documents). BestChoice totals have no pattern sub-directory.
    """
    tconf = conf['Totals']
    mconf = tconf[METHODS[method]]
    directory = Path(mconf['BaseDir']) / 'nc' / (conf.get('GridType') or conf.get('SystemType') or '')
    if pattern_field != 'BestChoice':
        directory = directory / pattern_field.lower()
    if tconf.get('MonthFlag', True):
        directory = directory / f'{t:%Y_%m}'
    return directory / f'{mconf["FilePrefix"]}{t:%Y_%m_%d_%H%M}.nc'


def load_grid(grid_file):
    """Load a two column (lon, lat) comma or whitespace delimited grid file into a GeoSeries."""
    lonlat = pd.read_csv(grid_file, header=None, names=['lon', 'lat'], sep=r'[,\s]+', engine='python')
    return gpd.GeoSeries([Point(xy) for xy in lonlat.to_numpy()], crs='EPSG:4326')


def load_polygon(polygon_file):
    """Load a two column (lon lat) polygon file, or return None when no file is configured."""
    if not polygon_file:
        return None
    poly = pd.read_csv(polygon_file, header=None, names=['lon', 'lat'], sep=r'[,\s]+', engine='python')
    return poly['lon'].to_numpy(), poly['lat'].to_numpy()


def radial_checksums(files):
    """SHA-1 checksum of the content of each file, keyed by file name."""
    return {f.name: hashlib.sha1(f.read_bytes()).hexdigest() for f in files}


# ----------------------------------------------------------------------------------------------------------------------
# Radials
# ----------------------------------------------------------------------------------------------------------------------
def load_radials(files):
    """
    Load the radial (and elliptical) files that exist and contain data (loadRDLFile_all.m).

    :return: dictionary {site: Radial}, list of the loaded file paths, and the radial_info strings of
             codar_driver_totals.m ('SITE(P): good/total', where good is the number of vectors not flagged PRIM == 4)
    """
    radials, loaded, info = {}, [], []
    for f in files:
        if not f.exists():
            continue
        try:
            r = Radial(f)
        except Exception:
            logger.exception(f'Could not load {f}')
            continue
        if r._iscorrupt or r.data.empty:
            logger.warning(f'Skipping empty or corrupt file {f.name}')
            continue
        site = r.metadata['Site'].split()[0].strip('"')
        radials[site] = r
        loaded.append(str(f))
        pattern = r.metadata.get('PatternType', '?').split()[0][:1].upper()
        good = int((r.data['PRIM'] != 4).sum()) if 'PRIM' in r.data else len(r.data)
        info.append(f'{site:>4s}({pattern}): {good:4d}/{len(r.data):4d}')
    return radials, loaded, info


def clean_radials(radials, remove_flagged, mask):
    """
    Remove radial vectors that cannot be used in the combination, in place:
        * vectors flagged as bad by the radial QC (PRIM == 4) when remove_flagged is true (cleanQCedRadials.m)
        * vectors without a valid velocity or heading (HFRProgs ignores NaN radial velocities)
        * vectors inside the mask polygon (maskRadials.m with InOrOut = false)
    """
    for site, r in radials.items():
        d = r.data
        keep = np.isfinite(d['VELO'].to_numpy(dtype=float)) & np.isfinite(d['HEAD'].to_numpy(dtype=float))
        if remove_flagged and 'PRIM' in d:
            keep &= (d['PRIM'] != 4).to_numpy()
        if mask is not None:
            keep &= ~inpolygon(d['LOND'].to_numpy(dtype=float), d['LATD'].to_numpy(dtype=float), *mask)
        r.data = d.loc[keep].reset_index(drop=True)
        logger.info(f'  {site}: {int(keep.sum())}/{len(keep)} radial vectors used')


# ----------------------------------------------------------------------------------------------------------------------
# Combination and QC
# ----------------------------------------------------------------------------------------------------------------------
def grid_resolution(grid):
    """Approximate grid resolution in meters (metadata only)."""
    lat = np.unique(grid.y.to_numpy())
    return float(np.median(np.diff(lat)) * 111195) if lat.size > 1 else 0.0


def combine(radials, grid, conf, method, t, driver):
    """
    Combine the radials into totals at every grid point (dropPoints=False keeps the grid points without vectors).
    NRAD is counted at every grid point, as HFRProgs TotalsNumRads (makeTotals.m / makeTotalsOI.m count the radials
    before checking MinNumSites and MinNumRads, so grid points without a total vector have a count too).
    """
    tconf = conf['Totals']
    rDF = pd.DataFrame({'Radial': pd.Series(radials, dtype=object)})
    kwargs = dict(gRes=grid_resolution(grid), tStp=t.to_pydatetime(), minContrSites=tconf.get('MinNumSites', 2),
                  minContrRads=tconf.get('MinNumRads', 3), useLandMask=False, dropPoints=False, engine='vectorized',
                  count_all_points=True)
    if method == 'uwls':
        T, warn = combineRadials(rDF, grid, sRad=tconf['UWLS']['spatthresh'] * 1000, method='uwls',
                                 distance=driver['uwls_search_distance'], **kwargs)
    else:
        oi = tconf['OI']
        T, warn = combineRadials(rDF, grid, sRad=2 * max(oi['sx'], oi['sy']) * 1000, method='oi',
                                 mdlvar=float(oi['mdlvar']), errvar=float(oi['errvar']), sx=float(oi['sx']),
                                 sy=float(oi['sy']), oi_option=oi.get('weighting_function') or 'exponential', **kwargs)
    if warn:
        logger.warning(f'{method.upper()}: {warn}')
    return T


def quality_control(T, conf, method, polygon, driver):
    """
    QC tests of codar_driver_totals.m. Flags are computed at every grid point, also where there is no vector.

    :return: names of the tests included in the primary flag, e.g. ['Q303', 'Q305', 'Q302']
    """
    qconf = conf['QC']
    T.initialize_qc()
    T.qc_qartod_maximum_velocity(max_speed=qconf['QC303_MaxTotSpeed'], high_speed=qconf['QC303_MaxTotSpeed'])
    # Without a QC305 mask file the polygon is empty and every grid point passes
    T.qc_valid_location_polygon(*(polygon if polygon is not None else ([], [])),
                                flag_inside=driver['qc305_polygon'] == 'land',
                                name=Path(qconf.get('QC305_MaskFile') or '').name)
    tests = ['Q303', 'Q305']

    if method == 'uwls':
        T.use_hfrprogs_gdop()  # QC302 uses the HFRProgs GDOP definition
        T.qc_qartod_gdop(max_GDOP=qconf['GDOP_UWLS']['QC302_GDOP']['MaxError'])
        tests.append('Q302')
    else:
        gdop_oi = qconf['GDOP_OI']
        T.qc_qartod_u_uncertainty(max_Uerr=gdop_oi['QC306_OI_Uerr']['MaxError'])
        T.qc_qartod_v_uncertainty(max_Verr=gdop_oi['QC307_OI_Verr']['MaxError'])
        tests += ['Q306', 'Q307']
        if 'QC308_OI_Terr' in gdop_oi:
            T.qc_threshold('Q308', 'GDOP', gdop_oi['QC308_OI_Terr']['MaxError'], name='qc_oi_total_errors')
            tests.append('Q308')
        if 'QC309_OI_UVcov' in gdop_oi:
            T.qc_threshold('Q309', 'CQAL', gdop_oi['QC309_OI_UVcov']['MaxError'], absolute=True,
                           name='qc_oi_uv_covariance')
            tests.append('Q309')

    T.qc_qartod_primary_flag()

    if driver['qc305_remove_vectors']:
        # Remove the vectors that fail the valid location test so users who ignore the QC flags do not use them
        bad = T.data['Q305'] == 4
        T.data.loc[bad, ['VELU', 'VELV', 'VELO', 'HEAD', 'UQAL', 'VQAL', 'CQAL', 'GDOP']] = np.nan
    return tests


def vector_stats(T):
    """numvectors, numvectors_passing_qc and num_overwater_qc_failures as computed by codar_driver_totals.m."""
    valid = T.data['VELU'].notna()
    fail = int((valid & (T.data['PRIM'] == 4)).sum())
    fail_location = int((valid & (T.data['Q305'] == 4)).sum())
    return {'numvectors': int(valid.sum()),
            'numvectors_passing_qc': int((valid & (T.data['PRIM'] == 1)).sum()),
            'num_overwater_qc_failures': fail - fail_location}


# ----------------------------------------------------------------------------------------------------------------------
# NetCDF (production format of realtime_qc_totals_to_netcdf.sh)
# ----------------------------------------------------------------------------------------------------------------------
def _flag_attrs(long_name, short_name, comment, northward=False):
    return {'standard_name': f'surface_{"northward" if northward else "eastward"}_sea_water_velocity status_flag',
            'long_name': long_name, 'short_name': short_name, 'units': '1', 'valid_min': np.int8(1),
            'valid_max': np.int8(9), 'coordinates': 'lon lat', 'grid_mapping': 'crs', 'flag_values': FLAG_VALUES,
            'flag_meanings': 'pass not_evaluated suspect fail missing_data', 'comment': comment}


def _err_attrs(short_name, long_name, comment, units='1', valid=(0, 1)):
    return {'short_name': short_name, 'units': units, 'valid_min': np.float32(valid[0]),
            'valid_max': np.float32(valid[1]), 'coordinates': 'lon lat', 'grid_mapping': 'crs',
            'long_name': long_name, 'comment': comment}


def to_production_dataset(T, method, conf, driver, tests, radial_files):
    """
    Build the CF-1.6/NCEI gridded dataset written by the production .mat to NetCDF conversion.

    Velocities are stored in whole cm/s with a 0.01 scale factor (decoded in m/s), QC flags as bytes at every grid
    point, and the processing parameters and radial file names as attributes of container variables.
    """
    qconf, tconf = conf['QC'], conf['Totals']
    d = T.data
    lon = np.unique(d['LOND'].to_numpy(dtype=np.float32))
    lat = np.unique(d['LATD'].to_numpy(dtype=np.float32))
    xi = np.searchsorted(lon, d['LOND'].to_numpy(dtype=np.float32))
    yi = np.searchsorted(lat, d['LATD'].to_numpy(dtype=np.float32))

    def grid(values):
        out = np.full((1, 1, lat.size, lon.size), np.nan)
        out[0, 0, yi, xi] = np.asarray(values, dtype=float)
        return out

    velu, velv = d['VELU'].astype(float), d['VELV'].astype(float)
    if method == 'uwls':
        # GDOP covariance matrix C = inv(A'A) (unitless): u_err = C(1,1), v_err = C(2,2), uv_covariance = C(1,2)
        u_err, v_err = d['UQAL'].astype(float) ** 2, d['VQAL'].astype(float) ** 2
    else:
        # Normalized OI uncertainties xi(1,1), xi(2,2)
        u_err, v_err = d['UQAL'].astype(float), d['VQAL'].astype(float)

    coords = ('time', 'z', 'lat', 'lon')
    ds = xr.Dataset(coords={'time': pd.DatetimeIndex([T.time]), 'z': np.array([0], dtype=np.float32),
                            'lat': lat, 'lon': lon})
    # Round in single precision like the production converter (np.round(U.astype(np.float32)))
    ds['u'] = (coords, grid(np.round(velu.to_numpy(dtype=np.float32))))
    ds['v'] = (coords, grid(np.round(velv.to_numpy(dtype=np.float32))))
    ds['u_err'] = (coords, grid(u_err))
    ds['v_err'] = (coords, grid(v_err))
    ds['total_errors'] = (coords, grid(d['GDOP'].astype(float)))
    ds['uv_covariance'] = (coords, grid(d['CQAL'].astype(float)))
    ds['number_of_radials'] = (coords, grid(d['NRAD'].astype(float)))

    names = {'Q302': 'qc302_gdop', 'Q303': 'qc303_maxspeed', 'Q305': 'qc305_validlocation', 'Q306': 'qc306_uerr',
             'Q307': 'qc307_verr', 'Q308': 'qc308_terr', 'Q309': 'qc309_uvcov'}
    for test in tests:
        ds[names[test]] = (coords, grid(d[test]))
    ds['qc_primary_flag'] = (coords, grid(d['PRIM']))
    ds['qc_operator_flag'] = (coords, grid(np.full(len(d), 2)))  # operator flag: not evaluated
    if driver['netcdf_missing_integers'] == 'zero':
        for name in ['number_of_radials', 'qc_primary_flag', 'qc_operator_flag', *(names[t] for t in tests)]:
            ds[name] = ds[name].fillna(0)
    ancillary = ' '.join(names[t] for t in tests)

    vel_attrs = {'units': 'm/s', 'scale_factor': np.float32(0.01), 'valid_min': np.float32(-300),
                 'valid_max': np.float32(300), 'coordinates': 'lon lat', 'grid_mapping': 'crs',
                 'ancillary_variables': f'{ancillary} qc_primary_flag qc_operator_flag'}
    ds['u'].attrs = {'long_name': 'Eastward Surface Current (m/s)', 'standard_name': 'surface_eastward_sea_water_velocity',
                     'short_name': 'u', **vel_attrs}
    ds['v'].attrs = {'long_name': 'Northward Surface Current (m/s)',
                     'standard_name': 'surface_northward_sea_water_velocity', 'short_name': 'v', **vel_attrs}

    if method == 'oi':
        uerr_max = qconf['GDOP_OI']['QC306_OI_Uerr']['MaxError']
        verr_max = qconf['GDOP_OI']['QC307_OI_Verr']['MaxError']
        ds['u_err'].attrs = _err_attrs('uerr', 'Normalized uncertainty associated with eastward velocity component, '
                                       '<(u_hat - u)^2>/<u^2>', 'normalized by constant signal covariance for all '
                                       f'gridpoints, velocity measurements with u_err values over {uerr_max} are of '
                                       'questionable quality')
        ds['v_err'].attrs = _err_attrs('verr', 'Normalized uncertainty associated with northward velocity component, '
                                       '<(v_hat - v)^2>/<v^2>', 'normalized by constant signal covariance for all '
                                       f'gridpoints, velocity measurements with v_err values over {verr_max} are of '
                                       'questionable quality')
        ds['total_errors'].attrs = _err_attrs('totalerr', 'Associated mapping error value, sqrt(U_err^2 + Verr^2)',
                                              'velocity measurements with high total_error values are of questionable '
                                              'quality', valid=(-300, 300))
        ds['uv_covariance'].attrs = _err_attrs('uvcov', 'Directional information of u and v, '
                                               '<(u_hat -u)(v_hat- v)>/sqrt(<u^2><v^2>)',
                                               'directional information of u and v', valid=(-300, 300))
    else:
        gdop_max = qconf['GDOP_UWLS']['QC302_GDOP']['MaxError']
        ds['u_err'].attrs = _err_attrs('uerr', 'Associated GDOP mapping error value associated with eastward velocity '
                                       'component', 'unitless; velocity measurements with high u_err values are of '
                                       'questionable quality', valid=(0, 300))
        ds['v_err'].attrs = _err_attrs('verr', 'Associated GDOP mapping error value associated with northward velocity '
                                       'component', 'unitless; velocity measurements with high v_err values are of '
                                       'questionable quality', valid=(0, 300))
        ds['total_errors'].attrs = _err_attrs('totalerr', 'Associated GDOP mapping error value, TotalErrors are the '
                                              'square-root of the norm covariance matrix (or the square root of the '
                                              'maximum eigenvalue of the 2x2 UV covariance matrix)',
                                              f'velocity measurements with total_error values over {gdop_max} are of '
                                              'questionable quality', valid=(0, 300))
        ds['uv_covariance'].attrs = _err_attrs('uvcov', 'Associated GDOP mapping error value, eastward and northward '
                                               'covariance', 'unitless; directional information of u and v',
                                               valid=(-300, 300))

    ds['number_of_radials'].attrs = {'long_name': 'Number of radial measurements used to calculate each totals velocity',
                                     'short_name': 'num_rads', 'units': '1',
                                     'comment': f'totals are not calculated with fewer than '
                                                f'{tconf.get("MinNumRads", 3)} contributing radial measurements from '
                                                f'{tconf.get("MinNumSites", 2)} sites',
                                     'coordinates': 'lon lat', 'grid_mapping': 'crs'}

    flag_comments = {
        'Q302': ('GDOP flag', f'qc302 GDOP Threshold [GDOP TotalErrors '
                              f'{qconf.get("GDOP_UWLS", {}).get("QC302_GDOP", {}).get("MaxError")}]'),
        'Q303': ('maximum speed test flag', f'qc303 Max Speed Threshold [max_vel = {qconf["QC303_MaxTotSpeed"]} (cm/s)]'),
        'Q305': ('valid location test flag', f'qc305 Valid Location [landmask file = {qconf.get("QC305_MaskFile")}]'),
        'Q306': ('Uerr flag', 'qc306 OI Uncertainty Threshold [OIuncert Uerr '
                              f'{qconf.get("GDOP_OI", {}).get("QC306_OI_Uerr", {}).get("MaxError")}]'),
        'Q307': ('Verr flag', 'qc307 OI Uncertainty Threshold [OIuncert Verr '
                              f'{qconf.get("GDOP_OI", {}).get("QC307_OI_Verr", {}).get("MaxError")}]'),
        'Q308': ('TotalErrors flag', 'qc308 OI Uncertainty Threshold [OIuncert TotalErrors '
                                     f'{qconf.get("GDOP_OI", {}).get("QC308_OI_Terr", {}).get("MaxError")}]'),
        'Q309': ('UVCovariance flag', 'qc309 OI Uncertainty Threshold [OIuncert abs(UVCovariance) '
                                      f'{qconf.get("GDOP_OI", {}).get("QC309_OI_UVcov", {}).get("MaxError")}]'),
    }
    for test in tests:
        long_name, comment = flag_comments[test]
        ds[names[test]].attrs = _flag_attrs(long_name, names[test].split('_')[0], QC_INFO + comment,
                                            northward=test == 'Q307')
    listed = ', '.join(names[t].split('_')[0] for t in tests)
    ds['qc_primary_flag'].attrs = _flag_attrs('Primary Flag', 'primaryflag',
                                              f'QCPrimaryFlagDefinition: Highest flag value of {listed}\nThis flag will '
                                              'be set to not_evaluated only if ALL individual tests were not_evaluated.')
    ds['qc_primary_flag'].attrs['ancillary_variables'] = ancillary
    ds['qc_operator_flag'].attrs = _flag_attrs('Operator Flag', 'opflag',
                                               QC_INFO + 'The qc_operator_flag follows QCFlagDefinitions and is set at '
                                               'discretion of the operator or data manager.')

    ds['time'].attrs['standard_name'] = 'time'
    ds['lon'].attrs = {'long_name': 'Longitude', 'standard_name': 'longitude', 'short_name': 'lon',
                       'units': 'degrees_east', 'axis': 'X', 'valid_min': np.float32(-180), 'valid_max': np.float32(180)}
    ds['lat'].attrs = {'long_name': 'Latitude', 'standard_name': 'latitude', 'short_name': 'lat',
                       'units': 'degrees_north', 'axis': 'Y', 'valid_min': np.float32(-90), 'valid_max': np.float32(90)}
    ds['z'].attrs = {'long_name': 'Average Depth of Sensor', 'standard_name': 'depth',
                     'comment': 'Derived from mean value of depth variable', 'units': 'm', 'axis': 'Z',
                     'positive': 'down'}

    # Container variables
    names_loaded = [Path(f).name for f in radial_files]
    ds['crs'] = ((), np.int8(0), {'grid_mapping_name': 'latitude_longitude', 'inverse_flattening': 298.257223563,
                                  'long_name': 'Coordinate Reference System', 'semi_major_axis': '6378137.0',
                                  'epsg_code': 'EPSG:4326', 'comment': 'http://www.opengis.net/def/crs/EPSG/0/4326'})
    ds['instrument'] = ((), np.int8(0), {'long_name': 'CODAR SeaSonde High Frequency Radar',
                                         'sensor_type': 'Direction-finding high frequency radar antenna',
                                         'make_model': 'CODAR SeaSonde', 'serial_number': 1})
    ds['radial_metadata'] = ((), np.int8(0), {'long_name': 'Metadata on radial velocities used to compute total solutions',
                                              'short_name': 'radmeta', 'number_files_loaded': str(len(names_loaded)),
                                              'number_files_loaded_description': 'Number of radial files loaded',
                                              'files_loaded': ','.join(names_loaded),
                                              'files_loaded_description': 'Radial file names loaded'})
    params = [(np.int8(tconf.get('MinNumSites', 2)), 'Minimum number of radial sites'),
              (np.int8(tconf.get('MinNumRads', 3)), 'Minimum number of radial vectors'),
              (np.float32(tconf[METHODS[method]].get('tempthresh', 1 / 48)),
               'Temporal search window for radial solutions (Fraction of a day)')]
    if method == 'oi':
        oi = tconf['OI']
        params += [(np.float32(oi['sx']), 'Decorrelation scales in the east direction'),
                   (np.float32(oi['sy']), 'Decorrelation scales in the north direction'),
                   (np.int16(oi['mdlvar']), 'Signal variance of the surface current fields (cm2 s-2)'),
                   (np.int16(oi['errvar']), 'Data error variance of the input radial velocities (cm2 s-2)')]
    else:
        params += [(np.float32(tconf['UWLS']['spatthresh']), 'Spatial search radius for radial solutions (km)')]
    pattrs = {'units': '1', 'long_name': 'General and method specific processing parameter information',
              'short_name': 'params'}
    for k, (value, description) in enumerate(params, start=1):
        pattrs[f'parameter{k}'] = value
        pattrs[f'parameter{k}_description'] = description
    ds['processing_parameters'] = ((), np.int8(0), pattrs)

    # Global attributes
    attrs = dict(conf.get('netcdf_attributes') or {})
    attrs.update(driver.get('netcdf_attributes') or {})
    attrs['method'] = 'Optimal Interpolation' if method == 'oi' else 'Unweighted Least Squares'
    if attrs:
        t = pd.Timestamp(T.time)
        gattrs = required_global_attributes(attrs, t, t)
        for key in ('date_created', 'date_modified', 'date_issued', 'date_metadata_modified', 'time_coverage_start',
                    'time_coverage_end'):  # production date format
            gattrs[key] = pd.Timestamp(gattrs[key]).strftime('%Y%m%dT%H%M%SZ')
        gattrs.update({k: v for k, v in attrs.items() if k not in gattrs})
    else:
        gattrs = {'method': attrs['method']}
    gattrs['geospatial_lat_min'], gattrs['geospatial_lat_max'] = lat.min(), lat.max()
    gattrs['geospatial_lon_min'], gattrs['geospatial_lon_max'] = lon.min(), lon.max()
    ds.attrs = gattrs
    return ds


def write_netcdf(ds, fname, extra_attrs):
    """Write a production dataset with the production encoding."""
    ds = ds.copy()
    ds.attrs.update(extra_attrs)
    comp = dict(zlib=True, complevel=4, shuffle=True)
    chunks = (1, 1, ds.sizes['lat'], ds.sizes['lon'])
    encoding = {}
    for name, var in ds.data_vars.items():
        if var.ndim == 0:
            encoding[name] = {'_FillValue': None}
        elif name.startswith('qc'):
            encoding[name] = dict(dtype='int8', _FillValue=FLAG_FILL, chunksizes=chunks, **comp)
        elif name == 'number_of_radials':
            encoding[name] = dict(dtype='int16', _FillValue=NRAD_FILL, chunksizes=chunks, **comp)
        else:
            encoding[name] = dict(dtype='float32', _FillValue=FLOAT_FILL, chunksizes=chunks, **comp)
    for name in ('lon', 'lat', 'z'):
        encoding[name] = dict(dtype='float32', _FillValue=COORD_FILL)
    encoding['time'] = dict(units='seconds since 1970-01-01', calendar='gregorian', dtype='float64', _FillValue=None)
    fname.parent.mkdir(parents=True, exist_ok=True)
    tmp = fname.with_suffix('.nc.tmp')
    ds.to_netcdf(tmp, encoding=encoding, format='NETCDF4', engine='netcdf4', unlimited_dims=['time'])
    tmp.replace(fname)  # never leave a partial file behind for the reprocessing check or for users


# ----------------------------------------------------------------------------------------------------------------------
# Reprocessing check and MongoDB bookkeeping
# ----------------------------------------------------------------------------------------------------------------------
def needs_processing(fnames, checksums, signature):
    """
    Decide whether an hour has to be (re)processed by comparing its total files with the radial files available now.

    :param fnames: total files of the hour (one per method)
    :param checksums: radial_checksums of the radial files of the hour that exist
    :param signature: processing_signature of the current configuration
    :return: reason for (re)processing, or an empty string if the total files are up to date
    """
    for fname in fnames:
        if not fname.exists():
            return f'{fname.name} does not exist'
        try:
            with xr.open_dataset(fname) as ds:
                recorded = json.loads(ds.attrs.get('radial_checksums', '{}'))
                old_signature = ds.attrs.get('processing_signature', '')
        except Exception:
            return f'{fname.name} could not be read'
        if old_signature != signature:
            return 'the processing configuration changed'
        names, old_names = set(checksums), set(recorded)
        if names != old_names:
            return f'available radials changed from {len(old_names)} to {len(names)} ' \
                   f'({", ".join(sorted(names ^ old_names))})'
        changed = sorted(n for n in names if checksums[n] != recorded[n])
        if changed:
            return f'radial files changed since the totals were made ({", ".join(changed)})'
    return ''


def read_totals_generated(fname):
    """Read the totals_generated record stored in a total file, or None if it has none."""
    with xr.open_dataset(fname) as ds:
        record = ds.attrs.get('totals_generated')
    return json.loads(record) if record else None


def mongo_connect(name):
    """
    Connect to the totals configuration collection (mongo_database.m) and read the keys of the existing
    totals_generated records of the configuration document.

    :return: (collection, records) where records is {pattern: set of record keys}, or (None, None) if MongoDB cannot
             be reached or the configuration document does not exist
    """
    try:
        collection = mongo_collection()
        doc = collection.find_one({'name': name}, {'totals_generated': 1})
    except Exception:
        logger.exception('Could not connect to MongoDB. Totals will be processed without updating totals_generated')
        return None, None
    if doc is None:
        logger.warning(f'MongoDB: no "{name}" document in {MONGO_DATABASE}.{MONGO_COLLECTION}. '
                       f'Totals will be processed without updating totals_generated')
        return None, None
    records = {pattern: set(recs) for pattern, recs in (doc.get('totals_generated') or {}).items()
               if isinstance(recs, dict)}
    return collection, records


def mongo_record(collection, name, pattern_field, t, generated):
    """Store the totals_generated record of one hour, with the fields written by create_hourly_totals.m."""
    field = f'totals_generated.{pattern_field}.{t:{MONGO_KEY_FORMAT}}'
    update = {f'{field}.num_radials': generated['numfiles'],
              f'{field}.files': generated['files'],
              f'{field}.radial_info': generated['radial_info']}
    for method in METHODS:
        stats = generated[method]
        if stats['numvectors'] > 0:
            update.update({f'{field}.{method}_num_vectors': stats['numvectors'],
                           f'{field}.{method}_num_vectors_passing_qc': stats['numvectors_passing_qc'],
                           f'{field}.{method}_num_overwater_qc_failures': stats['num_overwater_qc_failures']})
    collection.update_one({'name': name}, {'$set': update})


def mongo_cleanup(collection, name, records, patterns, t, hours_ago_to_check):
    """
    Remove the records more than HoursAgoToCheck + 24 hours older than t so the totals_generated field does not get
    too large. Unlike the MATLAB drivers, which remove a single record per run, every older record is removed so
    missed runs do not leave records behind.
    """
    cutoff = t - pd.Timedelta(hours=hours_ago_to_check + 24)
    fields = []
    for pattern in patterns:
        for key in sorted(records.get(pattern, ())):
            try:
                old = dt.datetime.strptime(key, MONGO_KEY_FORMAT) < cutoff
            except ValueError:
                continue
            if old:
                fields.append(f'totals_generated.{pattern}.{key}')
                records[pattern].discard(key)
    if fields:
        collection.update_one({'name': name}, {'$unset': {field: '' for field in fields}})
    logger.info(f'Cleaning totals_generated records older than {cutoff}: removed {len(fields)}.')


# ----------------------------------------------------------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------------------------------------------------------
def create_hourly_totals(t, conf, driver, files, outputs, grid, polygons, signature, checksums):
    """
    Generate, quality control and save the totals of one hour for every enabled method
    (create_hourly_totals.m + codar_driver_totals.m).

    :param outputs: {method: total file name}
    :return: totals_generated record of the hour
    """
    generated = {'files': [], 'numfiles': 0, 'radial_info': [],
                 **{m: {'numvectors': 0, 'numvectors_passing_qc': 0, 'num_overwater_qc_failures': 0}
                    for m in METHODS}}

    logger.info('Loading Radial Data.')
    radials, loaded, info = load_radials(files)
    generated.update(files=loaded, numfiles=len(loaded), radial_info=info)
    logger.info(f'{len(loaded)} radial files loaded: {", ".join(info)}')
    if len(radials) < conf['Totals'].get('MinNumSites', 2):
        logger.info('Not enough radials to generate totals.')
        return generated

    clean_radials(radials, remove_flagged=conf['QC'].get('Process', False), mask=polygons['radials'])

    datasets = {}
    for method, fname in outputs.items():
        logger.info(f'Generating {method.upper()} totals.')
        T = combine(copy.deepcopy(radials), grid, conf, method, t, driver)
        if T.data.empty or T.data['VELU'].notna().sum() == 0:
            logger.info(f'{method.upper()}: no total vectors')
            continue
        tests = quality_control(T, conf, method, polygons['totals'], driver)
        stats = vector_stats(T)
        generated[method] = stats
        logger.info(f'{method.upper()}: {stats["numvectors"]} vectors, {stats["numvectors_passing_qc"]} passing QC, '
                    f'{stats["num_overwater_qc_failures"]} over water QC failures')
        datasets[fname] = to_production_dataset(T, method, conf, driver, tests, loaded)

    # Write the files once every method is done, so each one keeps the complete record of the hour
    for fname, ds in datasets.items():
        write_netcdf(ds, fname, {'radial_checksums': json.dumps(checksums), 'processing_signature': signature,
                                 'totals_generated': json.dumps(generated)})
        logger.info(f'Saved {fname}')
    return generated


def main(conf, driver, start, end, force=False, use_mongo=False):
    init_time = time.time()
    logger.info(f'realtime_totals.py ({conf["name"]}) Start Time: {dt.datetime.now():%Y-%m-%d %H:%M:%S}')

    rconf = conf['Radials']
    methods = enabled_methods(conf)
    if not methods:
        logger.warning('Neither Totals.UWLS.Process nor Totals.OI.Process is enabled; nothing to do')
        return
    grid = load_grid(conf['Totals']['GridFile'])
    polygons = {'radials': load_polygon(rconf.get('MaskDir')), 'totals': load_polygon(conf['QC'].get('QC305_MaskFile'))}
    signature = processing_signature(conf, driver)
    sites, jobs = pattern_jobs(rconf)

    # Hourly timestamps, newest first
    time_steps = pd.date_range(start, end, freq='h')[::-1]

    collection, records = mongo_connect(conf['name']) if use_mongo else (None, None)
    if collection is not None:
        try:
            mongo_cleanup(collection, conf['name'], records, [j[0] for j in jobs], time_steps[0],
                          rconf.get('HoursAgoToCheck', 168))
        except Exception:
            logger.exception('Could not clean the totals_generated records')

    def save_record(pattern_field, t, generated):
        try:
            mongo_record(collection, conf['name'], pattern_field, t, generated)
            records.setdefault(pattern_field, set()).add(f'{t:{MONGO_KEY_FORMAT}}')
        except Exception:
            logger.exception(f'Could not save the totals_generated record of {pattern_field} {t} in MongoDB')

    summary = []
    for t in time_steps:
        logger.info('****************************************')
        logger.info(f'  Current time: {dt.datetime.now():%Y-%m-%d %H:%M:%S}')
        logger.info(f'  Processing data time: {t:%Y-%m-%d %H:%M:%S}')
        logger.info('****************************************')

        for pattern_field, types in jobs:
            files = radial_filenames(rconf['BaseDir'], sites, types, t, rconf.get('MonthFlag', False),
                                     rconf.get('TypeFlag', False))
            available = [f for f in files if f.exists()]
            outputs = {m: output_filename(conf, m, pattern_field, t) for m in methods}
            row = {'time': t, 'pattern': pattern_field, 'available': len(available)}

            if len(available) < conf['Totals'].get('MinNumSites', 2):
                logger.info(f'{pattern_field}: {len(available)} radial file(s) available. Not enough radials to '
                            f'generate totals. Continuing to next timestep.')
                summary.append({**row, 'action': 'not enough radials'})
                continue

            checksums = radial_checksums(available)
            reason = 'forced' if force else needs_processing(outputs.values(), checksums, signature)
            if not reason:
                logger.info(f'{pattern_field}: totals are up to date and do NOT need to be re-processed.')
                action = 'up to date'
                if collection is not None and f'{t:{MONGO_KEY_FORMAT}}' not in records.get(pattern_field, ()):
                    generated = read_totals_generated(next(iter(outputs.values())))
                    if generated:
                        logger.info(f'{pattern_field}: restoring the missing totals_generated record')
                        save_record(pattern_field, t, generated)
                        action = 'up to date, record restored'
                summary.append({**row, 'action': action})
                continue

            logger.info(f'{pattern_field}: making hourly totals, {reason}.')
            try:
                generated = create_hourly_totals(t, conf, driver, files, outputs, grid, polygons, signature, checksums)
            except Exception:
                logger.exception('Total Creation failed')
                summary.append({**row, 'action': 'failed'})
                continue
            if collection is not None:
                save_record(pattern_field, t, generated)
            summary.append({**row, 'action': 'processed', 'loaded': generated['numfiles'],
                            **{f'{m}_vectors': generated[m]['numvectors'] for m in methods},
                            **{f'{m}_passing': generated[m]['numvectors_passing_qc'] for m in methods}})

    logger.info('Total Creation Finished.')
    if summary:
        with pd.option_context('display.width', 250, 'display.max_columns', 20, 'display.max_rows', 1000):
            print('\n' + pd.DataFrame(summary).set_index(['time', 'pattern']).sort_index().to_string())

    elapsed = time.time() - init_time
    logger.info(f'Elapsed time is {int(elapsed // 60)} minutes and {elapsed % 60:.0f} seconds.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', help='Name of the totals_configs document (required with --mongo and no --config-file)')
    parser.add_argument('--config-file', type=Path, help='JSON file with the totals_configs document(s)')
    parser.add_argument('--overrides', type=Path,
                        help='JSON file merged over the document: {"conf": {...}, "driver": {...}, ...}')
    parser.add_argument('--start', help='First timestamp (UTC), e.g. 2026-10-06T00:00 (default: realtime window)')
    parser.add_argument('--end', help='Last timestamp (UTC) (default: the end of the current hour)')
    parser.add_argument('--hours-ago', type=int,
                        help='Process this many hours before --end (default: Radials.HoursAgoToCheck)')
    parser.add_argument('--force', action='store_true', help='Reprocess every hour, even if its totals are up to date')
    parser.add_argument('--mongo', action='store_true',
                        help='Use MongoDB (HFRADAR_MONGO_URI): read the configuration from it when there is no '
                             '--config-file, and save the totals_generated status records read by hfradmin')
    args = parser.parse_args()
    if not args.config_file and not (args.mongo and args.config):
        parser.error('use --config-file, or --mongo with --config to read the configuration from MongoDB')
    if args.mongo:
        try:
            import pymongo  # noqa: F401
        except ImportError:
            parser.error('--mongo requires pymongo: conda install pymongo')

    doc = load_config(args.config, args.config_file)
    driver_opts = copy.deepcopy(DRIVER_DEFAULTS)
    if args.overrides:
        overrides = json.loads(args.overrides.read_text())
        deep_update(driver_opts, overrides.pop('driver', {}))
        deep_update(doc, overrides)
    if driver_opts['qc305_polygon'] not in ('land', 'water'):
        parser.error('driver.qc305_polygon must be "land" or "water"')
    if driver_opts['netcdf_missing_integers'] not in ('zero', 'fill'):
        parser.error('driver.netcdf_missing_integers must be "zero" or "fill"')
    if driver_opts['uwls_search_distance'] not in ('flat', 'geodesic'):
        parser.error('driver.uwls_search_distance must be "flat" or "geodesic"')

    end = pd.Timestamp(args.end) if args.end else pd.Timestamp.now(tz='UTC').tz_localize(None).ceil('h')
    if args.start:
        start = pd.Timestamp(args.start)
    else:
        start = end - pd.Timedelta(hours=args.hours_ago or doc['conf']['Radials'].get('HoursAgoToCheck', 168))

    conf = doc['conf']
    conf['name'] = doc['name']
    conf['netcdf_attributes'] = doc.get('netcdf_attributes') or {}
    main(conf, driver_opts, start, end, force=args.force, use_mongo=args.mongo)
