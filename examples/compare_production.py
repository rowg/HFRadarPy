#!/usr/bin/env python
"""
@author Mike Smith
@email michaesm@marine.rutgers.edu
@purpose Check the output of realtime_totals.py against the production MARACOOS totals (HFRProgs makeTotals and
         makeTotalsOI run by codar_driver_totals.m, then realtime_qc_totals_to_netcdf.py), value for value.

Expected layout of the production directory (copied from boardwalk with rsync -rt):
    <prod>/oi_nc/hfr_rtv_midatl_06km_oi_maracoos_YYYY_MM_DD_HH00.nc      published OI NetCDF
    <prod>/lsq_mat/hfr_rtv_midatl_06km_uwls_maracoos_YYYY_MM_DD_HH00.mat  UWLS TUV structures (no UWLS NetCDF is made)
    <prod>/totals_configs_export.json                                     mongoexport of the totals_configs document,
                                                                          for the hourly totals_generated records

<ours> is the output of realtime_totals.py: the UWLS and OI BaseDir of the overrides file, <ours>/lsq and <ours>/oi.

Compared:
    OI    every variable of the production NetCDF as stored (velocities in whole cm/s), its attributes and the global
          attributes
    UWLS  u and v rounded as the production converter rounds them (np.round(float32)), the GDOP errors to 1 part in
          1e5, number_of_radials and the QC flags at every grid point
    totals_generated records (radial files, radial summary and vector counts) where production has one for the hour

With --baseline, every variable and attribute of each <ours> file is also compared exactly with the same file in an
earlier driver output, e.g. the run that matched production, to check that a code change did not alter the output.

Example usage:
    python compare_production.py ~/Documents/hfradar_parity/2026-10-06_week/prod ~/Desktop/prodcheck
    python compare_production.py PROD OURS --baseline ~/Documents/hfradar_parity/2026-10-06_week/baseline
    python compare_production.py PROD OURS --hours 2026_10_12   # only the hours starting with this stamp

Exit status is 0 when nothing differs, 1 otherwise.
"""
import argparse
import collections
import json
import re
import sys
from pathlib import Path

import netCDF4
import numpy as np
import scipy.io as sio

# Global attributes that change with every run
RUN_ATTRIBUTES = {'date_created', 'date_modified', 'date_issued', 'date_metadata_modified', 'history',
                  'time_coverage_start', 'time_coverage_end', 'radial_checksums', 'processing_signature',
                  'totals_generated'}

# Variable attributes that legitimately differ between machines: the comment quotes the local mask file path
MACHINE_ATTRIBUTES = {'qc305_validlocation.comment'}

# OI floating point values (e.g. uv_covariance) may differ in the last bits from the order of summation
OI_FLOAT_TOLERANCE = 1e-12

UWLS_QC = {'qc303_maxspeed': 'QC303', 'qc305_validlocation': 'QC305', 'qc302_gdop': 'QC302',
           'qc_primary_flag': 'PRIM', 'qc_operator_flag': 'qc_operator_flag'}
UWLS_GDOP = {'u_err': 'Uerr', 'v_err': 'Verr', 'uv_covariance': 'UVCovariance', 'total_errors': 'TotalErrors'}
RECORD_COUNTS = [('uwls', 'numvectors'), ('uwls', 'numvectors_passing_qc'), ('uwls', 'num_overwater_qc_failures'),
                 ('oi', 'numvectors'), ('oi', 'numvectors_passing_qc'), ('oi', 'num_overwater_qc_failures')]


class Tally:
    """Count the compared and differing values of each variable over all hours."""

    def __init__(self):
        self.stats = collections.defaultdict(lambda: dict(hours=0, values=0, differ=0, max_diff=0.0, first=[]))
        self.notes = collections.Counter()

    def values(self, key, ours, ref, stamp, tol=0.0):
        """Compare two arrays; NaN in both counts as equal, NaN in only one as a difference."""
        ours, ref = np.asarray(ours, dtype=float), np.asarray(ref, dtype=float)
        both_nan = np.isnan(ours) & np.isnan(ref)
        diff = np.where(both_nan, 0, np.abs(ours - ref))
        diff = np.where(np.isnan(diff), np.inf, diff)
        bad = diff > tol
        s = self.stats[key]
        s['hours'] += 1
        s['values'] += int((~both_nan).sum())
        s['differ'] += int(bad.sum())
        if bad.any():
            s['max_diff'] = max(s['max_diff'], float(diff[bad].max()))
            if len(s['first']) < 3:
                s['first'].append((stamp, int(bad.sum()), float(diff[bad].max())))

    def note(self, text):
        """Record a structural or attribute difference."""
        self.notes[text] += 1

    def ok(self):
        return not self.notes and all(s['differ'] == 0 for s in self.stats.values())

    def report(self, title):
        print(f'\n{title}')
        print(f'  {"variable":40s} {"hours":>5s} {"values":>10s} {"differ":>7s} {"max diff":>9s}  first (hour, values, max)')
        for (group, name), s in self.stats.items():
            print(f'  {group + " " + name:40s} {s["hours"]:5d} {s["values"]:10d} {s["differ"]:7d} '
                  f'{s["max_diff"]:9.3g}  {s["first"]}')
        if self.notes:
            print('  structure and attribute differences (number of hours):')
            for text, n in sorted(self.notes.items()):
                print(f'    {n:4d}  {text}')


def stored(nc, name):
    """A variable's values as stored in the file (no scale factor), with its fill value as NaN."""
    var = nc.variables[name]
    var.set_auto_maskandscale(False)
    values = np.asarray(var[:], dtype=float).squeeze()
    if '_FillValue' in var.ncattrs():
        values = np.where(values == var.getncattr('_FillValue'), np.nan, values)
    return values


def compare_attributes(tally, label, ref, ours, skip=()):
    """Note every attribute that is missing from one object or differs between them."""
    for k in sorted(set(ref.ncattrs()) | set(ours.ncattrs())):
        if k in skip:
            continue
        a = ref.getncattr(k) if k in ref.ncattrs() else None
        b = ours.getncattr(k) if k in ours.ncattrs() else None
        if str(a) != str(b):
            tally.note(f'{label}{k}: reference={str(a)[:60]!r} ours={str(b)[:60]!r}')


def compare_netcdf(tally, group, ref_file, our_file, stamp, tol=0.0):
    """Compare every variable and attribute of two NetCDF files, values as stored."""
    with netCDF4.Dataset(ref_file) as ref, netCDF4.Dataset(our_file) as ours:
        for name, rv in ref.variables.items():
            if name not in ours.variables:
                tally.note(f'{group} variable missing from ours: {name}')
                continue
            ov = ours.variables[name]
            if rv.dtype != ov.dtype or rv.shape != ov.shape:
                tally.note(f'{group} {name}: reference {rv.dtype}{rv.shape}, ours {ov.dtype}{ov.shape}')
                continue
            skip = {k.split('.', 1)[1] for k in MACHINE_ATTRIBUTES if k.startswith(f'{name}.')}
            compare_attributes(tally, f'{group} {name}.', rv, ov, skip)
            if rv.dtype.kind in 'fiu':
                tolerance = tol if rv.dtype.kind == 'f' else 0.0
                tally.values((group, name), stored(ours, name), stored(ref, name), stamp, tolerance)
        for name in ours.variables:
            if name not in ref.variables:
                tally.note(f'{group} extra variable in ours: {name}')
        compare_attributes(tally, f'{group} global ', ref, ours, RUN_ATTRIBUTES)


def compare_uwls_mat(tally, mat_file, our_file, stamp):
    """Compare the driver's UWLS NetCDF with a production TUV structure at the grid points of the structure."""
    tuv = sio.loadmat(mat_file, squeeze_me=True, struct_as_record=False)['TUV']
    with netCDF4.Dataset(our_file) as nc:
        lat = nc.variables['lat'][:].astype(float)
        lon = nc.variables['lon'][:].astype(float)
        ii = np.abs(lat[None, :] - tuv.LonLat[:, 1, None]).argmin(1)
        jj = np.abs(lon[None, :] - tuv.LonLat[:, 0, None]).argmin(1)
        if max(np.abs(lat[ii] - tuv.LonLat[:, 1]).max(), np.abs(lon[jj] - tuv.LonLat[:, 0]).max()) > 1e-3:
            tally.note('UWLS grid does not match the production grid')
            return
        at_grid = lambda name: stored(nc, name)[ii, jj]

        for name, ref in [('u', tuv.U), ('v', tuv.V)]:
            # realtime_qc_totals_to_netcdf.py rounds after converting to single precision
            tally.values(('uwls', name), at_grid(name), np.round(np.asarray(ref, np.float32)), stamp)
        gdop = {e.Type: e for e in np.atleast_1d(tuv.ErrorEstimates)}['GDOP']
        for name, field in UWLS_GDOP.items():
            ref = np.asarray(getattr(gdop, field), float)
            tally.values(('uwls', name), at_grid(name), ref, stamp, tol=1e-5 * max(1, np.nanmax(np.abs(ref))))
        num_rads = np.nan_to_num(np.asarray(tuv.OtherMatrixVars.makeTotals_TotalsNumRads, float))
        tally.values(('uwls', 'number_of_radials'), at_grid('number_of_radials'), num_rads, stamp)
        for name, field in UWLS_QC.items():
            tally.values(('uwls', name), at_grid(name), getattr(tuv, field), stamp)

        outside = np.ones((lat.size, lon.size), bool)
        outside[ii, jj] = False
        if any(np.isfinite(stored(nc, name)[outside]).any() for name in ('u', 'v')):
            tally.note('UWLS vectors outside the production grid')


def compare_records(tally, export_file, ours, stamps):
    """Compare the totals_generated record in each UWLS file with production's record for the same hour."""
    docs = json.loads(Path(export_file).read_text())
    docs = docs if isinstance(docs, list) else [docs]
    prod = docs[0]['totals_generated']['BestChoice']
    normalize = lambda info: sorted(re.sub(r'\s+', ' ', x) for x in info)
    for stamp in stamps:
        y, m, d, hm = stamp.split('_')
        key = f'{y}{m}{d}T{hm}00'
        if key not in prod:
            continue
        p = prod[key]
        with netCDF4.Dataset(find_file(ours / 'lsq', f'*_uwls_*_{stamp}.nc')) as nc:
            o = json.loads(nc.getncattr('totals_generated'))
        same = {
            'radial files read': sorted(Path(f).name for f in p['files']) == sorted(Path(f).name for f in o['files']),
            'number of radial files': p['num_radials'] == o['numfiles'],
            'radial summary per site': normalize(p['radial_info']) == normalize(o['radial_info']),
        }
        for method, field in RECORD_COUNTS:
            same[f'{method} {field}'] = p[f"{method}_{field.replace('numvectors', 'num_vectors')}"] == o[method][field]
        for name, equal in same.items():
            tally.values(('record', name), [0 if equal else 1], [0], stamp)


def find_file(directory, pattern):
    """The one file matching pattern under directory (in any month sub-directory)."""
    files = sorted(Path(directory).rglob(pattern))
    if len(files) != 1:
        raise FileNotFoundError(f'{len(files)} files match {pattern} in {directory}')
    return files[0]


def main(prod, ours, baseline=None, hours=None):
    stamps = sorted(re.search(r'(\d{4}_\d{2}_\d{2}_\d{4})\.nc$', f.name).group(1) for f in (prod / 'oi_nc').glob('*.nc'))
    if hours:
        stamps = [s for s in stamps if s.startswith(hours)]
    if not stamps:
        sys.exit(f'No production OI NetCDF files in {prod / "oi_nc"}' + (f' for {hours}' if hours else ''))
    print(f'{len(stamps)} hours, {stamps[0]} to {stamps[-1]}')

    production = Tally()
    for stamp in stamps:
        compare_netcdf(production, 'oi', find_file(prod / 'oi_nc', f'*_oi_*_{stamp}.nc'),
                       find_file(ours / 'oi', f'*_oi_*_{stamp}.nc'), stamp, OI_FLOAT_TOLERANCE)
        compare_uwls_mat(production, find_file(prod / 'lsq_mat', f'*_uwls_*_{stamp}.mat'),
                         find_file(ours / 'lsq', f'*_uwls_*_{stamp}.nc'), stamp)
    export = prod / 'totals_configs_export.json'
    if export.exists():
        compare_records(production, export, ours, stamps)
    else:
        print(f'No {export}: the totals_generated records are not compared')
    production.report('Against production')
    ok = production.ok()

    if baseline:
        previous = Tally()
        for stamp in stamps:
            for method, short in [('lsq', 'uwls'), ('oi', 'oi')]:
                compare_netcdf(previous, short, find_file(baseline / method, f'*_{short}_*_{stamp}.nc'),
                               find_file(ours / method, f'*_{short}_*_{stamp}.nc'), stamp)
        previous.report(f'Against the baseline run {baseline}')
        ok = ok and previous.ok()

    print('\nPASS: no differences' if ok else '\nFAIL: differences found')
    return ok


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('prod', type=Path, help='Directory with the production oi_nc and lsq_mat files')
    parser.add_argument('ours', type=Path, help='realtime_totals.py output directory (with lsq and oi)')
    parser.add_argument('--baseline', type=Path, help='Earlier realtime_totals.py output to compare with exactly')
    parser.add_argument('--hours', help='Only compare hours whose YYYY_MM_DD_HHMM stamp starts with this')
    args = parser.parse_args()
    sys.exit(0 if main(args.prod.expanduser(), args.ours.expanduser(),
                       args.baseline.expanduser() if args.baseline else None, args.hours) else 1)
