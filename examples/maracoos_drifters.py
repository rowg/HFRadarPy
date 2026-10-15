"""Independent validation of MARACOOS total configurations against surface drifters.

    python examples/maracoos_drifters.py --fetch
    python examples/maracoos_drifters.py --drifters drifters.csv

Drifter positions come from the OSMC ERDDAP real-time GTS feed (platform type DRIFTING BUOYS), downloaded with
--fetch, or from a CSV with columns platform_code, time, latitude, longitude. Following Nwankwo et al. (2026), drifter
velocities are central differences of the positions on the GRS80 ellipsoid (hfradarpy.validation.drifter_velocity),
speeds above 150 cm/s are removed, and a centered 3-hour moving average matches the radial averaging.

For every hour a drifter is within 4.3 km of a grid point, each configuration's total at that grid point is computed
from all sites (only those grid points are solved) and quality controlled. Configurations are those of
maracoos_study.py and selected OI settings of maracoos_sweep.py. Scores: complex correlation and veering angle
(Kundu 1976), skill (Fang et al. 2015), RMS vector difference and mean difference, on the samples every configuration
covers and on each configuration's own samples.

Global Drifter Program drifters are usually drogued at 15 m, while 5 MHz HF radar measures the top ~2 m, so part of
any difference is vertical shear (and windage if a drogue was lost; the real-time feed does not say).
"""
import argparse
import copy
import json
import logging
import urllib.parse
import urllib.request
import warnings
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

import realtime_totals as rt
from hfradarpy.calc import lonlat2km
from hfradarpy.totals import combineRadials
from hfradarpy.validation import complex_correlation, drifter_velocity, harmonize, skill_score
from maracoos_study import MAX_DISTANCE_KM, combine, configurations as study_configurations, quality_control, xyz
from maracoos_sweep import MAX_SPEED, MDLVAR, configurations as sweep_configurations, load_inputs, setup

EXAMPLES = Path(__file__).resolve().parent
ERDDAP = 'https://osmc.noaa.gov/erddap/tabledap/OSMC_flattened.csv'
SWEEP_DEFAULT = ['sx15_e264', 'sx10_e264', 'sx15_e132', 'sx10_e66', 'sx10_etmp', 'sx15_e66_window', 'sx15_e66_window_T6']


def fetch(start, end, bbox, out):
    """Download drifter positions from the OSMC ERDDAP (GTS) into out."""
    lon0, lon1, lat0, lat1 = bbox
    query = ','.join(['platform_code', 'platform_type', 'time', 'latitude', 'longitude'])
    constraints = [f'platform_type="DRIFTING BUOYS"', f'time>={start:%Y-%m-%dT%H:%M:%SZ}',
                   f'time<={end:%Y-%m-%dT%H:%M:%SZ}', f'latitude>={lat0}', f'latitude<={lat1}',
                   f'longitude>={lon0}', f'longitude<={lon1}']
    url = ERDDAP + '?' + query + ''.join('&' + urllib.parse.quote(c, safe='=') for c in constraints)
    out.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url, timeout=300) as response:
        out.write_bytes(response.read())
    return out


def drifter_samples(path, hours, max_speed=150.0, window='3h', min_travel_km=2.0):
    """
    Hourly drifter velocities (cm/s) and positions at the HF radar hours. Platforms that move less than
    min_travel_km over the period are not free drifters (e.g. a fixed platform reported as a drifting buoy) and are
    left out.
    """
    d = pd.read_csv(path, skiprows=[1], parse_dates=['time'])
    d['time'] = d['time'].dt.tz_localize(None)
    samples = []
    for code, g in d.groupby('platform_code'):
        g = g.sort_values('time').drop_duplicates('time')
        g = g[(g['time'] >= hours[0] - pd.Timedelta('3h')) & (g['time'] <= hours[-1] + pd.Timedelta('3h'))]
        if len(g) < 3:
            continue
        dx, dy = lonlat2km(g['longitude'].iloc[0], g['latitude'].iloc[0], g['longitude'], g['latitude'])
        if np.hypot(dx, dy).max() < min_travel_km:
            print(f'{code}: moved less than {min_travel_km} km, not a free drifter; left out', flush=True)
            continue
        u, v = drifter_velocity(g['longitude'], g['latitude'], g['time'])
        span = (g['time'].shift(-1) - g['time'].shift(1)).dt.total_seconds() / 3600
        bad = ~(span <= 2.5)  # central differences only across consecutive, roughly hourly fixes
        u, v = np.where(bad, np.nan, u * 100), np.where(bad, np.nan, v * 100)
        h = harmonize(u, v, g['time'], window=window, max_speed=max_speed)
        # Positions and velocities at the HF radar hours, from fixes within 30 minutes
        track = pd.DataFrame({'lon': g['longitude'].to_numpy(), 'lat': g['latitude'].to_numpy(),
                              'u': h['u'].to_numpy(), 'v': h['v'].to_numpy()}, index=h.index)
        at = track.reindex(track.index.union(hours)).interpolate('time', limit_area='inside').reindex(hours)
        nearest = pd.Series(track.index, index=track.index).reindex(hours, method='nearest',
                                                                     tolerance=pd.Timedelta('30min'))
        at = at[nearest.notna().to_numpy()].dropna()
        at['platform_code'] = str(code)
        samples.append(at)
    s = pd.concat(samples)
    s.index.name = 'time'
    return s.reset_index()


def evaluated_configurations(conf, sweep_names):
    prod = {k: conf['Totals']['OI'][k] for k in ('mdlvar', 'errvar', 'sx', 'sy')}
    prod['weighting_function'] = conf['Totals']['OI'].get('weighting_function') or 'exponential'
    configs = {}
    for key, combo in study_configurations(prod).items():
        for name, variant in combo['outputs'].items():
            configs[name] = ('study', combo, variant)
    sweep = sweep_configurations()
    for name in sweep_names:
        configs[f'sweep {name}'] = ('sweep', sweep[name], None)
    return configs


def totals_at(t, radials, points, kind, combo, variant, grid, land, conf, driver, polygon):
    """u, v of one configuration at grid points `points` for hour t, NaN where no vector passes QC."""
    solve = grid.iloc[points].reset_index(drop=True)
    nan = np.full(len(points), np.nan)
    if kind == 'study':
        total = combine(t, radials, combo, solve, conf, driver)
        if total is None:
            return nan, nan
        d = quality_control(total, combo, variant, conf, driver, polygon).data
        ok = (d['PRIM'] == 1).to_numpy()
        return np.where(ok, d['VELU'].to_numpy(float), np.nan), np.where(ok, d['VELV'].to_numpy(float), np.nan)

    rows, sites = [], []
    for k in range(-combo['window'], combo['window'] + 1):
        for s, r in radials.get(t + pd.Timedelta(hours=k), {}).items():
            rows.append(r)
            sites.append(s)
    rdf = pd.DataFrame({'Radial': pd.Series(rows, index=sites, dtype=object)})
    tc = conf['Totals']
    if rdf.index.nunique() < tc.get('MinNumSites', 2):
        return nan, nan
    kwargs = dict(gRes=6000, tStp=t.to_pydatetime(), minContrSites=tc.get('MinNumSites', 2),
                  minContrRads=tc.get('MinNumRads', 3), useLandMask=False, dropPoints=False, engine='vectorized',
                  method='oi', mdlvar=MDLVAR, errvar=combo['errvar'], sx=combo['sx'], sy=combo['sx'],
                  oi_option='exponential', sRad=2 * combo['sx'] * 1000, oi_time_scale=combo['time_scale'],
                  oi_errvar_column=combo['errcol'])
    if combo['window']:
        kwargs['tempthreshold'] = pd.Timedelta(hours=combo['window'])
    d = combineRadials(rdf, solve, **kwargs)[0].data
    u, v = d['VELU'].to_numpy(float), d['VELV'].to_numpy(float)
    ok = (np.isfinite(u) & (np.hypot(u, v) <= MAX_SPEED) & ~land[points]
          & (np.maximum(d['UQAL'].to_numpy(float), d['VQAL'].to_numpy(float)) <= 0.6))
    return np.where(ok, u, np.nan), np.where(ok, v, np.nan)


def scores(ref_u, ref_v, u, v):
    ok = np.isfinite(ref_u) & np.isfinite(ref_v) & np.isfinite(u) & np.isfinite(v)
    if ok.sum() < 3:
        return dict(samples=int(ok.sum()))
    r, angle = complex_correlation(ref_u[ok], ref_v[ok], u[ok], v[ok])
    skill, _, _ = skill_score(ref_u[ok], ref_v[ok], u[ok], v[ok])
    du, dv = u[ok] - ref_u[ok], v[ok] - ref_v[ok]
    return dict(samples=int(ok.sum()), complex_correlation=r, veering_deg=angle, skill=skill,
                rms_difference_cm_s=np.sqrt(np.mean(du ** 2 + dv ** 2)), mean_du_cm_s=du.mean(), mean_dv_cm_s=dv.mean())


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--data', type=Path, default=Path('~/Documents/hfradar_parity/2026-10-06_week'))
    parser.add_argument('--start', default='2026-10-06T00:00')
    parser.add_argument('--end', default='2026-10-12T23:00')
    parser.add_argument('--drifters', type=Path, help='CSV of drifter positions (platform_code, time, latitude, longitude)')
    parser.add_argument('--fetch', action='store_true', help='Download drifter positions from the OSMC ERDDAP')
    parser.add_argument('--sweep-configs', nargs='*', default=SWEEP_DEFAULT)
    parser.add_argument('--output', type=Path, default=EXAMPLES / 'output/maracoos_drifters')
    args = parser.parse_args()
    warnings.simplefilter('ignore', category=FutureWarning)
    logging.disable(logging.WARNING)

    out = args.output
    for sub in ('tables', 'figures'):
        (out / sub).mkdir(parents=True, exist_ok=True)
    hours = pd.date_range(args.start, args.end, freq='h')
    conf, driver, grid, land = setup(args.data)
    glon, glat = grid.x.to_numpy(), grid.y.to_numpy()
    path = args.drifters
    if args.fetch:
        bbox = (glon.min() - 0.5, glon.max() + 0.5, glat.min() - 0.5, glat.max() + 0.5)
        path = fetch(hours[0] - pd.Timedelta('1D'), hours[-1] + pd.Timedelta('1D'), bbox, out / 'osmc_drifters.csv')
    if path is None:
        parser.error('use --drifters or --fetch')

    samples = drifter_samples(path, hours)
    dist, point = cKDTree(xyz(glon, glat)).query(xyz(samples['lon'].to_numpy(), samples['lat'].to_numpy()))
    samples = samples[dist <= MAX_DISTANCE_KM].assign(point=point[dist <= MAX_DISTANCE_KM])
    print(f'{len(samples)} drifter hours within {MAX_DISTANCE_KM} km of the grid, '
          f'{samples.platform_code.nunique()} drifters', flush=True)

    inputs = load_inputs(conf, hours)
    polygon = rt.load_polygon(conf['QC'].get('QC305_MaskFile'))
    configs = evaluated_configurations(conf, args.sweep_configs)
    for name, (kind, combo, variant) in configs.items():
        u = np.full(len(samples), np.nan)
        v = np.full(len(samples), np.nan)
        for t, idx in samples.groupby('time').indices.items():
            radials = {t + pd.Timedelta(hours=k): inputs.get(t + pd.Timedelta(hours=k), {}) for k in (-1, 0, 1)}
            u[idx], v[idx] = totals_at(t, radials, samples['point'].to_numpy()[idx], kind, combo, variant, grid, land,
                                       conf, driver, polygon)
        samples[f'{name}::u'], samples[f'{name}::v'] = u, v
        print(f'{name}: {int(np.isfinite(u).sum())} samples', flush=True)
    samples.to_csv(out / 'tables' / 'drifter_samples.csv', index=False)

    names = list(configs)
    common = np.logical_and.reduce([np.isfinite(samples[f'{n}::u']).to_numpy() for n in names])
    rows = []
    for name in names:
        u, v = samples[f'{name}::u'].to_numpy(), samples[f'{name}::v'].to_numpy()
        own = scores(samples['u'].to_numpy(), samples['v'].to_numpy(), u, v)
        shared = scores(samples['u'].to_numpy()[common], samples['v'].to_numpy()[common], u[common], v[common])
        rows.append({'config': name, **{f'{k}': val for k, val in own.items()},
                     **{f'common_{k}': val for k, val in shared.items()}})
    table = pd.DataFrame(rows).set_index('config')
    table.to_csv(out / 'tables' / 'drifter_scores.csv')
    per_drifter = []
    for code, g in samples.groupby('platform_code'):
        for name in names:
            per_drifter.append({'platform_code': code, 'config': name,
                                **scores(g['u'].to_numpy(), g['v'].to_numpy(), g[f'{name}::u'].to_numpy(),
                                         g[f'{name}::v'].to_numpy())})
    pd.DataFrame(per_drifter).to_csv(out / 'tables' / 'drifter_scores_by_drifter.csv', index=False)

    # Time series of each drifter with production OI and production UWLS
    codes = sorted(samples.platform_code.unique())
    fig, axes = plt.subplots(2 * len(codes), 1, figsize=(13, 3.2 * 2 * len(codes)), sharex=True, constrained_layout=True)
    for k, code in enumerate(codes):
        g = samples[samples.platform_code == code]
        for j, comp in enumerate(('u', 'v')):
            ax = axes[2 * k + j]
            ax.plot(g['time'], g[comp], color='black', lw=1.6, label='drifter (3 h mean)')
            ax.plot(g['time'], g[f'OI production::{comp}'], color='#eb6834', lw=1.2, label='OI production')
            ax.plot(g['time'], g[f'UWLS production::{comp}'], color='#2a78d6', lw=1.0, alpha=0.8,
                    label='UWLS production')
            ax.set_ylabel(f'{comp} (cm/s)')
            ax.set_title(f'Drifter {code}, {"eastward" if comp == "u" else "northward"} velocity', fontsize=9)
            ax.grid(alpha=0.3)
    axes[0].legend(fontsize=8, frameon=False, ncol=3)
    fig.savefig(out / 'figures' / 'drifter_timeseries.png', dpi=120)
    plt.close(fig)

    summary = dict(drifters=codes, samples=int(len(samples)), common_samples=int(common.sum()),
                   source=str(path), note='GTS drifters, likely drogued at 15 m; drogue status unknown')
    (out / 'drifters.json').write_text(json.dumps(summary, indent=2))
    with pd.option_context('display.width', 250, 'display.max_columns', 30, 'display.float_format', '{:.3f}'.format):
        print(table)


if __name__ == '__main__':
    main()
