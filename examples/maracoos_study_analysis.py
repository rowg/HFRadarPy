"""Analysis and artifacts for maracoos_study.py; no numerical tuning on the evaluation week."""
import json
from pathlib import Path
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from hfradarpy.gapfill import fill_gaps, fill_gauss_markov
from maracoos_study import BASELINE, SEED, VARIABLES, atomic_json, digest


def load_totals(ctx):
    names = [n for c in ctx['combos'].values() for n in c['outputs']]
    results = {n: {v: [] for v in VARIABLES} for n in names}
    for key, combo in ctx['combos'].items():
        for t in ctx['hours']:
            path = ctx['run'] / 'totals' / f'{key}_{t:%Y%m%dT%H%M}.npz'
            with np.load(path) as arrays:
                for name in combo['outputs']:
                    for v in VARIABLES:
                        results[name][v].append(arrays[f'{name}::{v}'].astype(np.float32))
    return {n: {v: np.stack(a) for v, a in d.items()} for n, d in results.items()}


def load_loso(ctx):
    manifest = json.loads((ctx['run'] / 'manifest.json').read_text())
    step = manifest['loso_step_hours']
    counts = pd.read_csv(ctx['run'] / 'tables/input_counts.csv', parse_dates=['time'])
    expected_cases = {(t, s) for t, s in zip(counts.time, counts.site) if t in ctx['hours'][::step]}
    pieces = []
    for key, combo in ctx['combos'].items():
        paths = sorted((ctx['run'] / 'loso').glob(f'{key}_[0-9]*.npz'))
        if not paths:
            raise ValueError(f'No leave-one-site-out results for {key}')
        cases = set()
        for path in paths:
            stamp, site = path.stem[len(key) + 1:].split('_', 1)
            cases.add((pd.Timestamp(stamp), site))
            with np.load(path) as a:
                index = pd.MultiIndex.from_arrays(
                    [[pd.Timestamp(stamp)] * len(a['observed']), [site] * len(a['observed']), np.arange(len(a['observed']))],
                    names=['time', 'site', 'bin'])
                d = pd.DataFrame({n: a[n] for n in combo['outputs']}, index=index)
                d['observed'] = a['observed']
                d['eligible'] = a['eligible']
                pieces.append(d)
        if cases != expected_cases:
            raise ValueError(f'Incomplete LOSO site/time coverage for {key}')
    # Grouping joins the algorithms on explicit observation identities, preserving
    # rows for which every configuration failed. Never align a partial boolean index.
    joined = pd.concat(pieces).groupby(level=['time', 'site', 'bin'], sort=True).first()
    names = [n for c in ctx['combos'].values() for n in c['outputs']]
    if set(joined.index.get_level_values('time')) != set(ctx['hours'][::step]):
        raise ValueError('Incomplete LOSO time coverage')
    return joined.reindex(columns=['observed', 'eligible', *names])


def valid(r):
    return (r['PRIM'] == 1) & np.isfinite(r['VELU']) & np.isfinite(r['VELV'])


def coverage_table(ctx, results):
    rows, hourly = [], []
    for name, r in results.items():
        good = valid(r)
        computed = np.isfinite(r['VELU']) & np.isfinite(r['VELV'])
        rows.append(dict(config=name, computed_per_hour=computed.sum()/len(good),
                         pass_per_hour=good.sum()/len(good), qc_removed_pct=100*(1-good.sum()/max(1,computed.sum())),
                         ever_covered=int(good.any(axis=0).sum()),
                         covered_half_the_time=int((good.mean(axis=0) >= .5).sum())))
        hourly.append(pd.DataFrame({'time': ctx['hours'], 'config': name, 'passing': good.sum(axis=1)}))
    return pd.DataFrame(rows).set_index('config'), pd.concat(hourly, ignore_index=True)


def paired_interval(errors, reference, times, block_hours=24, draws=2000, seed=SEED):
    """Paired pooled-RMSE difference, resampling whole time blocks across all sites.

    Inputs must already be paired on identical observations. Neighboring bins are
    never treated as independent bootstrap draws. Returns candidate minus reference.
    """
    errors, reference = np.asarray(errors), np.asarray(reference)
    if len(errors) == 0:
        return np.nan, np.nan, 0
    times = pd.DatetimeIndex(times)
    blocks = ((times - times.min().normalize()) / pd.Timedelta(hours=block_hours)).astype(int)
    sums = pd.DataFrame({'e2': errors**2, 'r2': reference**2, 'n': 1, 'block': blocks}).groupby('block').sum()
    if len(sums) < 2:
        return np.nan, np.nan, len(sums)
    rng = np.random.default_rng(seed)
    sample = rng.integers(0, len(sums), size=(draws, len(sums)))
    totals = sums.to_numpy()[sample].sum(axis=1)
    delta = np.sqrt(totals[:,0]/totals[:,2]) - np.sqrt(totals[:,1]/totals[:,2])
    low, high = np.quantile(delta, [.025, .975])
    return low, high, len(sums)


def skill_tables(loso, names):
    eligible = loso.eligible.astype(bool) & np.isfinite(loso.observed)
    shared = eligible & np.isfinite(loso[names]).all(axis=1)
    rows, groups = [], []
    times = loso.index.get_level_values('time')
    for name in names:
        ok = eligible & np.isfinite(loso[name])
        pair = ok & np.isfinite(loso[BASELINE])
        error = loso[name] - loso.observed
        ref = loso[BASELINE] - loso.observed
        n = int(ok.sum())
        rmse = lambda a: float(np.sqrt(np.mean(np.asarray(a)**2))) if len(a) else np.nan
        row = dict(config=name, eligible_bins=int(eligible.sum()), predicted_bins=n,
                   coverage_pct=100*n/max(1,eligible.sum()), rmse_cm_s=rmse(error[ok]),
                   bias_cm_s=error[ok].mean(), median_abs_error_cm_s=error[ok].abs().median(),
                   all_common_bins=int(shared.sum()), all_common_rmse_cm_s=rmse(error[shared]),
                   paired_bins=int(pair.sum()), paired_rmse_cm_s=rmse(error[pair]),
                   baseline_paired_rmse_cm_s=rmse(ref[pair]),
                   delta_rmse_cm_s=rmse(error[pair])-rmse(ref[pair]))
        for block in (24, 48):
            low, high, count = paired_interval(error[pair], ref[pair], times[pair], block)
            row.update({f'delta_ci_low_{block}h': low, f'delta_ci_high_{block}h': high, f'blocks_{block}h': count})
        rows.append(row)
        paired = pd.DataFrame({'error2': error[pair]**2, 'reference2': ref[pair]**2,
                               'error': error[pair], 'day': times[pair].normalize()})
        for group in ('site', 'day'):
            gb = paired.groupby(level='site') if group == 'site' else paired.groupby('day')
            for label, g in gb:
                groups.append(dict(config=name, group=group, label=str(label), bins=len(g),
                                   rmse_cm_s=np.sqrt(g.error2.mean()), baseline_rmse_cm_s=np.sqrt(g.reference2.mean()),
                                   delta_rmse_cm_s=np.sqrt(g.error2.mean())-np.sqrt(g.reference2.mean()),
                                   bias_cm_s=g.error.mean()))
    return pd.DataFrame(rows).set_index('config'), pd.DataFrame(groups)


def ablation_table(loso, comparisons):
    rows = []
    eligible = loso.eligible.astype(bool) & np.isfinite(loso.observed)
    for label, reference, candidate in comparisons:
        rp = eligible & np.isfinite(loso[reference])
        cp = eligible & np.isfinite(loso[candidate])
        both = rp & cp
        re = loso.loc[both, reference] - loso.loc[both, 'observed']
        ce = loso.loc[both, candidate] - loso.loc[both, 'observed']
        low, high, _ = paired_interval(ce, re, loso.index.get_level_values('time')[both])
        rows.append(dict(change=label, reference=reference, candidate=candidate, paired_bins=int(both.sum()),
                         delta_rmse_cm_s=np.sqrt(np.mean(ce**2))-np.sqrt(np.mean(re**2)),
                         ci_low=low, ci_high=high,
                         delta_coverage_pp=100*(cp.sum()-rp.sum())/max(1,eligible.sum())))
    return pd.DataFrame(rows)


def daily_dataset(ctx, r):
    grid = ctx['grid']
    lon, lat = np.unique(grid.x.to_numpy()), np.unique(grid.y.to_numpy())
    ix, iy = np.searchsorted(lon, grid.x), np.searchsorted(lat, grid.y)
    good = valid(r)
    fields = {}
    for v in ('VELU', 'VELV'):
        hourly = pd.DataFrame(np.where(good, r[v], np.nan), index=ctx['hours'])
        daily = hourly.resample('1D').mean().where(hourly.notna().resample('1D').sum() >= 12)
        values = np.full((len(daily), len(lat), len(lon)), np.nan)
        values[:, iy, ix] = daily.to_numpy()
        fields[v] = (('time', 'lat', 'lon'), values)
    return xr.Dataset(fields, coords={'time': daily.index, 'lat': lat, 'lon': lon})


def holdout_mask(observed, kind, seed):
    """Fixed masks shared across configurations; only known common values are scored."""
    rng = np.random.default_rng(seed)
    nt, ny, nx = observed.shape
    hide = np.zeros(observed.shape, bool)
    if kind == 'random_10pct':
        hide = rng.random(observed.shape) < .1
    elif kind == 'spatial_patches':
        # 5x5 grid-cell patches, one center in observed water per day (~30 km wide).
        for t in range(nt):
            candidates = np.argwhere(observed[t])
            if len(candidates):
                for y, x in candidates[rng.integers(0, len(candidates), size=5)]:
                    hide[t, max(0,y-2):y+3, max(0,x-2):x+3] = True
    elif kind == 'two_day_outage':
        # Interior outage leaves observations on both sides for the temporal baseline.
        if nt >= 4:
            t = int(rng.integers(1, nt-2))
            hide[t:t+2] = True
    elif kind == 'space_time_patch':
        if nt >= 4:
            t = int(rng.integers(1, nt-2))
            candidates = np.argwhere(observed[t])
            if len(candidates):
                for y, x in candidates[rng.integers(0, len(candidates), size=10)]:
                    hide[t:t+2, max(0,y-2):y+3, max(0,x-2):x+3] = True
    else:
        raise ValueError(kind)
    return hide & observed


def linear_time_fill(ds):
    """Simple baseline: linear interpolation in time, between observations only."""
    return ds.interpolate_na(dim='time', method='linear', max_gap=pd.Timedelta('3D'))


def gap_experiments(ctx, results):
    names = ['UWLS production', BASELINE, 'Nwankwo (full)']
    daily = {n: daily_dataset(ctx, results[n]) for n in names}
    common = np.logical_and.reduce([np.isfinite(d.VELU.values) & np.isfinite(d.VELV.values) for d in daily.values()])
    if common.shape[0] < 4:
        raise ValueError('Gap experiments need at least four days')
    domain = np.zeros(common.shape[1:], bool)
    grid = ctx['grid']
    ds0 = next(iter(daily.values()))
    ix, iy = np.searchsorted(ds0.lon, grid.x), np.searchsorted(ds0.lat, grid.y)
    from hfradarpy.calc import inpolygon
    water = np.ones(len(grid), bool)
    if ctx['polygon'] is not None:
        inside = inpolygon(grid.x.to_numpy(), grid.y.to_numpy(), *ctx['polygon'])
        water = ~inside if ctx['driver']['qc305_polygon'] == 'land' else inside
    domain[iy[water],ix[water]] = True

    def mask(filled):
        for v in ('VELU', 'VELV'):
            filled[v] = filled[v].where(domain)
        if 'fill_flag' in filled:
            filled['fill_flag'] = filled.fill_flag.where(domain, 9).astype(np.int8)
        return filled

    rows = []
    for name, ds in daily.items():
        for kind in ('random_10pct', 'spatial_patches', 'two_day_outage', 'space_time_patch'):
            for repeat in range(5):
                hide = holdout_mask(common, kind, SEED+repeat)
                test = ds.copy(deep=True)
                for v in ('VELU', 'VELV'):
                    test[v] = test[v].where(~hide)
                nearest = mask(fill_gaps(test))
                estimates = {'linear_time': mask(linear_time_fill(test)),
                             'nearest_neighbor': nearest,
                             'nearest_plus_gauss_markov': mask(fill_gauss_markov(nearest))}
                # Bound comparisons to the same hidden points recovered by both methods.
                for method, filled in estimates.items():
                    ok = hide & np.isfinite(filled.VELU.values) & np.isfinite(filled.VELV.values)
                    error2 = (filled.VELU.values-ds.VELU.values)**2 + (filled.VELV.values-ds.VELV.values)**2
                    for reference in ('linear_time', 'nearest_neighbor'):
                        ref = estimates[reference]
                        paired = ok & np.isfinite(ref.VELU.values) & np.isfinite(ref.VELV.values)
                        ref2 = (ref.VELU.values-ds.VELU.values)**2 + (ref.VELV.values-ds.VELV.values)**2
                        rows.append(dict(config=name, gap=kind, repeat=repeat, method=method, reference=reference,
                            hidden=int(hide.sum()), filled=int(ok.sum()), recovered_pct=100*ok.sum()/max(1,hide.sum()),
                            squared_error_sum=float(np.nansum(error2[ok])),
                            paired=int(paired.sum()), paired_squared_error_sum=float(np.nansum(error2[paired])),
                            reference_squared_error_sum=float(np.nansum(ref2[paired]))))
            print(f'Gap tests: {name}, {kind}', flush=True)
    detail = pd.DataFrame(rows)
    keys = ['config', 'gap', 'method', 'reference']
    summary = detail.groupby(keys)[['hidden','filled','squared_error_sum','paired','paired_squared_error_sum',
                                    'reference_squared_error_sum']].sum()
    summary['recovered_pct'] = 100*summary.filled/summary.hidden
    summary['rmse_cm_s'] = np.sqrt(summary.squared_error_sum/summary.filled.replace(0,np.nan))
    summary['paired_rmse_cm_s'] = np.sqrt(summary.paired_squared_error_sum/summary.paired.replace(0,np.nan))
    summary['reference_rmse_cm_s'] = np.sqrt(summary.reference_squared_error_sum/summary.paired.replace(0,np.nan))
    summary['delta_rmse_cm_s'] = summary.paired_rmse_cm_s-summary.reference_rmse_cm_s
    # Export observed and filled fields with a domain mask, keeping land/outside-grid blank.
    for name, ds in daily.items():
        filled = mask(fill_gauss_markov(mask(fill_gaps(ds))))
        for v in ('VELU','VELV'):
            filled[v] = filled[v].where(domain)
            filled[v].attrs['units'] = 'cm s-1'
        filled['fill_flag'] = filled.fill_flag.where(domain, 9).astype(np.int8)
        filled.attrs.update(configuration=name, source_fingerprint=ctx['run'].name,
                            comment='Daily means require 12 hourly QC-passing vectors. Gap fills are estimates.')
        filled.to_netcdf(ctx['run'] / f'daily_{name.replace(" ","_")}.nc')
    return detail, summary.reset_index()


def parity(ctx, results):
    """Compare this study's production OI fields with archived published OI NetCDF."""
    settings = json.loads((ctx['run'] / 'manifest.json').read_text())['settings']
    export = next(Path(p) for p in settings if p.endswith('totals_configs_export.json'))
    prod = export.parent / 'oi_nc'
    rows = []
    grid = ctx['grid']
    for i, t in enumerate(ctx['hours']):
        files = list(prod.glob(f'*{t:%Y_%m_%d_%H%M}.nc'))
        if len(files) != 1:
            raise ValueError(f'Expected one production OI file for {t}, got {len(files)}')
        with xr.open_dataset(files[0], mask_and_scale=False) as ds:
            sample = ds.sel(lon=xr.DataArray(grid.x.to_numpy(), dims='point'),
                            lat=xr.DataArray(grid.y.to_numpy(), dims='point'), method='nearest').squeeze(drop=True)
            for source, target in [('VELU','u'), ('VELV','v'), ('PRIM','qc_primary_flag')]:
                theirs = sample[target].values.astype(float)
                fill = sample[target].attrs.get('_FillValue')
                if fill is not None:
                    theirs[theirs == fill] = np.nan
                ours = results[BASELINE][source][i].astype(float)
                if source != 'PRIM':
                    ours = np.round(ours.astype(np.float32)).astype(float)
                both = np.isfinite(ours) & np.isfinite(theirs)
                rows.append(dict(time=t, variable=target, paired=int(both.sum()),
                                 missing_mismatches=int((np.isfinite(ours)!=np.isfinite(theirs)).sum()),
                                 different=int((ours[both]!=theirs[both]).sum()),
                                 max_abs_difference=float(np.max(np.abs(ours[both]-theirs[both]))) if both.any() else 0))
    return pd.DataFrame(rows)


def make_figures(ctx, results, coverage, skill, groups, hourly, gaps):
    figdir = ctx['run'] / 'figures'
    names = list(skill.index)
    colors = plt.get_cmap('tab10')(np.arange(len(names)))
    fig, axes = plt.subplots(1, 2, figsize=(14,5), layout='constrained')
    y = np.arange(len(names))
    axes[0].barh(y, skill.all_common_rmse_cm_s, color=colors)
    axes[1].barh(y, skill.coverage_pct, color=colors)
    for ax in axes:
        ax.set_yticks(y, names, fontsize=8)
        ax.invert_yaxis()
        ax.grid(axis='x', alpha=.25)
    axes[0].set_xlabel('Withheld radial RMSE, all-common bins (cm/s)')
    axes[1].set_xlabel('Eligible withheld bins predicted (%)')
    fig.savefig(figdir/'accuracy_coverage.png', dpi=160)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(10,5), layout='constrained')
    for i,n in enumerate(names):
        row = skill.loc[n]
        ax.plot([row.delta_ci_low_24h,row.delta_ci_high_24h],[i,i], color=colors[i], lw=2)
        ax.plot(row.delta_rmse_cm_s,i,'o',color=colors[i])
    ax.axvline(0, color='black', lw=.8)
    ax.set_yticks(y,names,fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel('Paired RMSE change versus production OI (cm/s); negative is lower error')
    ax.set_title('Exploratory 95% bootstrap intervals; seven 24-hour blocks')
    fig.savefig(figdir/'paired_improvement.png',dpi=160)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(13,4), layout='constrained')
    for i,n in enumerate(names):
        d=hourly[hourly.config==n]
        ax.plot(pd.to_datetime(d.time),d.passing,label=n,color=colors[i])
    ax.legend(fontsize=7,ncol=2)
    ax.set_ylabel('Grid points passing QC')
    fig.savefig(figdir/'hourly_coverage.png',dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(2,4,figsize=(15,8),layout='constrained')
    grid=ctx['grid']
    for ax,n in zip(axes.flat,names):
        pct=100*valid(results[n]).mean(axis=0)
        good=pct>0
        sc=ax.scatter(grid.x[good],grid.y[good],c=pct[good],s=2,cmap='viridis',vmin=0,vmax=100)
        ax.set_title(n,fontsize=8)
        ax.set_aspect(1/np.cos(np.deg2rad(39)))
    fig.colorbar(sc,ax=axes,label='Hours with a QC-passing vector (%)',shrink=.7)
    fig.savefig(figdir/'coverage_maps.png',dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(1,2,figsize=(12,4),layout='constrained')
    for i,n in enumerate(names):
        r=results[n]
        computed=np.isfinite(r['VELU']) & np.isfinite(r['VELV'])
        if n != 'UWLS production':
            axes[0].hist(np.maximum(r['UQAL'],r['VQAL'])[computed],bins=np.linspace(0,1,41),
                         density=True,histtype='step',color=colors[i],label=n)
        axes[1].hist(r['MSFR'][computed],bins=np.linspace(0,1,41),density=True,histtype='step',color=colors[i],label=n)
    axes[0].set_xlabel('Normalized uncertainty, max(Uerr, Verr)')
    axes[1].set_xlabel('Largest contributing site share')
    for ax in axes:
        ax.set_ylabel('Density')
    axes[1].legend(fontsize=6)
    fig.savefig(figdir/'uncertainty_site_share.png',dpi=160)
    plt.close(fig)
    for group in ('site','day'):
        d=groups[groups.group==group].pivot(index='label',columns='config',values='delta_rmse_cm_s').reindex(columns=names)
        fig,ax=plt.subplots(figsize=(12,5),layout='constrained')
        im=ax.imshow(d.to_numpy().T,aspect='auto',cmap='RdBu_r',vmin=-5,vmax=5)
        ax.set_yticks(y,names,fontsize=8)
        ax.set_xticks(np.arange(len(d)),[s[:10] if group=='day' else s for s in d.index],rotation=45,ha='right')
        fig.colorbar(im,ax=ax,label='Paired RMSE change vs production OI (cm/s)')
        fig.savefig(figdir/f'error_by_{group}.png',dpi=160)
        plt.close(fig)


def independent_validation(ctx, results):
    """Optional vector comparisons on exactly the same in situ samples, with explicit units."""
    path = ctx.get('insitu_file')
    if path is None:
        atomic_json(ctx['run']/'independent_validation.json',
                    {'status': 'not_run', 'reason': 'No independent buoy/drifter observations supplied'})
        return
    from hfradarpy.validation import complex_correlation, skill_score
    from maracoos_study import xyz, MAX_DISTANCE_KM
    from scipy.spatial import cKDTree
    observations = pd.read_csv(path, parse_dates=['time'])
    required = ['time','lon','lat','u','v']
    if not all(c in observations for c in required):
        raise ValueError(f'In situ CSV needs columns {required}')
    observations = observations.dropna(subset=required).copy()
    observations['time'] = pd.to_datetime(observations.time, utc=True).dt.tz_localize(None)
    if ctx['insitu_units'] == 'm/s':
        observations[['u','v']] *= 100
    grid = ctx['grid']
    dist, point = cKDTree(xyz(grid.x.to_numpy(),grid.y.to_numpy())).query(
        xyz(observations.lon.to_numpy(),observations.lat.to_numpy()))
    ti = ctx['hours'].get_indexer(pd.DatetimeIndex(observations.time), method='nearest', tolerance=pd.Timedelta('30min'))
    eligible = (dist <= MAX_DISTANCE_KM) & (ti >= 0)
    predictions = {}
    for name,r in results.items():
        u,v = r['VELU'][np.maximum(ti,0),point].copy(), r['VELV'][np.maximum(ti,0),point].copy()
        ok = eligible & valid(r)[np.maximum(ti,0),point]
        u[~ok],v[~ok] = np.nan,np.nan
        predictions[name]=(u,v)
    common = np.logical_and.reduce([np.isfinite(u)&np.isfinite(v) for u,v in predictions.values()])
    rows=[]
    for name,(u,v) in predictions.items():
        ru,rv = observations.u.to_numpy()[common],observations.v.to_numpy()[common]
        hu,hv = u[common],v[common]
        corr,angle = complex_correlation(ru,rv,hu,hv)
        score,_,_ = skill_score(ru,rv,hu,hv)
        rows.append(dict(config=name, common_samples=int(common.sum()),
                         predicted_samples=int((np.isfinite(u)&np.isfinite(v)).sum()),
                         vector_rmse_cm_s=np.sqrt(np.mean((hu-ru)**2+(hv-rv)**2)),
                         u_bias_cm_s=np.mean(hu-ru),v_bias_cm_s=np.mean(hv-rv),
                         complex_correlation=corr,veering_degrees=angle,skill=score))
    pd.DataFrame(rows).to_csv(ctx['run']/'tables/independent_validation.csv',index=False)
    atomic_json(ctx['run']/'independent_validation.json',dict(status='complete',path=str(path),
        sha256=digest(path), input_units=ctx['insitu_units'], samples=len(observations),common_samples=int(common.sum()),
        note='No filtering or depth adjustment applied; input observations must be quality controlled and comparable.'))


def markdown_table(df, columns=None):
    if columns is not None:
        df=df[columns]
    d=df.reset_index()
    def fmt(v):
        return f'{v:.3f}' if isinstance(v,(float,np.floating)) else str(v)
    return '\n'.join(['| '+' | '.join(map(str,d.columns))+' |',
                      '| '+' | '.join(['---']*len(d.columns))+' |']+
                     ['| '+' | '.join(fmt(v) for v in row)+' |' for row in d.itertuples(index=False,name=None)])


def analyze(ctx):
    run=ctx['run']
    if not (run/'compute_complete.json').exists():
        raise ValueError('Computation is incomplete; resume --stage compute first')
    results=load_totals(ctx)
    loso=load_loso(ctx)
    loso.to_parquet(run/'withheld_predictions.parquet')
    coverage,hourly=coverage_table(ctx,results)
    skill,groups=skill_tables(loso,list(results))
    counts=pd.read_csv(run/'tables/input_counts.csv',parse_dates=['time']).groupby('time').size()
    availability=pd.DataFrame({'time':ctx['hours']})
    for label,offset in [('previous',-1),('current',0),('next',1)]:
        availability[label+'_sites']=counts.reindex(ctx['hours']+pd.Timedelta(hours=offset),fill_value=0).to_numpy()
    availability['complete_window']=(availability[['previous_sites','current_sites','next_sites']]>0).all(axis=1)
    availability.to_csv(run/'tables/window_availability.csv',index=False)
    interior=loso.index.get_level_values('time').isin(availability.loc[availability.complete_window,'time'])
    interior_skill,_=skill_tables(loso.loc[interior],list(results))
    interior_skill.to_csv(run/'tables/skill_complete_windows.csv')
    agreement=[]
    for name,r in results.items():
        for ref in (BASELINE,'UWLS production'):
            other=results[ref]
            both=valid(r)&valid(other)
            error2=(r['VELU'][both]-other['VELU'][both])**2+(r['VELV'][both]-other['VELV'][both])**2
            agreement.append(dict(config=name,reference=ref,paired_vectors=int(both.sum()),
                                  rms_vector_difference_cm_s=float(np.sqrt(np.mean(error2)))))
    pd.DataFrame(agreement).to_csv(run/'tables/production_agreement.csv',index=False)
    ablations=ablation_table(loso,[
        ('Add time window',BASELINE,'OI production + window'),
        ('Change OI parameters',BASELINE,'OI Nwankwo parameters'),
        ('Add window to revised parameters','OI Nwankwo parameters','OI Nwankwo parameters + window'),
        ('Relax uncertainty 0.6 to 0.7','OI Nwankwo parameters + window','OI Nwankwo + uncertainty QC'),
        ('Add site-share limit','OI Nwankwo parameters + window','OI Nwankwo + site QC'),
        ('Add share limit after relaxed uncertainty','OI Nwankwo + uncertainty QC','Nwankwo (full)'),
    ])
    for name,table in [('coverage',coverage),('hourly_coverage',hourly),('skill',skill),('skill_by_group',groups),
                       ('ablations',ablations)]:
        table.to_csv(run/'tables'/f'{name}.csv',index=name in ('coverage','skill'))
    print(skill[['coverage_pct','paired_rmse_cm_s','delta_rmse_cm_s','delta_ci_low_24h','delta_ci_high_24h']].to_string(),flush=True)
    parity_rows=parity(ctx,results)
    parity_rows.to_csv(run/'tables/production_parity.csv',index=False)
    gap_detail,gaps=gap_experiments(ctx,results)
    gap_detail.to_csv(run/'tables/gapfill_repeats.csv',index=False)
    gaps.to_csv(run/'tables/gapfill_summary.csv',index=False)
    make_figures(ctx,results,coverage,skill,groups,hourly,gaps)
    independent_validation(ctx,results)
    write_report(ctx,skill,coverage,ablations,parity_rows,gaps)
    atomic_json(run/'analysis_complete.json',dict(analysis_sha256=digest(__file__),
                input_fingerprint=run.name, artifacts=sorted(str(p.relative_to(run)) for p in (run/'tables').glob('*.csv'))))
    print(f'Completed study: {run / "report.md"}',flush=True)


def write_report(ctx,skill,coverage,ablations,parity_rows,gaps):
    run=ctx['run']
    mismatch=int(parity_rows.different.sum()+parity_rows.missing_mismatches.sum())
    full=skill.loc['Nwankwo (full)']
    baseline=skill.loc[BASELINE]
    independent=json.loads((run/'independent_validation.json').read_text())
    independent_status=('completed; see tables/independent_validation.csv' if independent['status']=='complete'
                        else 'pending: data not supplied')
    full_gaps=gaps[(gaps.config=='Nwankwo (full)') & (gaps.reference=='linear_time')].set_index(['gap','method'])
    text=f'''# MARACOOS processing comparison

Evaluation: {ctx['hours'][0]} through {ctx['hours'][-1]}. Fingerprint: `{run.name}`.

This is a fixed-configuration, one-week exploratory evaluation. No parameters were tuned on these observations.
Independent buoy/drifter validation: **{independent_status}**. An untouched additional week is **pending: data not supplied**.

## Measured result for the full candidate

- Withheld-radial coverage: **{full.coverage_pct:.2f}%**, versus **{baseline.coverage_pct:.2f}%** for production OI.
- On the same **{int(full.paired_bins):,}** withheld observations, RMSE is **{full.paired_rmse_cm_s:.3f} cm/s**, versus
  **{full.baseline_paired_rmse_cm_s:.3f} cm/s** for production OI.
- Candidate minus production RMSE: **{full.delta_rmse_cm_s:+.3f} cm/s**; exploratory 24-hour-block 95% interval
  **[{full.delta_ci_low_24h:+.3f}, {full.delta_ci_high_24h:+.3f}]**. The 48-hour-block interval is
  **[{full.delta_ci_low_48h:+.3f}, {full.delta_ci_high_48h:+.3f}] cm/s**.

These numbers describe this week and this withheld-radar test. A negative difference means lower error, a positive
difference means higher error; the interval and coverage tradeoff must be considered together.

The parameter candidates follow the existing study notebook, informed by
[Nwankwo et al. (2026)](https://doi.org/10.1007/s10236-026-01813-8).
Values 625 and 100 are used in the native squared-cm/s units of these radials, as in the supplied notebook;
the publication labels its variances m²/s², so this is an explicit study assumption, not a verified unit correction.
The gap-filling scales come from [Nwankwo et al. (2023), Appendix A](https://doi.org/10.1029/2022JC019055).

## Production baseline check

Compared this run's production OI u, v and primary QC with the published NetCDF at the study grid,
including matching missingness. Velocities are compared after the production float32/whole-cm/s rounding.
Total mismatched entries: **{mismatch}**. See `tables/production_parity.csv` for all hours.

## Withheld-site prediction

Each available site is excluded from **every hour** of the input window. Its cleaned radial observations are predicted
from the remaining sites at the nearest grid point within 4.3 km. The interval between tested hours is six hours by
default (see manifest). A whole site/time/bin identity is used to align every comparison. These targets are other
radar observations, not independent truth; shared errors and representativeness differences can remain.

All-common RMSE uses observations predicted by all eight configurations. Paired RMSE and its difference use each
candidate's overlap with production OI. Coverage uses all finite withheld observations within the grid-distance limit.
Intervals resample whole 24-hour blocks across all sites; 48-hour sensitivity intervals are in the CSV. With only seven
days, both are exploratory, not evidence of general performance across seasons or conditions.
The intervals are not adjusted for the multiple prespecified comparisons.
The archive has no radials outside this week: the first and last full-network hours have shortened windows.
`tables/window_availability.csv` records the available sites in every neighboring hour;
`tables/skill_complete_windows.csv` repeats the skill comparison only where all three input hours are available.

{markdown_table(skill,['coverage_pct','all_common_rmse_cm_s','paired_rmse_cm_s','baseline_paired_rmse_cm_s',
                      'delta_rmse_cm_s','delta_ci_low_24h','delta_ci_high_24h'])}

![Accuracy and coverage](figures/accuracy_coverage.png)
![Paired differences](figures/paired_improvement.png)

## Which change matters?

Negative RMSE differences mean lower prediction error on the *same* observations. QC-only changes can leave common-bin
RMSE unchanged while changing which observations are accepted; that is selection, not improved velocity estimates.
Raising the normalized-uncertainty limit from 0.6 to 0.7 relaxes it. Adding a site-share limit is a separate restriction.

{markdown_table(ablations.set_index('change'),['delta_rmse_cm_s','ci_low','ci_high','delta_coverage_pp','paired_bins'])}

## Full-network coverage

{markdown_table(coverage,['computed_per_hour','pass_per_hour','qc_removed_pct','ever_covered'])}

![Coverage through time](figures/hourly_coverage.png)
![Coverage maps](figures/coverage_maps.png)
![Errors by site](figures/error_by_site.png)
![Errors by day](figures/error_by_day.png)

## Gap reconstruction

Daily fields require at least 12 QC-passing hourly vectors. Identical masks hide values observed by all three fields.
Five deterministic repeats test random 10% holes, 5x5-cell spatial patches, two-day outages, and two-day spatial patches.
Original hidden values are removed **before** means or fills are computed. Methods: linear temporal interpolation
(bracketing observations, maximum three-day span), nearest-neighbor averaging, and nearest-neighbor plus Gauss-Markov.
Recovery and vector RMSE are both reported; paired differences use the same recovered points.

These tests measure recovery of radar-derived daily values, not accuracy against the ocean. Repeats overlap and are
not independent samples. Gauss-Markov uses the toolbox's unchanged 4-hour scale and 8-hour search limit: on daily data,
its neighbor search uses the same day only. No temporal-scale tuning was performed. Exported fills are flagged and
masked outside the study's water grid. See `tables/gapfill_summary.csv` and `tables/gapfill_repeats.csv`.

For the full candidate's daily field, each row below compares a fill method with linear temporal interpolation
on exactly the same recovered hidden values. Errors are vector RMSE in cm/s; recovery is a percentage of hidden values.

{markdown_table(full_gaps,['recovered_pct','paired_rmse_cm_s','reference_rmse_cm_s','delta_rmse_cm_s'])}

## Interpretation and operational limits

The study does not assume that a candidate wins. Assess error, retained coverage, consistency by site/day, and both
block sizes together. A centered t-1/t/t+1 window needs one additional hour of radial availability. The gap-filling
methods and temporal interpolation use future observations and are retrospective estimates.

The additional-week test remains unrun; select any candidate using this week, then freeze the parameters before
examining another week. Independent current validation is {independent_status}. WLS and EHN algorithms were not changed.

## Reproduction

Run `python examples/maracoos_study.py --data /path/to/2026-10-06_week`. Use `--stage analysis` to rebuild tables/figures
from completed checkpoints. `--start`, `--end`, `--loso-step`, and `--output` select another archived period and output.
The manifest records input checksums, numerical source checksums, configurations and package versions. Computation
checkpoints use atomic replacement and are isolated by fingerprint. Do not compare cached runs without their manifests.
'''
    (run/'report.md').write_text(text)
