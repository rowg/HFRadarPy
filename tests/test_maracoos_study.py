"""Check study alignment, leakage prevention, and paired uncertainty calculations."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'examples'))
from maracoos_study import BASELINE
from maracoos_study_analysis import paired_interval, skill_tables, holdout_mask


def test_withheld_site_is_excluded_from_every_hour_and_grid_subsetting_is_equivalent():
    from maracoos_study import combine, configurations
    from tests.test_totals import _radial_dataframe
    from geopandas import GeoSeries
    from shapely.geometry import Point
    t = pd.Timestamp('2025-01-01')
    inputs = {}
    for offset in (-1, 0, 1):
        stamp = t + pd.Timedelta(hours=offset)
        rdf = _radial_dataframe(uniform=(12., -7.))
        for r in rdf.Radial:
            r.time = stamp.to_pydatetime()
        inputs[stamp] = rdf.Radial.to_dict()
    lon, lat = np.meshgrid(np.arange(-74.6, -73.2, .1), np.arange(38.8, 40.3, .1))
    grid = GeoSeries([Point(x,y) for x,y in zip(lon.ravel(),lat.ravel())], crs='EPSG:4326')
    combo = configurations(dict(mdlvar=420.,errvar=66.,sx=10.,sy=10.,weighting_function='exponential'))['oi_prod_window']
    conf = {'Totals': {'MinNumSites':2,'MinNumRads':3}}
    full = combine(t, inputs, combo, grid, conf, {}, exclude='SITC')
    assert full is not None and full.data.VELU.notna().any()
    select = np.flatnonzero(full.data.VELU.notna())[::3]
    subset = grid.iloc[select].reset_index(drop=True)
    ref = full.data[['VELU','VELV']].iloc[select].to_numpy(float)
    for radials in inputs.values():
        radials['SITC'].data['VELO'] = 1e6
    result = combine(t, inputs, combo, subset, conf, {}, exclude='SITC')
    np.testing.assert_allclose(result.data[['VELU','VELV']].to_numpy(float), ref, atol=1e-9, rtol=0)


def test_paired_bootstrap_preserves_identical_predictions():
    times = pd.date_range('2026-10-06', periods=28, freq='6h').repeat(3)
    errors = np.arange(len(times)) / 10
    lo, hi, blocks = paired_interval(errors, errors.copy(), times)
    assert (lo, hi, blocks) == (0., 0., 7)


def test_bootstrap_resamples_days_not_individual_bins():
    times = pd.date_range('2026-10-06', periods=7, freq='1D')
    errors = np.array([1., 2., 3., 1., 2., 3., 1.])
    ref = errors + 1
    first = paired_interval(errors, ref, times)
    repeated = paired_interval(np.repeat(errors, 100), np.repeat(ref, 100), times.repeat(100))
    np.testing.assert_allclose(first, repeated, atol=1e-12)
    assert first[1] < 0


def test_skill_keeps_all_missing_rows_in_coverage_denominator():
    index = pd.MultiIndex.from_tuples([(pd.Timestamp('2026-10-06'), 'A', i) for i in range(4)],
                                     names=['time','site','bin'])
    data = pd.DataFrame({'observed': [10., 20., 30., 40.], 'eligible': [True]*4,
                         BASELINE: [12., 22., np.nan, np.nan],
                         'candidate': [11., np.nan, 31., np.nan],
                         'never predicts': [np.nan]*4}, index=index)
    skill, _ = skill_tables(data, [BASELINE, 'candidate', 'never predicts'])
    assert skill.loc['candidate', 'coverage_pct'] == 50
    assert skill.loc['candidate', 'paired_bins'] == 1
    assert skill.loc['candidate', 'delta_rmse_cm_s'] == -1
    assert skill.loc['candidate', 'all_common_bins'] == 0
    assert np.isnan(skill.loc['candidate', 'all_common_rmse_cm_s'])
    assert skill.loc['never predicts', 'coverage_pct'] == 0


def test_checkpoint_join_distinguishes_algorithm_prefixes_and_requires_all_cases(tmp_path):
    import json
    from maracoos_study_analysis import load_loso
    (tmp_path/'loso').mkdir()
    (tmp_path/'tables').mkdir()
    (tmp_path/'manifest.json').write_text(json.dumps({'loso_step_hours':6}))
    hours = pd.date_range('2026-10-06', periods=1)
    pd.DataFrame({'time':hours,'site':['A'],'bins':[2]}).to_csv(tmp_path/'tables/input_counts.csv',index=False)
    combos={'oi_prod': {'outputs':{BASELINE:None}}, 'oi_prod_window': {'outputs':{'window':None}}}
    for key,c in combos.items():
        np.savez(tmp_path/'loso'/f'{key}_20261006T0000_A.npz', observed=[1.,2.],eligible=[True,True],
                 **{next(iter(c['outputs'])): [3.,np.nan]})
    ctx={'run':tmp_path,'hours':hours,'combos':combos}
    joined=load_loso(ctx)
    assert len(joined)==2
    assert list(joined.columns)==['observed','eligible',BASELINE,'window']
    assert joined[BASELINE].isna().sum()==1
    (tmp_path/'loso/oi_prod_window_20261006T0000_A.npz').unlink()
    with pytest.raises(ValueError, match='No leave-one-site-out'):
        load_loso(ctx)


@pytest.mark.parametrize('kind', ['random_10pct','spatial_patches','two_day_outage','space_time_patch'])
def test_holdouts_are_reproducible_and_only_hide_observations(kind):
    observed = np.ones((7,20,20), dtype=bool)
    observed[:, :2] = False
    a = holdout_mask(observed, kind, seed=42)
    np.testing.assert_array_equal(a, holdout_mask(observed, kind, seed=42))
    assert a.any()
    assert not (a & ~observed).any()
    if kind == 'two_day_outage':
        days = np.flatnonzero(a.any(axis=(1,2)))
        assert len(days) == 2 and days[1] == days[0]+1
        assert days[0] > 0 and days[-1] < 6
