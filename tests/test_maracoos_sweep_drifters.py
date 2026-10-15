"""Check the sweep's QC scoring and the drifter processing of the MARACOOS example scripts."""
import sys
from pathlib import Path

import numpy as np
import numpy.testing as npt
import pandas as pd
from pyproj import Geod

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'examples'))
from maracoos_drifters import drifter_samples
from maracoos_sweep import accepted, configurations, curves


def test_sweep_includes_production_and_varies_one_setting_at_a_time():
    configs = configurations()
    assert configs['sx15_e66'] == dict(sx=15.0, errvar=66.0, window=0, time_scale=None, errcol=None,
                                       family='length scale x variance ratio')
    assert configs['sx15_e66_window_T6']['time_scale'] == 6.0
    assert configs['sx10_etmp']['errcol'] == 'ETMP'


def test_qc_thresholds_are_scored_from_stored_predictions():
    d = pd.DataFrame({'config': 'a', 'observed': [0.0, 0.0, 0.0, 0.0], 'pred': [1.0, 2.0, 3.0, np.nan],
                      'uerr': [0.2, 0.5, 0.9, 0.1], 'msfr': [0.6, 0.8, 0.5, 0.5], 'ror': [1.0, 5.0, 1.0, 1.0]})
    assert accepted(d, 0.6).tolist() == [True, True, False, False]
    assert accepted(d, 1.0, share=0.75).tolist() == [True, False, True, False]
    assert accepted(d, 1.0, ror=4).tolist() == [True, False, True, False]
    c = curves(d, n_eligible=4).set_index('uncertainty')
    assert c.loc[0.6, 'coverage_pct'] == 50
    assert c.loc[0.6, 'rmse_cm_s'] == np.sqrt((1 + 4) / 2)


def test_drifter_velocities_and_fixed_platforms(tmp_path):
    hours = pd.date_range('2026-10-06', periods=12, freq='h')
    g = Geod(ellps='GRS80')
    lon, lat = [-70.4], [40.8]
    for _ in range(len(hours) + 5):  # 20 cm/s toward the northeast, fixes every hour
        lo, la, _ = g.fwd(lon[-1], lat[-1], 45.0, 0.2 * 3600)
        lon.append(lo)
        lat.append(la)
    times = pd.date_range(hours[0] - pd.Timedelta('3h'), periods=len(lon), freq='h')
    rows = [dict(platform_code='moving', platform_type='DRIFTING BUOYS', time=f'{t:%Y-%m-%dT%H:%M:%SZ}',
                 latitude=la, longitude=lo) for t, lo, la in zip(times, lon, lat)]
    rows += [dict(platform_code='fixed', platform_type='DRIFTING BUOYS', time=f'{t:%Y-%m-%dT%H:%M:%SZ}',
                  latitude=36.18, longitude=-75.75) for t in times]
    path = tmp_path / 'drifters.csv'
    pd.concat([pd.DataFrame([dict(platform_code='', platform_type='', time='UTC', latitude='degrees_north',
                                  longitude='degrees_east')]), pd.DataFrame(rows)]).to_csv(path, index=False)

    s = drifter_samples(path, hours)
    assert set(s['platform_code']) == {'moving'}
    npt.assert_allclose(s['u'], 20 * np.sin(np.pi / 4), atol=0.05)
    npt.assert_allclose(s['v'], 20 * np.cos(np.pi / 4), atol=0.05)
