"""Tests for hfradarpy.validation (comparison of HF radar totals with drifters and buoys)."""
import numpy as np
import numpy.testing as npt
import pandas as pd
import pytest
import xarray as xr
from pyproj import Geod

from hfradarpy.validation import (complex_correlation, drifter_velocity, harmonize, sample_totals, skill_score,
                                  variance_ellipse)


def _rotate(u, v, degrees):
    """Rotate vectors counterclockwise."""
    a = np.deg2rad(degrees)
    return u * np.cos(a) - v * np.sin(a), u * np.sin(a) + v * np.cos(a)


@pytest.fixture
def series():
    rng = np.random.default_rng(7)
    t = np.arange(500)
    u = 20 * np.sin(2 * np.pi * t / 12.42) + rng.normal(0, 5, t.size)
    v = 10 * np.cos(2 * np.pi * t / 24) + rng.normal(0, 5, t.size)
    return u, v


class TestComplexCorrelation:
    def test_identical_series(self, series):
        r, angle = complex_correlation(*series, *series)
        assert r == pytest.approx(1.0)
        assert angle == pytest.approx(0.0, abs=1e-12)

    @pytest.mark.parametrize("degrees", [-30.0, 7.0, 45.0])
    def test_rotated_series(self, series, degrees):
        u2, v2 = _rotate(*series, degrees)
        r, angle = complex_correlation(*series, u2, v2)
        assert r == pytest.approx(1.0)
        assert angle == pytest.approx(degrees)

    def test_matches_kundu_formula(self, series):
        u1, v1 = series
        u2, v2 = _rotate(u1 + 3, v1 - 2, 10)
        norm = np.sqrt(np.mean(u1 ** 2 + v1 ** 2) * np.mean(u2 ** 2 + v2 ** 2))
        re = np.mean(u1 * u2 + v1 * v2) / norm
        im = np.mean(u1 * v2 - v1 * u2) / norm
        r, angle = complex_correlation(u1, v1, u2, v2)
        assert r == pytest.approx(np.hypot(re, im))
        assert angle == pytest.approx(np.rad2deg(np.arctan2(im, re)))

    def test_missing_samples_are_ignored(self, series):
        u, v = series
        u2 = u.copy()
        u2[::10] = np.nan
        r, _ = complex_correlation(u, v, u2, v)
        assert r == pytest.approx(1.0)


class TestSkillScore:
    def test_identical_series(self, series):
        skill, su, sv = skill_score(*series, *series)
        assert skill == su == sv == 1.0

    @pytest.mark.parametrize("mean_of", ["estimate", "reference"])
    def test_matches_hand_computation(self, mean_of):
        u_ref, v_ref = np.array([1.0, 2.0, 3.0]), np.array([0.0, 1.0, -1.0])
        u_est, v_est = np.array([1.5, 1.5, 2.0]), np.array([0.5, 0.5, -0.5])

        def by_hand(ref, est):
            mean = (est if mean_of == "estimate" else ref).mean()
            return 1 - np.sum((ref - est) ** 2) / np.sum((np.abs(ref - mean) + np.abs(est - mean)) ** 2)

        skill, su, sv = skill_score(u_ref, v_ref, u_est, v_est, mean_of=mean_of)
        assert su == pytest.approx(by_hand(u_ref, u_est))
        assert sv == pytest.approx(by_hand(v_ref, v_est))
        assert skill == pytest.approx((su + sv) / 2)

    def test_fang_uses_the_mean_of_the_estimate(self):
        # Fang et al. (2015) put the mean of the estimate in the denominator, unlike Willmott (1981)
        u_ref, u_est = np.array([0.0, 1.0, 2.0, 3.0]), np.array([1.0, 1.0, 3.0, 5.0])
        fang = skill_score(u_ref, u_ref, u_est, u_est)[0]
        willmott = skill_score(u_ref, u_ref, u_est, u_est, mean_of="reference")[0]
        assert fang != pytest.approx(willmott)

    def test_invalid_mean_of(self, series):
        with pytest.raises(ValueError):
            skill_score(*series, *series, mean_of="both")


class TestVarianceEllipse:
    def test_known_covariance(self):
        rng = np.random.default_rng(1)
        major, minor = rng.normal(0, 4, 20000), rng.normal(0, 1, 20000)
        u, v = _rotate(major, minor, 30.0)
        e = variance_ellipse(u + 5, v - 2)
        assert e["u_mean"] == pytest.approx(5, abs=0.1)
        assert e["v_mean"] == pytest.approx(-2, abs=0.1)
        assert e["major"] == pytest.approx(4, rel=0.03)
        assert e["minor"] == pytest.approx(1, rel=0.03)
        assert e["angle"] == pytest.approx(30, abs=1)
        assert e["n"] == 20000

    def test_angle_range(self):
        u, v = _rotate(np.array([-3.0, 3.0, -1.0, 1.0]), np.array([0.0, 0.0, 0.0, 0.0]), 120.0)
        assert variance_ellipse(u, v)["angle"] == pytest.approx(-60)


class TestDrifterVelocity:
    def test_constant_velocity_track(self):
        # A drifter moving 0.3 m/s east and 0.4 m/s north, positions every 30 minutes
        time = pd.date_range("2024-01-01", periods=10, freq="30min")
        g = Geod(ellps="GRS80")
        lon, lat = [-95.0], [28.0]
        for _ in range(9):
            lo, la, _ = g.fwd(lon[-1], lat[-1], np.rad2deg(np.arctan2(0.3, 0.4)), 0.5 * 1800)
            lon.append(lo)
            lat.append(la)
        u, v = drifter_velocity(lon, lat, time)
        assert np.isnan(u[0]) and np.isnan(u[-1])
        npt.assert_allclose(u[1:-1], 0.3, atol=1e-4)  # legs of fixed bearing are not one exact geodesic
        npt.assert_allclose(v[1:-1], 0.4, atol=1e-4)


class TestHarmonize:
    def test_speed_limit_and_moving_average(self):
        time = pd.date_range("2024-01-01", periods=7, freq="30min")
        u = np.array([1.0, 1.0, 1.0, 9.0, 1.0, 1.0, 1.0])
        v = np.zeros(7)
        out = harmonize(u, v, time, window="3h", max_speed=5)
        npt.assert_allclose(out["u"], 1.0)
        unfiltered = harmonize(u, v, time, window="3h")
        assert unfiltered["u"].iloc[3] > 1


class TestSampleTotals:
    def test_nearest_grid_point_and_time(self):
        time = pd.date_range("2024-01-01", periods=3, freq="h")
        lat, lon = np.array([28.0, 28.1]), np.array([-95.0, -94.9])
        velu = np.arange(12, dtype=float).reshape(3, 2, 2)
        ds = xr.Dataset({"VELU": (("time", "lat", "lon"), velu), "VELV": (("time", "lat", "lon"), -velu)},
                        coords=dict(time=time, lat=lat, lon=lon))
        u, v = sample_totals(ds, lon=[-94.91, -95.0], lat=[28.09, 28.0],
                             time=[time[1] + pd.Timedelta("10min"), time[2]], max_distance_km=3, max_dt="30min")
        npt.assert_array_equal(u, [velu[1, 1, 1], velu[2, 0, 0]])
        npt.assert_array_equal(v, -u)

    def test_limits(self):
        time = pd.date_range("2024-01-01", periods=2, freq="h")
        ds = xr.Dataset({"VELU": (("time", "lat", "lon"), np.ones((2, 1, 1))),
                         "VELV": (("time", "lat", "lon"), np.ones((2, 1, 1)))},
                        coords=dict(time=time, lat=[28.0], lon=[-95.0]))
        u, _ = sample_totals(ds, lon=[-95.0, -94.0, -95.0], lat=[28.0, 28.0, 28.0],
                             time=[time[0], time[0], time[1] + pd.Timedelta("5h")], max_distance_km=6, max_dt="1h")
        assert u[0] == 1 and np.isnan(u[1]) and np.isnan(u[2])
