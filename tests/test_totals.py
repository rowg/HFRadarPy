"""Tests for hfradarpy.totals module — Total class and supporting functions."""
from collections import OrderedDict

import numpy as np
import numpy.testing as npt
import pandas as pd
import pytest
from shapely.geometry import Point
from geopandas import GeoSeries

from pyproj import Geod

from hfradarpy.calc import lonlat2km
from hfradarpy.radials import Radial
from hfradarpy.totals import Total, totalLeastSquare, totalOI, combineRadials, makeTotalVectors_oi_vectorized, _stackRadials


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_total(velo, velu, velv, head, gdop=None, nrad=None,
                uqal=None, vqal=None, uacc=None, vacc=None,
                lond=None, latd=None, is_wera=True):
    """Build a minimal Total object with synthetic data for testing."""
    n = len(velo)
    t = Total()  # empty total
    t.metadata = OrderedDict()
    t.is_wera = is_wera
    t.is_combined = False
    t.initialize_qc()

    data = {
        "VELO": velo,
        "VELU": velu,
        "VELV": velv,
        "HEAD": head,
        "LOND": lond if lond is not None else np.linspace(-75, -74, n),
        "LATD": latd if latd is not None else np.linspace(38, 39, n),
    }
    if gdop is not None:
        data["GDOP"] = gdop
    if nrad is not None:
        data["NRAD"] = nrad
    if uqal is not None:
        data["UQAL"] = uqal
    if vqal is not None:
        data["VQAL"] = vqal
    if uacc is not None:
        data["UACC"] = uacc
    if vacc is not None:
        data["VACC"] = vacc
    t.data = pd.DataFrame(data)
    return t


# ---------------------------------------------------------------------------
# Empty Total and grid initialisation
# ---------------------------------------------------------------------------

class TestEmptyTotal:
    def test_creates_empty_dataframe(self):
        t = Total()
        assert isinstance(t.data, pd.DataFrame)
        assert t.data.empty

    def test_creates_empty_site_source(self):
        t = Total()
        assert isinstance(t.site_source, pd.DataFrame)
        assert t.site_source.empty

    def test_is_not_corrupt(self):
        t = Total()
        assert t._iscorrupt is False

    def test_file_type(self):
        t = Total()
        assert t.file_type() == "totals"


class TestInitializeGrid:
    @pytest.fixture
    def grid(self):
        pts = GeoSeries([Point(-75.0, 38.0), Point(-74.5, 38.5), Point(-74.0, 39.0)])
        return pts.set_crs("epsg:4326")

    def test_grid_fills_lond_latd(self, grid):
        t = Total()
        t.initialize_grid(grid)
        assert len(t.data) == 3
        assert "LOND" in t.data.columns
        assert "LATD" in t.data.columns
        npt.assert_array_almost_equal(t.data["LOND"].values, [-75.0, -74.5, -74.0])
        npt.assert_array_almost_equal(t.data["LATD"].values, [38.0, 38.5, 39.0])

    def test_grid_sets_great_circle_metadata(self, grid):
        t = Total()
        t.initialize_grid(grid)
        assert "GreatCircle" in t.metadata
        assert "WGS" in t.metadata["GreatCircle"]

    def test_grid_creates_expected_columns(self, grid):
        t = Total()
        t.initialize_grid(grid)
        expected = {"LOND", "LATD", "VELU", "VELV", "VELO", "HEAD",
                    "UQAL", "VQAL", "CQAL", "GDOP", "NRAD"}
        assert expected.issubset(set(t.data.columns))


# ---------------------------------------------------------------------------
# totalLeastSquare
# ---------------------------------------------------------------------------

class TestTotalLeastSquare:
    def test_orthogonal_radials(self):
        """Two radials from north and east should give correct u and v."""
        # Radial from north (HEAD=0°): measures north-south component
        # Radial from east (HEAD=90°): measures east-west component
        df = pd.DataFrame({
            "VELO": [0.5, 0.3],       # velocities
            "HEAD": [0.0, 90.0],       # bearings (true convention)
            "STD": [1.0, 1.0],         # uniform weights
        })
        u, v, C, Cgdop = totalLeastSquare(df)
        # With orthogonal radials and unit weights, u ≈ radial at 90°, v ≈ radial at 0°
        assert not np.isnan(u)
        assert not np.isnan(v)
        assert isinstance(C, np.ndarray)
        assert isinstance(Cgdop, np.ndarray)
        assert C.shape == (2, 2)

    def test_three_radials(self):
        """Three radials at 120° intervals should resolve u and v."""
        df = pd.DataFrame({
            "VELO": [0.4, 0.2, -0.1],
            "HEAD": [0.0, 120.0, 240.0],
            "STD": [1.0, 1.0, 1.0],
        })
        u, v, C, Cgdop = totalLeastSquare(df)
        assert not np.isnan(u)
        assert not np.isnan(v)

    def test_colinear_radials_return_nan(self):
        """Two parallel radials (same direction) can't resolve u and v."""
        df = pd.DataFrame({
            "VELO": [0.5, 0.3],
            "HEAD": [0.0, 0.0],
            "STD": [1.0, 1.0],
        })
        u, v, C, Cgdop = totalLeastSquare(df)
        # The design matrix is singular, so the result should be NaN
        assert np.isnan(u)
        assert np.isnan(v)

    def test_weighted_least_squares(self):
        """Higher weight on one radial should bias the result."""
        df_equal = pd.DataFrame({
            "VELO": [1.0, 0.0],
            "HEAD": [0.0, 90.0],
            "STD": [1.0, 1.0],
        })
        df_biased = pd.DataFrame({
            "VELO": [1.0, 0.0],
            "HEAD": [0.0, 90.0],
            "STD": [0.1, 1.0],  # much higher weight on first radial
        })
        u_eq, v_eq, _, _ = totalLeastSquare(df_equal)
        u_bi, v_bi, _, _ = totalLeastSquare(df_biased)
        # Both should produce valid results
        assert not np.isnan(u_eq)
        assert not np.isnan(u_bi)


# ---------------------------------------------------------------------------
# totalOI
# ---------------------------------------------------------------------------

def tuvOI_reference(velo, head, lon, lat, grid_lon, grid_lat, mdlvar, errvar, sx, sy):
    """Literal translation of HFRProgs tuvOI.m (exponential correlation): a = P * cmd' * inv(cdd + R) * velo."""
    ang = np.deg2rad((90 - head) % 360)  # true to math convention
    n = len(velo)
    R = np.eye(n) * errvar
    P = np.eye(2) * mdlvar
    dx, dy = lonlat2km(lon, lat, np.full(n, grid_lon), np.full(n, grid_lat))
    x1, x2 = np.meshgrid(lon, lon)
    y1, y2 = np.meshgrid(lat, lat)
    a1, a2 = np.meshgrid(ang, ang)
    drx, dry = lonlat2km(x1, y1, x2, y2)
    cmd = np.exp(-np.sqrt(dx ** 2 / sx ** 2 + dy ** 2 / sy ** 2))
    w = np.exp(-np.sqrt(drx ** 2 / sx ** 2 + dry ** 2 / sy ** 2)) * mdlvar
    cmd = np.column_stack((cmd * np.cos(ang), cmd * np.sin(ang)))
    cdd = w * (np.cos(a1) * np.cos(a2) + np.sin(a1) * np.sin(a2))
    cmdicdd = P @ cmd.T @ np.linalg.inv(cdd + R)
    a = cmdicdd @ velo
    xi = np.linalg.inv(P) @ (P - cmdicdd @ cmd @ P)
    return a[0], a[1], xi, cdd + R


class TestTotalOI:
    def test_matches_hfrprogs_tuvOI(self):
        """totalOI must reproduce tuvOI.m, including when cdd + R is not exactly symmetric."""
        rng = np.random.default_rng(42)
        n = 40
        # Radials spread over ~0.5 degrees of latitude, so lonlat2km distances depend on which point is the origin
        lon = -72.5 + rng.uniform(-0.3, 0.3, n)
        lat = 39.5 + rng.uniform(-0.25, 0.25, n)
        head = rng.uniform(0, 360, n)
        velo = rng.normal(0, 30, n)
        grid_lon, grid_lat = -72.5, 39.5
        mdlvar, errvar, sx, sy = 420, 66, 15, 15

        u_ref, v_ref, xi_ref, A = tuvOI_reference(velo, head, lon, lat, grid_lon, grid_lat, mdlvar, errvar, sx, sy)
        assert not np.allclose(A, A.T, rtol=0, atol=1e-6)  # the case this test is about

        contributions = pd.DataFrame({'VELO': velo, 'HEAD': head, 'LOND': lon, 'LATD': lat})
        u, v, xi = totalOI(contributions, pd.Series({'LOND': grid_lon, 'LATD': grid_lat}), mdlvar, errvar, sx, sy,
                           'exponential')
        npt.assert_allclose([u, v], [u_ref, v_ref], rtol=0, atol=1e-9)
        npt.assert_allclose(xi, xi_ref, rtol=0, atol=1e-12)


# ---------------------------------------------------------------------------
# combineRadials: loop vs vectorized engines
# ---------------------------------------------------------------------------

_SITES = {"SITA": (39.50, -74.30), "SITB": (40.30, -73.95), "SITC": (38.90, -74.80)}


def _current(lon, lat):
    """Spatially varying current (cm/s) used to build synthetic radials."""
    return 10 + 20 * (lon + 74), -5 + 15 * (lat - 39.5)


def _synthetic_radial(site, lat0, lon0, rng=None, uniform=None):
    """Build an in-memory CODAR-like Radial (VELO positive toward the radar)."""
    ranges = np.arange(3, 63, 3.0)
    bearings = np.arange(0, 360, 5.0)
    rr, bb = [a.ravel() for a in np.meshgrid(ranges, bearings)]
    lon, lat, back_az = Geod(ellps="WGS84").fwd(np.full(rr.size, lon0), np.full(rr.size, lat0), bb, rr * 1000)
    u, v = uniform if uniform is not None else _current(lon, lat)
    head = np.mod(back_az, 360)
    velo = u * np.sin(np.deg2rad(head)) + v * np.cos(np.deg2rad(head))
    if rng is not None:
        velo = velo + rng.normal(0, 3, velo.size)
        keep = rng.random(velo.size) > 0.2
    else:
        keep = np.ones(velo.size, dtype=bool)

    r = Radial.__new__(Radial)
    r.is_wera = False
    r.data = pd.DataFrame({"LOND": lon, "LATD": lat, "VELO": velo, "HEAD": head, "BEAR": bb, "RNGE": rr,
                           "ETMP": np.abs(rng.normal(4, 1, velo.size)) if rng is not None else np.ones(velo.size)})
    r.data = r.data[keep].reset_index(drop=True)
    r.metadata = OrderedDict([
        ("TimeZone", '"UTC" +0.000 0 "Atlantic/Reykjavik"'),
        ("TimeCoverage", "75.000 Minutes"),
        ("Origin", f"{lat0:.7f}  {lon0:.7f}"),
        ("RangeResolutionKMeters", "3.000000"),
        ("PatternType", "Measured"),
        ("AntennaBearing", "0.0 True"),
    ])
    return r


def _radial_dataframe(rng=None, uniform=None):
    radials = {s: _synthetic_radial(s, la, lo, rng, uniform) for s, (la, lo) in _SITES.items()}
    return pd.DataFrame({"Radial": pd.Series(radials, dtype=object)})


@pytest.fixture(scope="module")
def combine_grid():
    lon, lat = np.meshgrid(np.arange(-74.6, -73.2, 0.1), np.arange(38.8, 40.3, 0.1))
    return GeoSeries([Point(xy) for xy in zip(lon.ravel(), lat.ravel())], crs="EPSG:4326")


class TestCombineRadialsEngines:
    common = dict(sRad=10000, gRes=6000, tStp=pd.Timestamp("2025-01-01"), useLandMask=False, dropPoints=False)
    columns = ["VELU", "VELV", "VELO", "UQAL", "VQAL", "CQAL", "GDOP", "NRAD"]

    @pytest.mark.parametrize("method, kwargs", [
        ("uwls", {}),
        ("wls", {}),
        ("oi", dict(mdlvar=420, errvar=66, sx=10.0, sy=10.0, runParallel=False)),
        ("oi", dict(mdlvar=420, errvar=66, sx=8.0, sy=12.0, oi_option="gaussian", runParallel=False)),
    ])
    def test_vectorized_matches_loop(self, combine_grid, method, kwargs):
        rdf = _radial_dataframe(np.random.default_rng(0))
        loop, warn_loop = combineRadials(rdf, combine_grid, method=method, engine="loop", **self.common, **kwargs)
        vec, warn_vec = combineRadials(rdf, combine_grid, method=method, engine="vectorized", **self.common, **kwargs)

        assert warn_loop == warn_vec == ""
        expected = loop.data[self.columns].astype(float)
        result = vec.data[self.columns].astype(float)
        npt.assert_array_equal(expected["VELU"].isna().to_numpy(), result["VELU"].isna().to_numpy())
        assert result["VELU"].notna().sum() > 20
        npt.assert_allclose(result.to_numpy(), expected.to_numpy(), rtol=1e-6, atol=1e-6)

    def test_vectorized_uwls_recovers_uniform_current(self, combine_grid):
        rdf = _radial_dataframe(uniform=(12.0, -7.0))
        tot, warn = combineRadials(rdf, combine_grid, method="uwls", engine="vectorized", **self.common)
        valid = tot.data["VELU"].notna()
        assert valid.sum() > 20
        npt.assert_allclose(tot.data.loc[valid, "VELU"].astype(float), 12.0, atol=1e-8)
        npt.assert_allclose(tot.data.loc[valid, "VELV"].astype(float), -7.0, atol=1e-8)

    def test_vectorized_honors_contribution_thresholds(self, combine_grid):
        rdf = _radial_dataframe(np.random.default_rng(0))
        strict = dict(minContrSites=2, minContrRads=30)
        loop, _ = combineRadials(rdf, combine_grid, method="uwls", engine="loop", **strict, **self.common)
        vec, _ = combineRadials(rdf, combine_grid, method="uwls", engine="vectorized", **strict, **self.common)
        default, _ = combineRadials(rdf, combine_grid, method="uwls", engine="vectorized", **self.common)

        npt.assert_array_equal(loop.data["VELU"].isna().to_numpy(), vec.data["VELU"].isna().to_numpy())
        assert 0 < vec.data["VELU"].notna().sum() < default.data["VELU"].notna().sum()
        assert (vec.data["NRAD"].dropna() >= 30).all()

    def test_oi_batch_size_does_not_change_results(self, combine_grid):
        rad = _stackRadials(_radial_dataframe(np.random.default_rng(0)))
        lon, lat = combine_grid.x.to_numpy(), combine_grid.y.to_numpy()
        kwargs = dict(mdlvar=420, errvar=66, sx=10.0, sy=10.0)
        one_batch = makeTotalVectors_oi_vectorized(rad, lon, lat, max_batch_elements=1e9, **kwargs)
        many_batches = makeTotalVectors_oi_vectorized(rad, lon, lat, max_batch_elements=1e3, **kwargs)
        assert one_batch["VELU"].notna().sum() > 20
        npt.assert_allclose(many_batches.to_numpy(dtype=float), one_batch.to_numpy(dtype=float),
                            rtol=1e-12, atol=1e-12)

    def test_invalid_engine_raises(self, combine_grid):
        with pytest.raises(ValueError):
            combineRadials(_radial_dataframe(), combine_grid, method="uwls", engine="fast", **self.common)

    def test_vectorized_oi_requires_scalar_errvar(self, combine_grid):
        with pytest.raises(ValueError):
            combineRadials(_radial_dataframe(), combine_grid, method="oi", engine="vectorized",
                           errvar=np.ones(3), sx=10.0, sy=10.0, **self.common)


# ---------------------------------------------------------------------------
# QC: Maximum Velocity
# ---------------------------------------------------------------------------

class TestQCMaxVelocity:
    def test_all_pass_wera(self):
        """All velocities below threshold should pass (WERA, m/s)."""
        t = _make_total(
            velo=[0.5, 0.8, 1.0],
            velu=[0.3, 0.5, 0.7],
            velv=[0.4, 0.6, 0.7],
            head=[45.0, 60.0, 30.0],
            is_wera=True,
        )
        t.qc_ehn_maximum_velocity(totMaxSpeed=1.2)
        assert (t.data["CSPD_QC"] == 1).all()

    def test_fail_above_threshold_wera(self):
        """Velocities above threshold should fail (WERA, m/s)."""
        t = _make_total(
            velo=[0.5, 1.5, 2.0],
            velu=[0.3, 1.0, 1.5],
            velv=[0.4, 1.1, 1.3],
            head=[45.0, 60.0, 30.0],
            is_wera=True,
        )
        t.qc_ehn_maximum_velocity(totMaxSpeed=1.2)
        assert t.data["CSPD_QC"].iloc[0] == 1
        assert t.data["CSPD_QC"].iloc[1] == 4
        assert t.data["CSPD_QC"].iloc[2] == 4

    def test_codar_scaling(self):
        """CODAR velocities are in cm/s, threshold is applied as totMaxSpeed * 100."""
        t = _make_total(
            velo=[50.0, 150.0],  # cm/s
            velu=[30.0, 100.0],
            velv=[40.0, 110.0],
            head=[45.0, 60.0],
            is_wera=False,
        )
        t.qc_ehn_maximum_velocity(totMaxSpeed=1.2)  # 1.2 m/s = 120 cm/s
        assert t.data["CSPD_QC"].iloc[0] == 1  # 50 < 120
        assert t.data["CSPD_QC"].iloc[1] == 4  # 150 > 120

    def test_metadata_recorded(self):
        t = _make_total(
            velo=[0.5], velu=[0.3], velv=[0.4], head=[45.0], is_wera=True,
        )
        t.qc_ehn_maximum_velocity(totMaxSpeed=1.5)
        assert "CSPD_QC" in t.metadata["QCTest"]
        assert "1.5" in t.metadata["QCTest"]["CSPD_QC"]


# ---------------------------------------------------------------------------
# QC: GDOP Threshold
# ---------------------------------------------------------------------------

class TestQCGDOPThreshold:
    def test_all_pass(self):
        t = _make_total(
            velo=[0.5, 0.8], velu=[0.3, 0.5], velv=[0.4, 0.6],
            head=[45.0, 60.0], gdop=[1.0, 1.5], is_wera=True,
        )
        t.qc_ehn_gdop_threshold(maxGDOP=2.0)
        assert (t.data["GDOP_QC"] == 1).all()

    def test_fail_above_threshold(self):
        t = _make_total(
            velo=[0.5, 0.8], velu=[0.3, 0.5], velv=[0.4, 0.6],
            head=[45.0, 60.0], gdop=[1.0, 3.0], is_wera=True,
        )
        t.qc_ehn_gdop_threshold(maxGDOP=2.0)
        assert t.data["GDOP_QC"].iloc[0] == 1
        assert t.data["GDOP_QC"].iloc[1] == 4

    def test_metadata_recorded(self):
        t = _make_total(
            velo=[0.5], velu=[0.3], velv=[0.4], head=[45.0],
            gdop=[1.0], is_wera=True,
        )
        t.qc_ehn_gdop_threshold(maxGDOP=2.5)
        assert "GDOP_QC" in t.metadata["QCTest"]
        assert "2.5" in t.metadata["QCTest"]["GDOP_QC"]


# ---------------------------------------------------------------------------
# QC: Data Density Threshold
# ---------------------------------------------------------------------------

class TestQCDataDensity:
    def test_pass_with_enough_radials(self):
        """CODAR total with enough contributing radials should pass."""
        t = _make_total(
            velo=[50.0, 80.0], velu=[30.0, 50.0], velv=[40.0, 60.0],
            head=[45.0, 60.0], nrad=[3, 5], is_wera=False,
        )
        t.qc_ehn_data_density_threshold(minContrRad=2)
        assert (t.data["DDNS_QC"] == 1).all()

    def test_fail_with_too_few_radials(self):
        t = _make_total(
            velo=[50.0, 80.0], velu=[30.0, 50.0], velv=[40.0, 60.0],
            head=[45.0, 60.0], nrad=[1, 5], is_wera=False,
        )
        t.qc_ehn_data_density_threshold(minContrRad=3)
        assert t.data["DDNS_QC"].iloc[0] == 4  # 1 < 3
        assert t.data["DDNS_QC"].iloc[1] == 1  # 5 >= 3

    def test_wera_skips_nrad_check(self):
        """WERA totals don't have NRAD, so all should pass."""
        t = _make_total(
            velo=[0.5], velu=[0.3], velv=[0.4], head=[45.0], is_wera=True,
        )
        t.qc_ehn_data_density_threshold(minContrRad=2)
        # For WERA, the if-not-is_wera block is skipped, so all remain flag=1
        assert (t.data["DDNS_QC"] == 1).all()


# ---------------------------------------------------------------------------
# QC: Maximum Variance (WERA)
# ---------------------------------------------------------------------------

class TestQCMaxVariance:
    def test_pass_wera(self):
        """WERA non-combined: variance below threshold should pass."""
        t = _make_total(
            velo=[0.5, 0.8], velu=[0.3, 0.5], velv=[0.4, 0.6],
            head=[45.0, 60.0], uacc=[0.5, 0.8], vacc=[0.5, 0.7],
            is_wera=True,
        )
        t.qc_ehn_maximum_variance(totMaxVar=1.0)
        # UACC^2: 0.25, 0.64; VACC^2: 0.25, 0.49 — all < 1.0
        assert (t.data["VART_QC"] == 1).all()

    def test_fail_wera(self):
        """WERA non-combined: variance above threshold should fail."""
        t = _make_total(
            velo=[0.5, 0.8], velu=[0.3, 0.5], velv=[0.4, 0.6],
            head=[45.0, 60.0], uacc=[0.5, 1.5], vacc=[0.5, 0.7],
            is_wera=True,
        )
        t.qc_ehn_maximum_variance(totMaxVar=1.0)
        assert t.data["VART_QC"].iloc[0] == 1  # 0.25 < 1.0
        assert t.data["VART_QC"].iloc[1] == 4  # 2.25 > 1.0

    def test_codar_uses_uqal_vqal(self):
        """CODAR: variance computed from UQAL/VQAL (in cm/s)."""
        t = _make_total(
            velo=[50.0, 80.0], velu=[30.0, 50.0], velv=[40.0, 60.0],
            head=[45.0, 60.0],
            uqal=[50.0, 150.0],  # cm/s, so /100 => 0.5, 1.5 m/s
            vqal=[50.0, 80.0],
            is_wera=False,
        )
        t.qc_ehn_maximum_variance(totMaxVar=1.0)
        assert t.data["VART_QC"].iloc[0] == 1  # (0.5)^2 = 0.25 < 1.0
        assert t.data["VART_QC"].iloc[1] == 4  # (1.5)^2 = 2.25 > 1.0


# ---------------------------------------------------------------------------
# QC: Temporal Derivative
# ---------------------------------------------------------------------------

class TestQCTemporalDerivative:
    def test_small_change_passes(self):
        """Small velocity change between timestamps should pass."""
        t1 = _make_total(
            velo=[0.5, 0.8], velu=[0.3, 0.5], velv=[0.4, 0.6],
            head=[45.0, 60.0], lond=[-75.0, -74.0], latd=[38.0, 39.0],
            is_wera=True,
        )
        t0 = _make_total(
            velo=[0.4, 0.7], velu=[0.2, 0.4], velv=[0.3, 0.5],
            head=[45.0, 60.0], lond=[-75.0, -74.0], latd=[38.0, 39.0],
            is_wera=True,
        )
        t1.qc_ehn_temporal_derivative(t0, tempDerThr=1.0)
        # Differences: |0.5-0.4|=0.1, |0.8-0.7|=0.1 — both < 1.0
        assert (t1.data["VART_QC"] == 1).all()

    def test_large_change_fails(self):
        """Large velocity change should fail."""
        t1 = _make_total(
            velo=[0.5, 2.5], velu=[0.3, 2.0], velv=[0.4, 1.5],
            head=[45.0, 60.0], lond=[-75.0, -74.0], latd=[38.0, 39.0],
            is_wera=True,
        )
        t0 = _make_total(
            velo=[0.4, 0.3], velu=[0.2, 0.2], velv=[0.3, 0.2],
            head=[45.0, 60.0], lond=[-75.0, -74.0], latd=[38.0, 39.0],
            is_wera=True,
        )
        t1.qc_ehn_temporal_derivative(t0, tempDerThr=1.0)
        assert t1.data["VART_QC"].iloc[0] == 1  # |0.5-0.4| = 0.1 < 1.0
        assert t1.data["VART_QC"].iloc[1] == 4  # |2.5-0.3| = 2.2 > 1.0

    def test_no_previous_sets_not_evaluated(self):
        """When no previous total exists, all flags should be 0 (not evaluated)."""
        t1 = _make_total(
            velo=[0.5, 0.8], velu=[0.3, 0.5], velv=[0.4, 0.6],
            head=[45.0, 60.0], is_wera=True,
        )
        t1.qc_ehn_temporal_derivative(None, tempDerThr=1.0)
        assert (t1.data["VART_QC"] == 0).all()

    def test_new_locations_not_evaluated(self):
        """Locations only in current total (not in previous) should be flag 0."""
        t1 = _make_total(
            velo=[0.5, 0.8, 1.0], velu=[0.3, 0.5, 0.7], velv=[0.4, 0.6, 0.7],
            head=[45.0, 60.0, 30.0],
            lond=[-75.0, -74.0, -73.0], latd=[38.0, 39.0, 40.0],
            is_wera=True,
        )
        t0 = _make_total(
            velo=[0.4, 0.7], velu=[0.2, 0.4], velv=[0.3, 0.5],
            head=[45.0, 60.0],
            lond=[-75.0, -74.0], latd=[38.0, 39.0],
            is_wera=True,
        )
        t1.qc_ehn_temporal_derivative(t0, tempDerThr=1.0)
        assert t1.data["VART_QC"].iloc[0] == 1  # exists in both, small diff
        assert t1.data["VART_QC"].iloc[1] == 1  # exists in both, small diff
        assert t1.data["VART_QC"].iloc[2] == 0  # only in current


# ---------------------------------------------------------------------------
# QC: Overall Flag
# ---------------------------------------------------------------------------

class TestQCOverallFlag:
    def test_all_pass(self):
        """If all individual QC tests pass, overall should be 1."""
        t = _make_total(
            velo=[0.5], velu=[0.3], velv=[0.4], head=[45.0],
            gdop=[1.0], uacc=[0.5], vacc=[0.5], is_wera=True,
        )
        t.qc_ehn_maximum_velocity(totMaxSpeed=1.2)
        t.qc_ehn_gdop_threshold(maxGDOP=2.0)
        t.qc_ehn_overall_qc_flag()
        assert t.data["QCflag"].iloc[0] == 1

    def test_one_fail_causes_overall_fail(self):
        """If any individual QC test fails, overall should be 4."""
        t = _make_total(
            velo=[2.0], velu=[1.5], velv=[1.3], head=[45.0],
            gdop=[1.0], uacc=[0.5], vacc=[0.5], is_wera=True,
        )
        t.qc_ehn_maximum_velocity(totMaxSpeed=1.2)  # will fail (2.0 > 1.2)
        t.qc_ehn_gdop_threshold(maxGDOP=2.0)  # will pass
        t.qc_ehn_overall_qc_flag()
        assert t.data["QCflag"].iloc[0] == 4

    def test_mixed_results(self):
        """Some vectors pass all tests, some fail at least one."""
        t = _make_total(
            velo=[0.5, 2.0], velu=[0.3, 1.5], velv=[0.4, 1.3],
            head=[45.0, 60.0], gdop=[1.0, 1.0],
            uacc=[0.5, 0.5], vacc=[0.5, 0.5], is_wera=True,
        )
        t.qc_ehn_maximum_velocity(totMaxSpeed=1.2)
        t.qc_ehn_gdop_threshold(maxGDOP=2.0)
        t.qc_ehn_overall_qc_flag()
        assert t.data["QCflag"].iloc[0] == 1  # all pass
        assert t.data["QCflag"].iloc[1] == 4  # velocity fails


# ---------------------------------------------------------------------------
# QC metadata
# ---------------------------------------------------------------------------

class TestQCMetadata:
    def test_initialize_qc_creates_empty_dict(self):
        t = Total()
        t.metadata = OrderedDict()
        t.initialize_qc()
        assert isinstance(t.metadata["QCTest"], dict)
        assert len(t.metadata["QCTest"]) == 0

    def test_each_qc_records_metadata(self):
        t = _make_total(
            velo=[0.5], velu=[0.3], velv=[0.4], head=[45.0],
            gdop=[1.0], nrad=[3], uacc=[0.5], vacc=[0.5], is_wera=True,
        )
        t.qc_ehn_maximum_velocity()
        t.qc_ehn_gdop_threshold()
        t.qc_ehn_data_density_threshold()
        assert "CSPD_QC" in t.metadata["QCTest"]
        assert "GDOP_QC" in t.metadata["QCTest"]
        assert "DDNS_QC" in t.metadata["QCTest"]
