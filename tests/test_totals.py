"""Tests for hfradarpy.totals module — Total class and supporting functions."""
import datetime as dt
from collections import OrderedDict

import numpy as np
import numpy.testing as npt
import pandas as pd
import pytest
import xarray as xr
from shapely.geometry import Point
from geopandas import GeoSeries

from pyproj import Geod

from hfradarpy.calc import lonlat2km
from hfradarpy.radials import Radial
from hfradarpy.totals import Total, totalLeastSquare, totalUnWeightedLeastSquare, totalOI, combineRadials, \
    makeTotalVectors_oi_vectorized, _stackRadials, _siteContribution


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

def _radial_components(u, v, head):
    """Project a (u, v) current onto radial headings (degrees true): VELO = u*sin(HEAD) + v*cos(HEAD)."""
    head = np.deg2rad(np.asarray(head, dtype=float))
    return u * np.sin(head) + v * np.cos(head)


class TestTotalXarray:
    @pytest.mark.parametrize("lon, lat", [
        ([-74.0, -74.3, -74.1], [39.3, 39.0, 39.0]),
        ([-74.0, -74.0], [39.0, 39.2]),
        ([-74.0, -74.2], [39.0, 39.0]),
        ([-74.0], [39.0]),
        ([179.9, -179.9], [39.0, 39.2]),
    ])
    def test_default_grid_preserves_coordinates_and_values(self, lon, lat):
        total = Total(grid=GeoSeries([Point(x, y) for x, y in zip(lon, lat)], crs="EPSG:4326"))
        total.time = dt.datetime(2025, 1, 1)
        total.data["VELU"] = np.arange(len(lon), dtype=float) + 1
        total.data["VELV"] = -total.data["VELU"]
        before = total.data.copy(deep=True)

        ds = total.to_xarray()

        assert dict(ds.sizes) == {"time": 1, "depth": 1, "lat": len(set(lat)), "lon": len(set(lon))}
        assert pd.Timestamp(ds.time.values[0]) == total.time
        assert ds.VELU.dims == ("time", "depth", "lat", "lon")
        assert int(ds.VELU.notnull().sum()) == len(lon)
        for i, (x, y) in enumerate(zip(lon, lat)):
            cell = ds.sel(lon=x, lat=y)
            assert cell.VELU.item() == i + 1
            assert cell.VELV.item() == -(i + 1)
        pd.testing.assert_frame_equal(total.data, before)

    @pytest.mark.parametrize("grid_kind", ["geoseries", "file"])
    def test_explicit_grid_keeps_unobserved_cells(self, grid_kind, tmp_path):
        grid = GeoSeries([Point(x, y) for y in (39.0, 39.5) for x in (-74.5, -74.0)], crs="EPSG:4326")
        total = Total(grid=grid)
        total.time = dt.datetime(2025, 1, 1)
        total.data["VELU"] = [1., 2., 3., 4.]
        expected = total.to_xarray()
        total.data = total.data.iloc[[0, 3]].copy()
        expected["VELU"] = expected.VELU.where(expected.VELU.isin([1., 4.]))
        if grid_kind == "file":
            path = tmp_path / "grid.txt"
            np.savetxt(path, np.column_stack((grid.x, grid.y)))
            grid = str(path)

        actual = total.to_xarray(grid=grid)

        xr.testing.assert_allclose(actual, expected)


class TestTotalLeastSquare:
    def test_orthogonal_radials(self):
        """A radial heading north measures v and one heading east measures u."""
        df = pd.DataFrame({"VELO": [0.5, 0.3], "HEAD": [0.0, 90.0], "STD": [1.0, 1.0]})
        u, v, C, Cgdop = totalLeastSquare(df)
        assert u == pytest.approx(0.3, abs=1e-12)
        assert v == pytest.approx(0.5, abs=1e-12)
        npt.assert_allclose(C, np.eye(2), atol=1e-12)
        npt.assert_allclose(Cgdop, np.eye(2), atol=1e-12)

    def test_recovers_known_current(self):
        """Noise-free radials from any non-degenerate geometry recover u and v exactly, whatever the weights."""
        head = np.array([0.0, 90.0, 200.0, 300.0])
        std = np.array([1.0, 2.0, 0.5, 1.0])
        df = pd.DataFrame({"VELO": _radial_components(0.3, -0.2, head), "HEAD": head, "STD": std})
        u, v, C, Cgdop = totalLeastSquare(df)
        assert u == pytest.approx(0.3, abs=1e-12)
        assert v == pytest.approx(-0.2, abs=1e-12)

        # Covariances computed independently from the design matrix rows [sin(HEAD), cos(HEAD)]
        G = np.column_stack((np.sin(np.deg2rad(head)), np.cos(np.deg2rad(head))))
        npt.assert_allclose(C, np.linalg.inv(G.T @ np.diag(1 / std ** 2) @ G), atol=1e-12)
        npt.assert_allclose(Cgdop, np.linalg.inv(G.T @ G), atol=1e-12)

    def test_covariance_scales_with_std_but_gdop_does_not(self):
        df = pd.DataFrame({"VELO": [0.5, 0.3], "HEAD": [0.0, 90.0], "STD": [2.0, 2.0]})
        _, _, C, Cgdop = totalLeastSquare(df)
        npt.assert_allclose(C, 4 * np.eye(2), atol=1e-12)
        npt.assert_allclose(Cgdop, np.eye(2), atol=1e-12)

    @pytest.mark.parametrize("std, expected_v", [
        ([1.0, 1.0, 1.0], 0.5),           # equal weights: v is the plain mean of 1.0 and 0.0
        ([0.1, 1.0, 1.0], 100 / 101),     # weights 1/STD^2 = 100 and 1: v is the weighted mean
    ])
    def test_weights_are_inverse_variance(self, std, expected_v):
        """Two conflicting northward radials are averaged with weights 1/STD^2."""
        df = pd.DataFrame({"VELO": [1.0, 0.0, 0.0], "HEAD": [0.0, 0.0, 90.0], "STD": std})
        u, v, _, _ = totalLeastSquare(df)
        assert u == pytest.approx(0.0, abs=1e-12)
        assert v == pytest.approx(expected_v, abs=1e-12)

    def test_colinear_radials_return_nan(self):
        """Two parallel radials (same direction) can't resolve u and v."""
        df = pd.DataFrame({"VELO": [0.5, 0.3], "HEAD": [0.0, 0.0], "STD": [1.0, 1.0]})
        u, v, C, Cgdop = totalLeastSquare(df)
        assert np.isnan(u)
        assert np.isnan(v)


class TestTotalUnWeightedLeastSquare:
    def test_recovers_known_current(self):
        head = np.array([10.0, 100.0, 230.0])
        df = pd.DataFrame({"VELO": _radial_components(25.0, -12.0, head), "HEAD": head})
        u, v, C, Cgdop = totalUnWeightedLeastSquare(df)
        assert u == pytest.approx(25.0, abs=1e-10)
        assert v == pytest.approx(-12.0, abs=1e-10)
        G = np.column_stack((np.sin(np.deg2rad(head)), np.cos(np.deg2rad(head))))
        npt.assert_allclose(C, np.linalg.inv(G.T @ G), atol=1e-12)
        npt.assert_array_equal(C, Cgdop)

    def test_conflicting_radials_are_averaged(self):
        df = pd.DataFrame({"VELO": [1.0, 0.0, 0.0], "HEAD": [0.0, 0.0, 90.0]})
        u, v, _, _ = totalUnWeightedLeastSquare(df)
        assert u == pytest.approx(0.0, abs=1e-12)
        assert v == pytest.approx(0.5, abs=1e-12)

    def test_colinear_radials_return_nan(self):
        df = pd.DataFrame({"VELO": [0.5, 0.3], "HEAD": [45.0, 225.0]})
        u, v, C, Cgdop = totalUnWeightedLeastSquare(df)
        assert np.isnan(u)
        assert np.isnan(v)


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


def _synthetic_radial(site, lat0, lon0, rng=None, uniform=None, time=None):
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
    r.time = (time if time is not None else pd.Timestamp("2025-01-01")).to_pydatetime()
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

    @pytest.mark.parametrize("engine", ["loop", "vectorized"])
    @pytest.mark.parametrize("method", ["uwls", "wls"])
    def test_recovers_uniform_current(self, combine_grid, method, engine):
        """Noise-free radials of a uniform current must combine back to that current at every grid point."""
        rdf = _radial_dataframe(uniform=(12.0, -7.0))
        tot, warn = combineRadials(rdf, combine_grid, method=method, engine=engine, **self.common)
        assert warn == ""
        valid = tot.data["VELU"].notna()
        assert valid.sum() > 20
        npt.assert_allclose(tot.data.loc[valid, "VELU"].astype(float), 12.0, atol=1e-8)
        npt.assert_allclose(tot.data.loc[valid, "VELV"].astype(float), -7.0, atol=1e-8)
        npt.assert_allclose(tot.data.loc[valid, "VELO"].astype(float), np.hypot(12.0, -7.0), atol=1e-8)
        npt.assert_allclose(tot.data.loc[valid, "HEAD"].astype(float),
                            np.rad2deg(np.arctan2(12.0, -7.0)) % 360, atol=1e-6)

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

    @staticmethod
    def _hfrprogs_flat_counts(rdf, grid, radius_km=10):
        """Radials within radius_km of each grid point by HFRProgs lonlat2dist (makeTotals TotalsNumRads)."""
        rad = _stackRadials(rdf)
        # HFRProgs lonlat2km (Bowditch) metres per degree at the grid point's latitude
        glon, glat = grid.x.to_numpy()[:, None], grid.y.to_numpy()[:, None]
        con = np.deg2rad(glat)
        ymetr = 111132.92 - 559.82 * np.cos(2 * con) + 1.175 * np.cos(4 * con) - 0.0023 * np.cos(6 * con)
        xmetr = 111412.84 * np.cos(con) - 93.50 * np.cos(3 * con) + 0.0118 * np.cos(5 * con)
        dist_km = np.hypot((rad["LOND"] - glon) * xmetr, (rad["LATD"] - glat) * ymetr) / 1000
        return (dist_km < radius_km).sum(axis=1)

    @pytest.mark.parametrize("engine", ["loop", "vectorized"])
    def test_flat_distance_selects_hfrprogs_radials(self, combine_grid, engine):
        """distance='flat' counts the radials HFRProgs makeTotals does: lonlat2dist(grid point, radial) < spatthresh."""
        rdf = _radial_dataframe(np.random.default_rng(0))
        tot, _ = combineRadials(rdf, combine_grid, method="uwls", engine=engine, distance="flat", **self.common)
        expected = self._hfrprogs_flat_counts(rdf, combine_grid)

        valid = tot.data["VELU"].notna().to_numpy()
        assert valid.sum() > 20
        npt.assert_array_equal(tot.data.loc[valid, "NRAD"].to_numpy(dtype=float), expected[valid])

    def test_count_all_points_matches_hfrprogs(self, combine_grid):
        """count_all_points=True gives HFRProgs TotalsNumRads at every grid point, also without a total vector."""
        rdf = _radial_dataframe(np.random.default_rng(0))
        strict = dict(minContrSites=2, minContrRads=30)
        tot, _ = combineRadials(rdf, combine_grid, method="uwls", engine="vectorized", distance="flat",
                                count_all_points=True, **strict, **self.common)
        expected = self._hfrprogs_flat_counts(rdf, combine_grid)

        nrad = tot.data["NRAD"].to_numpy(dtype=float)
        no_vector = tot.data["VELU"].isna().to_numpy()
        npt.assert_array_equal(nrad, expected)
        assert (no_vector & (nrad > 0)).sum() > 20  # counts where too few radials gave no vector
        assert (nrad == 0).sum() > 0  # grid points out of range are 0, not NaN

    @pytest.mark.parametrize("method, kwargs", [
        ("uwls", {}),
        ("uwls", dict(distance="flat")),
        ("oi", dict(mdlvar=420, errvar=66, sx=10.0, sy=10.0)),
    ])
    def test_count_all_points_keeps_the_vectors(self, combine_grid, method, kwargs):
        rdf = _radial_dataframe(np.random.default_rng(0))
        default, _ = combineRadials(rdf, combine_grid, method=method, engine="vectorized", **kwargs, **self.common)
        counted, _ = combineRadials(rdf, combine_grid, method=method, engine="vectorized", count_all_points=True,
                                    **kwargs, **self.common)

        others = [c for c in self.columns if c != "NRAD"]
        npt.assert_array_equal(counted.data[others].to_numpy(dtype=float), default.data[others].to_numpy(dtype=float))
        valid = default.data["VELU"].notna()
        assert valid.sum() > 20
        npt.assert_array_equal(counted.data.loc[valid, "NRAD"].to_numpy(dtype=float),
                               default.data.loc[valid, "NRAD"].to_numpy(dtype=float))
        assert counted.data["NRAD"].notna().all()
        assert (counted.data.loc[~valid, "NRAD"] > 0).any()

    def test_count_all_points_without_radials_in_range(self, combine_grid):
        far = combine_grid.translate(xoff=20)  # grid far away from every radial
        tot, _ = combineRadials(_radial_dataframe(), far, method="uwls", engine="vectorized", count_all_points=True,
                                **self.common)
        assert tot.data["VELU"].isna().all()
        assert (tot.data["NRAD"] == 0).all()

    def test_count_all_points_needs_vectorized(self, combine_grid):
        with pytest.raises(ValueError):
            combineRadials(_radial_dataframe(), combine_grid, method="uwls", engine="loop", count_all_points=True,
                           **self.common)

    def test_flat_distance_matches_loop(self, combine_grid):
        rdf = _radial_dataframe(np.random.default_rng(1))
        loop, _ = combineRadials(rdf, combine_grid, method="uwls", engine="loop", distance="flat", **self.common)
        vec, _ = combineRadials(rdf, combine_grid, method="uwls", engine="vectorized", distance="flat", **self.common)
        npt.assert_allclose(vec.data[self.columns].to_numpy(dtype=float), loop.data[self.columns].to_numpy(dtype=float),
                            rtol=1e-6, atol=1e-6)

    def test_invalid_engine_raises(self, combine_grid):
        with pytest.raises(ValueError):
            combineRadials(_radial_dataframe(), combine_grid, method="uwls", engine="fast", **self.common)

    def test_invalid_distance_raises(self, combine_grid):
        with pytest.raises(ValueError):
            combineRadials(_radial_dataframe(), combine_grid, method="uwls", distance="haversine", **self.common)

    def test_vectorized_oi_requires_scalar_errvar(self, combine_grid):
        with pytest.raises(ValueError):
            combineRadials(_radial_dataframe(), combine_grid, method="oi", engine="vectorized",
                           errvar=np.ones(3), sx=10.0, sy=10.0, **self.common)


class TestSiteContribution:
    def test_counts_shares_and_ratio(self):
        # grid point 0: 9 radials from site 0 and 3 from site 1; point 1: 5 from site 0 only;
        # point 2: 4, 4 and 2 radials from three sites; point 3: no radials
        gi = np.array([0] * 12 + [1] * 5 + [2] * 10)
        site = np.array([0] * 9 + [1] * 3 + [0] * 5 + [0] * 4 + [1] * 4 + [2] * 2)
        rng = np.random.default_rng(3)
        order = rng.permutation(gi.size)  # pairs come in any order
        n_sites, ror, share = _siteContribution(gi[order], site[order], nsite=3, npts=4)

        npt.assert_array_equal(n_sites, [2, 1, 3, 0])
        npt.assert_allclose(ror[:3], [3.0, np.inf, 1.0])
        npt.assert_allclose(share[:3], [0.75, 1.0, 0.4])
        assert np.isnan(ror[3]) and np.isnan(share[3])


class TestCombineRadialsWindow:
    """Several radials per site, the time window (tempthreshold) and the site contribution columns."""
    t0 = pd.Timestamp("2025-01-01 12:00")
    common = dict(sRad=10000, gRes=6000, useLandMask=False, dropPoints=False)

    def _window(self, currents):
        """rDF with one radial per site for each (hour offset, uniform current) item."""
        radials, sites = [], []
        for hours, uv in currents:
            for s, (la, lo) in _SITES.items():
                radials.append(_synthetic_radial(s, la, lo, uniform=uv, time=self.t0 + pd.Timedelta(hours=hours)))
                sites.append(s)
        return pd.DataFrame({"Radial": radials}, index=sites)

    def test_repeated_site_codes_count_as_one_site(self, combine_grid):
        single = self._window([(0, (12.0, -7.0))])
        double = self._window([(0, (12.0, -7.0)), (1, (12.0, -7.0))])
        one, _ = combineRadials(single, combine_grid, tStp=self.t0, method="uwls", engine="vectorized", **self.common)
        two, _ = combineRadials(double, combine_grid, tStp=self.t0, method="uwls", engine="vectorized", **self.common)

        valid = one.data["VELU"].notna()
        npt.assert_array_equal(valid, two.data["VELU"].notna())
        npt.assert_array_equal(two.data.loc[valid, "NSIT"], one.data.loc[valid, "NSIT"])
        npt.assert_array_equal(two.data.loc[valid, "NRAD"], 2 * one.data.loc[valid, "NRAD"])
        npt.assert_allclose(two.data.loc[valid, "VELU"].astype(float), 12.0, atol=1e-8)
        assert two.site_source["Name"].tolist() == list(_SITES)

    def test_one_site_at_several_times_is_not_enough(self, combine_grid):
        rdf = self._window([(0, (12.0, -7.0)), (1, (12.0, -7.0))])
        rdf = rdf[rdf.index == "SITA"]
        tot, warn = combineRadials(rdf, combine_grid, tStp=self.t0, method="uwls", engine="vectorized", **self.common)
        assert warn == "No combination performed: not enough contributing radial sites"

    @pytest.mark.parametrize("window", [pd.Timedelta(hours=1), dt.timedelta(hours=1), 1])
    def test_tempthreshold_keeps_radials_within_the_window(self, combine_grid, window):
        # Radials 3 hours away describe a different current and must be left out
        rdf = self._window([(-1, (12.0, -7.0)), (0, (12.0, -7.0)), (1, (12.0, -7.0)), (3, (40.0, 30.0))])
        tot, _ = combineRadials(rdf, combine_grid, tStp=self.t0, method="uwls", engine="vectorized",
                                tempthreshold=window, **self.common)
        ref, _ = combineRadials(self._window([(0, (12.0, -7.0))]), combine_grid, tStp=self.t0, method="uwls",
                                engine="vectorized", **self.common)
        valid = tot.data["VELU"].notna()
        assert valid.sum() > 20
        npt.assert_allclose(tot.data.loc[valid, "VELU"].astype(float), 12.0, atol=1e-8)
        npt.assert_allclose(tot.data.loc[valid, "VELV"].astype(float), -7.0, atol=1e-8)
        npt.assert_array_equal(tot.data.loc[valid, "NRAD"], 3 * ref.data.loc[valid, "NRAD"])

    def test_without_tempthreshold_every_radial_is_used(self, combine_grid):
        rdf = self._window([(0, (12.0, -7.0)), (3, (40.0, 30.0))])
        tot, _ = combineRadials(rdf, combine_grid, tStp=self.t0, method="uwls", engine="vectorized", **self.common)
        valid = tot.data["VELU"].notna()
        assert not np.allclose(tot.data.loc[valid, "VELU"].astype(float), 12.0)

    @pytest.mark.parametrize("method, kwargs", [("uwls", {}), ("oi", dict(mdlvar=625, errvar=100, sx=10.0, sy=10.0))])
    def test_site_contribution_columns(self, combine_grid, method, kwargs):
        rdf = _radial_dataframe(np.random.default_rng(0))
        tot, _ = combineRadials(rdf, combine_grid, tStp=self.t0, method=method, engine="vectorized", **self.common, **kwargs)
        d = tot.data[tot.data["VELU"].notna()]
        assert len(d) > 20
        assert ((d["NSIT"] >= 2) & (d["NSIT"] <= len(_SITES))).all()
        assert (d["ROR"] >= 1).all()
        assert ((d["MSFR"] >= 1 / d["NSIT"]) & (d["MSFR"] < 1)).all()
        assert tot.data.loc[tot.data["VELU"].isna(), ["NSIT", "ROR", "MSFR"]].isna().all().all()

    def test_oi_search_factor_widens_the_search(self, combine_grid):
        rdf = _radial_dataframe(np.random.default_rng(0))
        oi = dict(method="oi", engine="vectorized", mdlvar=625, errvar=100, sx=10.0, sy=10.0)
        narrow, _ = combineRadials(rdf, combine_grid, tStp=self.t0, oi_search_factor=1, **oi, **self.common)
        default, _ = combineRadials(rdf, combine_grid, tStp=self.t0, **oi, **self.common)
        both = narrow.data["VELU"].notna() & default.data["VELU"].notna()
        assert both.sum() > 20
        assert (narrow.data.loc[both, "NRAD"] < default.data.loc[both, "NRAD"]).all()

    OI = dict(method="oi", engine="vectorized", mdlvar=625.0, errvar=100.0, sx=10.0, sy=10.0)

    def _oi(self, rdf, grid, **kwargs):
        tot, _ = combineRadials(rdf, grid, tStp=self.t0, **self.OI, **self.common, **kwargs)
        return tot.data[["VELU", "VELV", "UQAL", "VQAL"]].to_numpy(dtype=float)

    def test_time_scale_without_other_hours_changes_nothing(self, combine_grid):
        rdf = self._window([(0, (12.0, -7.0))])
        npt.assert_array_equal(self._oi(rdf, combine_grid, oi_time_scale=3.0), self._oi(rdf, combine_grid))

    def test_long_time_scale_is_the_unweighted_window(self, combine_grid):
        rdf = self._window([(-1, (12.0, -7.0)), (0, (12.0, -7.0)), (1, (12.0, -7.0))])
        window = dict(tempthreshold=pd.Timedelta(hours=1))
        # Repeated radials at the same bins make the system nearly singular, so allow for rounding
        npt.assert_allclose(self._oi(rdf, combine_grid, oi_time_scale=1e9, **window),
                            self._oi(rdf, combine_grid, **window), rtol=1e-7, atol=1e-7)

    def test_short_time_scale_ignores_other_hours(self, combine_grid):
        # Radials at t-1 and t+1 describe another current; with T -> 0 they are uncorrelated with the totals
        rdf = self._window([(-1, (40.0, 30.0)), (0, (12.0, -7.0)), (1, (40.0, 30.0))])
        only_t = self._oi(self._window([(0, (12.0, -7.0))]), combine_grid)
        npt.assert_allclose(self._oi(rdf, combine_grid, oi_time_scale=1e-6, tempthreshold=pd.Timedelta(hours=1)),
                            only_t, rtol=1e-9, atol=1e-9)

    def test_other_hours_count_less_with_a_shorter_time_scale(self, combine_grid):
        rdf = self._window([(-1, (40.0, 30.0)), (0, (12.0, -7.0)), (1, (40.0, 30.0))])
        window = dict(tempthreshold=pd.Timedelta(hours=1))
        u = {T: self._oi(rdf, combine_grid, oi_time_scale=T, **window)[:, 0] for T in (1.0, 6.0)}
        u[None] = self._oi(rdf, combine_grid, **window)[:, 0]
        ok = np.isfinite(u[1.0]) & np.isfinite(u[6.0]) & np.isfinite(u[None])
        assert ok.sum() > 20
        # The t radials say u = 12 and the others u = 40: a shorter time scale stays closer to 12
        assert np.median(u[1.0][ok]) < np.median(u[6.0][ok]) < np.median(u[None][ok])

    def test_radial_error_column(self, combine_grid):
        rdf = self._window([(0, (12.0, -7.0))])
        for r in rdf["Radial"]:
            r.data["ETMP"] = 5.0
        npt.assert_allclose(self._oi(rdf, combine_grid, oi_errvar_column="ETMP"),
                            combineRadials(rdf, combine_grid, tStp=self.t0, **dict(self.OI, errvar=25.0),
                                           **self.common)[0].data[["VELU", "VELV", "UQAL", "VQAL"]].to_numpy(float),
                            rtol=1e-12, atol=1e-12)
        # Radials without the column use errvar
        for r in rdf["Radial"]:
            r.data = r.data.drop(columns="ETMP")
        npt.assert_array_equal(self._oi(rdf, combine_grid, oi_errvar_column="ETMP"), self._oi(rdf, combine_grid))

    def test_new_oi_options_need_the_vectorized_engine(self, combine_grid):
        rdf = _radial_dataframe()
        for kwargs in (dict(oi_time_scale=3.0), dict(oi_errvar_column="ETMP")):
            with pytest.raises(ValueError):
                combineRadials(rdf, combine_grid, tStp=self.t0, method="oi", engine="loop", **kwargs, **self.common)
        with pytest.raises(ValueError):
            combineRadials(rdf, combine_grid, tStp=self.t0, oi_time_scale=0, **self.OI, **self.common)

    def test_loop_engine_rejects_new_options(self, combine_grid):
        rdf = self._window([(0, (12.0, -7.0)), (1, (12.0, -7.0))])
        with pytest.raises(ValueError):
            combineRadials(rdf, combine_grid, tStp=self.t0, method="uwls", engine="loop", **self.common)
        with pytest.raises(ValueError):
            combineRadials(_radial_dataframe(), combine_grid, tStp=self.t0, method="oi", engine="loop",
                           oi_search_factor=1, **self.common)


class TestQCSiteContribution:
    def _total(self):
        t = _make_total(velo=np.ones(5), velu=np.ones(5), velv=np.zeros(5), head=np.full(5, 90.0))
        t.data["ROR"] = [1.0, 3.0, 5.0, np.inf, np.nan]
        t.data["MSFR"] = [0.5, 0.75, 0.83, 1.0, np.nan]
        return t

    def test_max_share(self):
        t = self._total()
        t.qc_site_contribution(max_share=0.75)
        assert t.data["SCON"].tolist() == [1, 1, 4, 4, 2]

    def test_max_ror(self):
        t = self._total()
        t.qc_site_contribution(max_ror=4)
        assert t.data["SCON"].tolist() == [1, 1, 4, 4, 2]

    def test_metadata_and_primary_flag(self):
        t = self._total()
        t.qc_site_contribution(max_share=0.75)
        assert "qc_site_contribution (SCON)" in t.metadata["QCTest"]["SCON"]
        t.qc_qartod_primary_flag()
        assert t.data["PRIM"].tolist() == [1, 1, 4, 4, 1]

    def test_requires_a_threshold_and_columns(self):
        with pytest.raises(ValueError):
            self._total().qc_site_contribution()
        t = _make_total(velo=np.ones(2), velu=np.ones(2), velv=np.zeros(2), head=np.full(2, 90.0))
        with pytest.raises(ValueError):
            t.qc_site_contribution(max_share=0.75)

    @pytest.mark.parametrize("with_header", [False, True])
    def test_recalculate_primary_flag_includes_new_test(self, with_header):
        t = self._total()
        t.qc_qartod_maximum_velocity()
        t.qc_qartod_primary_flag()
        assert (t.data.PRIM == 1).all()
        if with_header:
            t._tables[1] = {"_TableHeader": [["VELU", "PRIM", "Q303"], ["(cm/s)", "(flag)", "(flag)"]]}
        t.qc_site_contribution(max_share=0.75)
        t.qc_qartod_primary_flag(recalculate=True)
        assert t.data.PRIM.tolist() == [1, 1, 4, 4, 1]
        assert "SCON" in t.metadata["QCTest"]["PRIM"]
        if with_header:
            assert t._tables[1]["_TableHeader"] == [["VELU", "Q303"], ["(cm/s)", "(flag)"]]
        t.qc_qartod_primary_flag(recalculate=True)
        assert t.data.PRIM.tolist() == [1, 1, 4, 4, 1]


class TestQCValidLocationPolygon:
    # Square polygon from -74.8 to -74.2 E and 38.2 to 38.8 N
    polygon = ([-74.8, -74.2, -74.2, -74.8, -74.8], [38.2, 38.2, 38.8, 38.8, 38.2])

    def _total(self):
        # inside, outside, on the edge, outside
        return _make_total(velo=np.ones(4), velu=np.ones(4), velv=np.zeros(4), head=np.full(4, 90.0),
                           lond=[-74.5, -74.0, -74.2, -74.5], latd=[38.5, 38.5, 38.5, 39.0])

    def test_land_polygon_flags_inside_and_edge(self):
        t = self._total()
        t.qc_valid_location_polygon(*self.polygon)
        assert t.data["Q305"].tolist() == [4, 1, 4, 1]

    def test_water_polygon_flags_outside(self):
        t = self._total()
        t.qc_valid_location_polygon(*self.polygon, flag_inside=False)
        assert t.data["Q305"].tolist() == [1, 4, 1, 4]

    def test_empty_polygon_passes_everything(self):
        t = self._total()
        t.qc_valid_location_polygon([], [])
        assert t.data["Q305"].tolist() == [1, 1, 1, 1]

    def test_metadata_and_primary_flag(self):
        t = self._total()
        t.qc_valid_location_polygon(*self.polygon, name="mask.txt")
        assert "qc_valid_location_polygon (Q305)" in t.metadata["QCTest"]["Q305"]
        assert "mask.txt" in t.metadata["QCTest"]["Q305"]
        t.qc_qartod_primary_flag()
        assert t.data["PRIM"].tolist() == [4, 1, 4, 1]

    def test_does_not_run_twice(self):
        t = self._total()
        t.qc_valid_location_polygon(*self.polygon)
        t.qc_valid_location_polygon(*self.polygon, flag_inside=False)
        assert t.data["Q305"].tolist() == [4, 1, 4, 1]


class TestQCThreshold:
    def _total(self):
        t = _make_total(velo=np.ones(4), velu=np.ones(4), velv=np.zeros(4), head=np.full(4, 90.0))
        t.data["CQAL"] = [0.1, -0.6, 0.5, np.nan]
        return t

    def test_maximum(self):
        t = self._total()
        t.qc_threshold("Q309", "CQAL", 0.5)
        assert t.data["Q309"].tolist() == [1, 1, 1, 1]  # 0.5 is not above the threshold; NaN passes

    def test_absolute(self):
        t = self._total()
        t.qc_threshold("Q309", "CQAL", 0.5, absolute=True)
        assert t.data["Q309"].tolist() == [1, 4, 1, 1]

    def test_metadata_and_primary_flag(self):
        t = self._total()
        t.qc_threshold("Q309", "CQAL", 0.5, absolute=True, name="qc_oi_uv_covariance")
        assert "qc_oi_uv_covariance (Q309)" in t.metadata["QCTest"]["Q309"]
        t.qc_qartod_primary_flag()
        assert t.data["PRIM"].tolist() == [1, 4, 1, 1]


class TestHFRProgsGDOP:
    def test_largest_eigenvalue_of_the_covariance(self):
        rng = np.random.default_rng(0)
        a = rng.normal(size=(50, 6, 2))  # 50 grid points, 6 radials each
        cov = np.linalg.inv(np.swapaxes(a, 1, 2) @ a)  # C = inv(A'A)
        t = _make_total(velo=np.ones(50), velu=np.ones(50), velv=np.zeros(50), head=np.full(50, 90.0),
                        gdop=np.sqrt(cov[:, 0, 0] + cov[:, 1, 1]),
                        uqal=np.sqrt(cov[:, 0, 0]), vqal=np.sqrt(cov[:, 1, 1]))
        t.data["CQAL"] = cov[:, 0, 1]
        trace_gdop = t.data["GDOP"].to_numpy()

        t.use_hfrprogs_gdop()
        npt.assert_allclose(t.data["GDOP"], np.sqrt(np.linalg.eigvalsh(cov)[:, -1]), rtol=1e-12)
        assert (t.data["GDOP"].to_numpy() < trace_gdop).all()

    def test_missing_vectors_stay_missing(self):
        t = _make_total(velo=[1.0, np.nan], velu=[1.0, np.nan], velv=[0.0, np.nan], head=[90.0, np.nan],
                        uqal=[1.0, np.nan], vqal=[2.0, np.nan])
        t.data["CQAL"] = [0.0, np.nan]
        t.use_hfrprogs_gdop()
        npt.assert_allclose(t.data["GDOP"].to_numpy(), [2.0, np.nan])


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
