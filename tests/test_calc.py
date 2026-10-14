"""Tests for hfradarpy.calc utility functions."""
import datetime as dt

import numpy as np
import numpy.testing as npt
import pytest
from pyproj import Geod

from hfradarpy.calc import (
    dms2dd,
    find_nearest,
    inverse_transformation,
    reckon,
    roundToNearest20,
    scircle1,
    spdir2uv,
    true2mathAngle,
    uv2spdir,
)


class TestRoundToNearest20:
    def test_round_down_to_0(self):
        result = roundToNearest20(dt.datetime(2024, 1, 1, 12, 5, 0))
        assert result == dt.datetime(2024, 1, 1, 12, 0, 0)

    def test_round_to_20(self):
        result = roundToNearest20(dt.datetime(2024, 1, 1, 12, 15, 0))
        assert result == dt.datetime(2024, 1, 1, 12, 20, 0)

    def test_round_to_40(self):
        result = roundToNearest20(dt.datetime(2024, 1, 1, 12, 35, 0))
        assert result == dt.datetime(2024, 1, 1, 12, 40, 0)

    def test_round_up_to_next_hour(self):
        result = roundToNearest20(dt.datetime(2024, 1, 1, 12, 55, 0))
        assert result == dt.datetime(2024, 1, 1, 13, 0, 0)

    def test_string_input(self):
        result = roundToNearest20("2024-01-01 12:05:00")
        assert result == dt.datetime(2024, 1, 1, 12, 0, 0)

    def test_exact_boundary(self):
        result = roundToNearest20(dt.datetime(2024, 1, 1, 12, 10, 0))
        assert result == dt.datetime(2024, 1, 1, 12, 0, 0)


class TestFindNearest:
    def test_exact_match(self):
        arr = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        val, idx = find_nearest(arr, 3.0)
        assert val == 3.0
        assert idx == 2

    def test_between_values(self):
        arr = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        val, idx = find_nearest(arr, 2.3)
        assert val == 2.0
        assert idx == 1

    def test_negative_values(self):
        arr = np.array([-5.0, -2.0, 0.0, 2.0, 5.0])
        val, idx = find_nearest(arr, -1.0)
        assert val == -2.0
        assert idx == 1


class TestDms2dd:
    def test_basic_conversion(self):
        # 40 degrees, 30 minutes, 0 seconds = 40.5 degrees
        result = dms2dd((40, 30, 0))
        assert result == pytest.approx(40.5)

    def test_with_seconds(self):
        # 40 degrees, 30 minutes, 36 seconds = 40.51 degrees
        result = dms2dd((40, 30, 36))
        assert result == pytest.approx(40.51)

    def test_zero(self):
        result = dms2dd((0, 0, 0))
        assert result == 0.0


class TestTrue2MathAngle:
    def test_north(self):
        # North (0 degrees true) = 90 degrees math
        result = true2mathAngle(np.array([0.0]))
        npt.assert_array_almost_equal(result, [90.0])

    def test_east(self):
        # East (90 degrees true) = 0 degrees math
        result = true2mathAngle(np.array([90.0]))
        npt.assert_array_almost_equal(result, [0.0])

    def test_south(self):
        # South (180 degrees true) = 270 degrees math
        result = true2mathAngle(np.array([180.0]))
        npt.assert_array_almost_equal(result, [270.0])

    def test_west(self):
        # West (270 degrees true) = 180 degrees math
        result = true2mathAngle(np.array([270.0]))
        npt.assert_array_almost_equal(result, [180.0])

    def test_radians_input(self):
        # pi/2 radians = 90 degrees true = 0 degrees math
        result = true2mathAngle(np.array([np.pi / 2]), radians=True)
        npt.assert_array_almost_equal(result, [0.0])


class TestSpdir2uv:
    def test_northward_flow(self):
        # Speed=1, direction=0 (north in math convention = radians)
        u, v = spdir2uv(np.array([1.0]), np.array([0.0]))
        npt.assert_array_almost_equal(u, [0.0])
        npt.assert_array_almost_equal(v, [1.0])

    def test_eastward_flow(self):
        # Speed=1, direction=pi/2 (east)
        u, v = spdir2uv(np.array([1.0]), np.array([90.0]), degrees=True)
        npt.assert_array_almost_equal(u, [1.0])
        npt.assert_array_almost_equal(v, [0.0], decimal=5)

    def test_zero_speed(self):
        u, v = spdir2uv(np.array([0.0]), np.array([45.0]), degrees=True)
        npt.assert_array_almost_equal(u, [0.0])
        npt.assert_array_almost_equal(v, [0.0])


class TestUv2spdir:
    def test_pure_eastward(self):
        speed, direction = uv2spdir(np.array([1.0]), np.array([0.0]))
        npt.assert_array_almost_equal(speed, [1.0])
        npt.assert_array_almost_equal(direction, [90.0])

    def test_pure_northward(self):
        speed, direction = uv2spdir(np.array([0.0]), np.array([1.0]))
        npt.assert_array_almost_equal(speed, [1.0])
        npt.assert_array_almost_equal(direction, [0.0])

    def test_zero_velocity(self):
        speed, direction = uv2spdir(np.array([0.0]), np.array([0.0]))
        npt.assert_array_almost_equal(speed, [0.0])


class TestReckon:
    def test_northward_displacement(self):
        # Starting at 0,0, heading north (0 degrees), 111.32 km ~ 1 degree latitude
        lon, lat = reckon(0.0, 0.0, 0.0, 111.32)
        assert lat == pytest.approx(1.0, abs=0.01)
        assert lon == pytest.approx(0.0, abs=0.01)

    def test_eastward_displacement(self):
        # Starting at 0,0, heading east (90 degrees), ~111.32 km
        lon, lat = reckon(0.0, 0.0, 90.0, 111.32)
        assert lon == pytest.approx(1.0, abs=0.01)
        assert lat == pytest.approx(0.0, abs=0.01)


class TestInverseTransformation:
    def test_known_distance(self):
        # Two points 1 degree of latitude apart at the equator
        fwd_az, back_az, dist = inverse_transformation(
            np.array([0.0]), np.array([0.0]),
            np.array([0.0]), np.array([1.0])
        )
        # Distance should be ~110.57 km (WGS84 at equator)
        assert dist[0] == pytest.approx(110.57, abs=0.5)
        # Forward azimuth should be ~0 (due north)
        assert fwd_az[0] == pytest.approx(0.0, abs=0.01)


class TestScircle1:
    # Reference vertices produced by the previous polycircles-based implementation
    def test_matches_polycircles_reference(self):
        lat, lon = scircle1(41.2504500, -69.9740667, 210.0)
        assert len(lat) == len(lon) == 100
        expected = {
            0: (43.14102876259626, -69.9740667),
            25: (41.19322187855507, -67.47072842159889),
            50: (39.360174318748605, -70.05136117002404),
            75: (41.31317820714869, -72.47946459324108),
            99: (43.14102876259626, -69.9740667),
        }
        for i, (exp_lat, exp_lon) in expected.items():
            assert lat[i] == pytest.approx(exp_lat, abs=1e-9)
            assert lon[i] == pytest.approx(exp_lon, abs=1e-9)

    def test_custom_vertex_count_matches_reference(self):
        lat, lon = scircle1(36.6917167, -75.9226333, 10.0, num=36)
        assert len(lat) == 37
        expected = {
            0: (36.78182894007643, -75.9226333),
            9: (36.69166412390019, -75.81073864673829),
            18: (36.60160309293067, -75.9226333),
            27: (36.69166412390019, -76.03452795326172),
        }
        for i, (exp_lat, exp_lon) in expected.items():
            assert lat[i] == pytest.approx(exp_lat, abs=1e-9)
            assert lon[i] == pytest.approx(exp_lon, abs=1e-9)

    def test_vertices_lie_on_the_circle(self):
        lat0, lon0, radius_km = 38.5, -74.0, 75.0
        lat, lon = scircle1(lat0, lon0, radius_km)
        assert lat[0] == lat[-1] and lon[0] == lon[-1]  # closed polygon
        az, _, dist = Geod(ellps="WGS84").inv(np.full(lat.size, lon0), np.full(lat.size, lat0), lon, lat)
        npt.assert_allclose(dist, radius_km * 1e3, rtol=1e-9)
        npt.assert_allclose(np.mod(az[:-1], 360), 360 / 99 * np.arange(99), atol=1e-6)
