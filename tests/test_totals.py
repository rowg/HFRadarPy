"""Tests for hfradarpy.totals"""
import numpy as np
import numpy.testing as npt
import pandas as pd

from hfradarpy.calc import lonlat2km
from hfradarpy.totals import totalOI


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
