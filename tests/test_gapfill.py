"""Tests for hfradarpy.gapfill (gap filling of gridded HF radar totals, Nwankwo et al. 2023)."""
import numpy as np
import numpy.testing as npt
import pandas as pd
import pytest
import xarray as xr

from hfradarpy.gapfill import (GAUSS_MARKOV, MISSING, NEAREST_NEIGHBOR, OBSERVED, fill_gauss_markov, fill_gaps,
                               fill_nearest_neighbor)


def _dataset(u, v=None, depth=False):
    nt, ny, nx = u.shape
    coords = dict(time=pd.date_range("2024-01-01", periods=nt, freq="h"),
                  lat=28 + 0.054 * np.arange(ny), lon=-95 + 0.061 * np.arange(nx))
    v = -u if v is None else v
    if depth:
        dims = ("time", "depth", "lat", "lon")
        coords["depth"] = [0.0]
        u, v = u[:, None], v[:, None]
    else:
        dims = ("time", "lat", "lon")
    return xr.Dataset({"VELU": (dims, u), "VELV": (dims, v)}, coords=coords)


class TestNearestNeighbor:
    @pytest.mark.parametrize("depth", [False, True])
    def test_hole_in_uniform_field(self, depth):
        u = np.full((3, 5, 5), 12.0)
        u[1, 2, 2] = np.nan
        out = fill_nearest_neighbor(_dataset(u, depth=depth))
        assert out["VELU"].dims == (("time", "depth", "lat", "lon") if depth else ("time", "lat", "lon"))
        filled = out["VELU"].values.reshape(3, 5, 5)
        assert filled[1, 2, 2] == pytest.approx(12.0)
        assert out["fill_flag"].values[1, 2, 2] == NEAREST_NEIGHBOR
        assert (out["fill_flag"].values.sum() == NEAREST_NEIGHBOR)

    def test_average_of_neighbors(self):
        u = np.full((3, 3, 3), np.nan)
        u[0, 0, 0], u[1, 1, 0], u[2, 2, 2] = 1.0, 2.0, 6.0
        out = fill_nearest_neighbor(_dataset(u))
        assert out["VELU"].values[1, 1, 1] == pytest.approx(3.0)
        assert out["VELV"].values[1, 1, 1] == pytest.approx(-3.0)

    def test_too_few_neighbors_stays_missing(self):
        u = np.full((3, 3, 3), np.nan)
        u[0, 0, 0], u[2, 2, 2] = 1.0, 2.0
        out = fill_nearest_neighbor(_dataset(u))
        assert np.isnan(out["VELU"].values[1, 1, 1])
        assert out["fill_flag"].values[1, 1, 1] == MISSING
        assert out["fill_flag"].values[0, 0, 0] == OBSERVED

    def test_cross_neighborhood_ignores_diagonals(self):
        u = np.full((3, 3, 3), np.nan)
        u[0, 0, 0], u[2, 2, 2], u[0, 2, 0] = 1.0, 2.0, 3.0  # all diagonal from the center
        box = fill_nearest_neighbor(_dataset(u))
        cross = fill_nearest_neighbor(_dataset(u), neighborhood="cross")
        assert np.isfinite(box["VELU"].values[1, 1, 1])
        assert np.isnan(cross["VELU"].values[1, 1, 1])

    def test_passes(self):
        u = np.full((1, 1, 7), np.nan)
        u[0, 0, :3] = 5.0
        one = fill_nearest_neighbor(_dataset(u), min_neighbors=1, neighborhood="cross")
        three = fill_nearest_neighbor(_dataset(u), min_neighbors=1, neighborhood="cross", passes=3)
        assert np.isfinite(one["VELU"].values[0, 0]).sum() == 4
        assert np.isfinite(three["VELU"].values[0, 0]).sum() == 6


class TestGaussMarkov:
    def _smooth_field(self, nt=12, ny=12, nx=12):
        t, y, x = np.meshgrid(np.arange(nt), np.arange(ny), np.arange(nx), indexing="ij")
        u = 20 + 5 * np.sin(x / 4) * np.cos(y / 5) + 2 * np.sin(2 * np.pi * t / 12)
        v = -10 + 4 * np.cos(x / 5) + 1.5 * np.cos(2 * np.pi * t / 12)
        return u, v

    def test_recovers_held_out_points(self):
        u, v = self._smooth_field()
        rng = np.random.default_rng(0)
        hold = rng.random(u.shape) < 0.15
        ds = _dataset(np.where(hold, np.nan, u), np.where(hold, np.nan, v))
        out = fill_gauss_markov(ds, length_scale_km={"VELU": 8.0, "VELV": 8.0}, noise=0.01)
        flags = out["fill_flag"].values
        assert (flags[hold] == GAUSS_MARKOV).mean() > 0.95
        assert (flags[~hold] == OBSERVED).all()
        filled = hold & (flags == GAUSS_MARKOV)
        assert np.abs(out["VELU"].values[filled] - u[filled]).max() < 1.0
        assert np.abs(out["VELV"].values[filled] - v[filled]).max() < 1.0

    def test_gaps_without_enough_data_stay_missing(self):
        u = np.full((5, 6, 6), np.nan)
        u[:, 0, 0] = 1.0
        u[2, 0, 0] = np.nan
        out = fill_gauss_markov(_dataset(u), min_points=10)
        assert np.isnan(out["VELU"].values[2, 0, 0])

    def test_never_observed_points_are_not_filled(self):
        u, v = self._smooth_field(6, 6, 6)
        u[:, 3, 3] = np.nan
        v[:, 3, 3] = np.nan
        out = fill_gauss_markov(_dataset(u, v))
        assert np.isnan(out["VELU"].values[:, 3, 3]).all()


class TestFillGaps:
    def test_both_stages(self):
        u = np.full((4, 8, 8), 15.0)
        u[1, 4, 4] = np.nan          # surrounded: nearest neighbor
        u[:, 0:3, 0:3] = np.nan      # never observed, and (0, 0) and (1, 1) have no valid neighbors: stay missing
        u[2:4, 6:8, 6:8] = np.nan    # corner block: needs Gauss-Markov
        nn_only = fill_gaps(_dataset(u))
        both = fill_gaps(_dataset(u), gauss_markov=True)
        assert nn_only["fill_flag"].values[1, 4, 4] == NEAREST_NEIGHBOR
        assert (both["fill_flag"].values == GAUSS_MARKOV).any()
        npt.assert_allclose(both["VELU"].values[both["fill_flag"].values == GAUSS_MARKOV], 15.0, atol=1e-6)
        assert (both["fill_flag"].values[:, 0, 0] == MISSING).all()
        assert (both["fill_flag"].values[:, 1, 1] == MISSING).all()
