"""
Gap filling of gridded HF radar total vectors.

Implements the two-step gap filling of Nwankwo et al. (2023), section 2.1 and appendix A, which Nwankwo et al.
(2026) applied to the Texas shelf HF radar fields:

1. Nearest-neighbor averaging: a missing value is replaced by the mean of the valid values one grid step away
   in space and time, when at least three of them exist.
2. Gauss-Markov interpolation (optional): remaining gaps are estimated from nearby data anomalies,
   u' = R_ud R_dd^-1 d', with the space-time correlation R(r, tau) = 0.5 exp(-r^2/delta^2) exp(-tau^2/T^2) + 0.5.
   Gaps without enough nearby data are left empty.

The fields are an xarray.Dataset with time, lat and lon dimensions, for example the output of
Total.to_xarray(model='gridded') concatenated over time. Other dimensions of length one (such as depth) are allowed.

References:
    Nwankwo, U., Howden, S., Nechaev, D. and Dzwonkowski, B., 2023. Subinertial sea surface heights anomalies
        estimated using high frequency radar surface current data in the Mississippi Bight. J. Geophys. Res. Oceans,
        128, e2022JC019055.
    Nwankwo, U., Nechaev, D., DiMarco, S., Knap, A. and Forbes, C., 2026. Review of the Texas shelf currents using
        high frequency radar surface currents data. Ocean Dynamics, 76, 57.
"""
import warnings

import numpy as np
import pandas as pd
import xarray as xr
from scipy import ndimage
from scipy.spatial import cKDTree

from hfradarpy.calc import lonlat2km

# Values of the fill_flag variable
OBSERVED = 0
NEAREST_NEIGHBOR = 1
GAUSS_MARKOV = 2
MISSING = 9


def _as_tll(da):
    """Array of a variable ordered (time, lat, lon), with other length-one dimensions squeezed out."""
    extra = [d for d in da.dims if d not in ('time', 'lat', 'lon')]
    for d in extra:
        if da.sizes[d] != 1:
            raise ValueError(f"Dimension '{d}' must have length 1 for gap filling")
    return da.squeeze(extra, drop=True).transpose('time', 'lat', 'lon').values.astype(float)


def _put_back(ds, name, values):
    """Write a (time, lat, lon) array back into ds[name] with its original dimensions."""
    da = ds[name]
    extra = [d for d in da.dims if d not in ('time', 'lat', 'lon')]
    template = da.squeeze(extra, drop=True).transpose('time', 'lat', 'lon')
    ds[name] = template.copy(data=values).expand_dims({d: da[d].values for d in extra}).transpose(*da.dims)


def _flags(ds, variables):
    """Existing fill_flag (time, lat, lon), or a new one from the observed values."""
    if 'fill_flag' in ds:
        return _as_tll(ds['fill_flag']).astype(np.int8)
    valid = np.logical_and.reduce([np.isfinite(_as_tll(ds[v])) for v in variables])
    return np.where(valid, OBSERVED, MISSING).astype(np.int8)


def _store_flags(ds, flags):
    da = xr.DataArray(flags, dims=('time', 'lat', 'lon'), coords={d: ds[d] for d in ('time', 'lat', 'lon')})
    da.attrs = dict(long_name='Gap filling flag', flag_values=np.array([OBSERVED, NEAREST_NEIGHBOR, GAUSS_MARKOV, MISSING],
                                                                          dtype=np.int8),
                    flag_meanings='observed nearest_neighbor_average gauss_markov missing')
    ds['fill_flag'] = da


def fill_nearest_neighbor(ds, variables=('VELU', 'VELV'), min_neighbors=3, passes=1, neighborhood='box'):
    """
    Fill gaps with the average of the valid neighbors one grid step away in space and time
    (Nwankwo et al. 2023, section 2.1).

    The variables are filled together: a neighbor counts only when all variables are valid there.

    INPUT:
        ds: xarray.Dataset with time, lat and lon dimensions
        variables: names of the variables to fill (default the eastward and northward velocities)
        min_neighbors: minimum number of valid neighbors needed to fill a gap (default 3)
        passes: number of times to repeat the fill; values filled in one pass can fill gaps in the next (default 1)
        neighborhood: 'box' (default) uses the 26 neighbors of the 3x3x3 time-lat-lon box; 'cross' uses only the
                      6 neighbors one step along a single axis

    OUTPUT:
        a filled copy of ds with a fill_flag variable (0 observed, 1 nearest-neighbor average, 9 missing)
    """
    if neighborhood == 'box':
        kernel = np.ones((3, 3, 3))
        kernel[1, 1, 1] = 0
    elif neighborhood == 'cross':
        kernel = ndimage.generate_binary_structure(3, 1).astype(float)
        kernel[1, 1, 1] = 0
    else:
        raise ValueError("neighborhood must be 'box' or 'cross'")

    ds = ds.copy()
    variables = list(variables)
    values = {v: _as_tll(ds[v]) for v in variables}
    flags = _flags(ds, variables)

    for _ in range(passes):
        valid = np.logical_and.reduce([np.isfinite(values[v]) for v in variables])
        count = ndimage.convolve(valid.astype(float), kernel, mode='constant', cval=0.0)
        fill = ~valid & (count >= min_neighbors)
        if not fill.any():
            break
        for v in variables:
            total = ndimage.convolve(np.where(valid, values[v], 0.0), kernel, mode='constant', cval=0.0)
            values[v] = np.where(fill, total / np.maximum(count, 1), values[v])
        flags[fill] = NEAREST_NEIGHBOR

    for v in variables:
        _put_back(ds, v, values[v])
    _store_flags(ds, flags)
    return ds


def fill_gauss_markov(ds, variables=('VELU', 'VELV'), length_scale_km=None, time_scale_hours=4.0, amplitude=0.5,
                      offset=0.5, max_distance_km=18.0, max_dt_hours=8.0, max_points=30, min_points=10,
                      noise=0.1, batch_size=5000):
    """
    Fill gaps by Gauss-Markov interpolation of nearby data (Nwankwo et al. 2023, appendix A).

    For each gap, u' = R_ud R_dd^-1 d', where d' are the anomalies of the nearby data from the time-mean field,
    R_dd their correlation matrix and R_ud their correlation with the gap. The correlation of two points r km and
    tau hours apart is amplitude * exp(-r^2 / delta^2) * exp(-tau^2 / T^2) + offset. The estimate u' is added to
    the time mean at the gap.

    Nwankwo et al. (2023) used delta = 3.6 km for u and 4.72 km for v, T = 4 hours, amplitude = offset = 0.5.
    They do not give the size of the data window, the minimum number of data or an error term, so these are
    settings here: the nearest max_points valid values within max_distance_km and max_dt_hours are used, gaps
    with fewer than min_points of them are left empty, and noise is added to the diagonal of R_dd (as a fraction
    of the variance) to keep it invertible.

    INPUT:
        ds: xarray.Dataset with time, lat and lon dimensions (e.g. after fill_nearest_neighbor)
        variables: names of the variables to fill
        length_scale_km: dict {variable: delta in km}; default {first variable: 3.6, second variable: 4.72}
        time_scale_hours: T in hours (default 4)
        amplitude, offset: correlation model parameters (default 0.5 and 0.5)
        max_distance_km, max_dt_hours: size of the data window around each gap (default 18 km and 8 hours)
        max_points: largest number of data used for one gap (default 30)
        min_points: gaps with fewer data in the window are not filled (default 10)
        noise: error variance added to the diagonal of R_dd, as a fraction of the variance (default 0.1)
        batch_size: number of gaps solved together (limits memory use)

    OUTPUT:
        a filled copy of ds with a fill_flag variable (0 observed, 1 nearest-neighbor average, 2 Gauss-Markov,
        9 missing)
    """
    variables = list(variables)
    if length_scale_km is None:
        length_scale_km = dict(zip(variables, [3.6, 4.72]))
    missing_scales = [v for v in variables if v not in length_scale_km]
    if missing_scales:
        raise ValueError(f"No length scale for {', '.join(missing_scales)}")

    ds = ds.copy()
    values = {v: _as_tll(ds[v]) for v in variables}
    flags = _flags(ds, variables)

    # Positions in km (east, north of the domain center) and time in hours
    lat = ds['lat'].values.astype(float)
    lon = ds['lon'].values.astype(float)
    hours = (pd.to_datetime(ds['time'].values) - pd.to_datetime(ds['time'].values[0])) / pd.Timedelta(hours=1)
    hours = np.asarray(hours, dtype=float)
    lon2, lat2 = np.meshgrid(lon, lat)
    x, y = lonlat2km(lon.mean(), lat.mean(), lon2, lat2)
    nt, ny, nx = flags.shape
    X = np.broadcast_to(x, (nt, ny, nx)).ravel()
    Y = np.broadcast_to(y, (nt, ny, nx)).ravel()
    T = np.broadcast_to(hours[:, None, None], (nt, ny, nx)).ravel()

    valid = np.logical_and.reduce([np.isfinite(values[v]) for v in variables]).ravel()
    with warnings.catch_warnings():
        # Grid points that are never observed have no mean; they are not filled
        warnings.simplefilter('ignore', category=RuntimeWarning)
        means = {v: np.nanmean(values[v], axis=0) for v in variables} if valid.any() else {}
    has_mean = np.logical_and.reduce([np.isfinite(means[v]) for v in variables]) if means else np.zeros((ny, nx), bool)
    gaps = np.flatnonzero(~valid & np.broadcast_to(has_mean, (nt, ny, nx)).ravel())
    data = np.flatnonzero(valid)
    if gaps.size == 0 or data.size < min_points:
        for v in variables:
            _put_back(ds, v, values[v])
        _store_flags(ds, flags)
        return ds

    # Nearest data inside the space-time window (an ellipsoid in normalized coordinates)
    def scaled(idx):
        return np.column_stack((X[idx] / max_distance_km, Y[idx] / max_distance_km, T[idx] / max_dt_hours))

    tree = cKDTree(scaled(data))
    k = min(max_points, data.size)
    _, nearest = tree.query(scaled(gaps), k=k, distance_upper_bound=1.0)
    nearest = nearest.reshape(gaps.size, k)
    present = nearest < data.size
    enough = present.sum(axis=1) >= min_points
    gaps, nearest, present = gaps[enough], nearest[enough], present[enough]
    nb = data[np.where(present, nearest, 0)]  # flat indices of the neighbor data

    anomalies = {v: (values[v].ravel() - np.broadcast_to(means[v], (nt, ny, nx)).ravel()) for v in variables}
    filled = {v: np.full(gaps.size, np.nan) for v in variables}

    for start in range(0, gaps.size, batch_size):
        sl = slice(start, start + batch_size)
        g, n, m = gaps[sl], nb[sl], present[sl]
        # Gap-data and data-data separations
        r_ud = np.hypot(X[n] - X[g][:, None], Y[n] - Y[g][:, None])
        t_ud = T[n] - T[g][:, None]
        r_dd = np.hypot(X[n][:, :, None] - X[n][:, None, :], Y[n][:, :, None] - Y[n][:, None, :])
        t_dd = T[n][:, :, None] - T[n][:, None, :]
        m2 = m[:, :, None] & m[:, None, :]
        eye = np.eye(n.shape[1])[None, :, :]
        for v in variables:
            delta = length_scale_km[v]
            R_ud = np.where(m, amplitude * np.exp(-(r_ud / delta) ** 2) * np.exp(-(t_ud / time_scale_hours) ** 2)
                            + offset, 0.0)
            R_dd = amplitude * np.exp(-(r_dd / delta) ** 2) * np.exp(-(t_dd / time_scale_hours) ** 2) + offset
            # Padded neighbors form an identity block that does not touch the real data
            R_dd = np.where(m2, R_dd, 0.0) + eye * np.where(m, noise, 1.0)[:, :, None]
            d = np.where(m, anomalies[v][n], 0.0)
            weights = np.linalg.solve(R_dd, R_ud[:, :, None])[:, :, 0]
            filled[v][sl] = np.einsum('bk,bk->b', weights, d)

    gap_t, gap_y, gap_x = np.unravel_index(gaps, (nt, ny, nx))
    for v in variables:
        out = values[v].copy()
        out[gap_t, gap_y, gap_x] = filled[v] + means[v][gap_y, gap_x]
        values[v] = out
        _put_back(ds, v, out)
    flags[gap_t, gap_y, gap_x] = GAUSS_MARKOV
    _store_flags(ds, flags)
    return ds


def fill_gaps(ds, variables=('VELU', 'VELV'), min_neighbors=3, passes=1, neighborhood='box', gauss_markov=False,
              **gauss_markov_options):
    """
    Gap filling of Nwankwo et al. (2023): nearest-neighbor averaging, then optionally Gauss-Markov interpolation.

    INPUT:
        ds: xarray.Dataset with time, lat and lon dimensions
        variables: names of the variables to fill
        min_neighbors, passes, neighborhood: see fill_nearest_neighbor
        gauss_markov: also fill the remaining gaps with fill_gauss_markov (default False)
        gauss_markov_options: options passed to fill_gauss_markov

    OUTPUT:
        a filled copy of ds with a fill_flag variable (0 observed, 1 nearest-neighbor average, 2 Gauss-Markov,
        9 missing)
    """
    ds = fill_nearest_neighbor(ds, variables, min_neighbors=min_neighbors, passes=passes, neighborhood=neighborhood)
    if gauss_markov:
        ds = fill_gauss_markov(ds, variables, **gauss_markov_options)
    return ds
