"""
Validation of HF radar total vectors against in situ currents (drifters and buoys).

The functions follow the validation of Nwankwo et al. (2026), "Review of the Texas shelf currents using high
frequency radar surface currents data", Ocean Dynamics 76:57:
    * drifter velocities from central differences of drifter positions (their eq. 1)
    * harmonization of the time series (speed limit and a 3-hour moving average)
    * complex correlation and veering angle (Kundu 1976)
    * skill score (Fang et al. 2015, eqs. 15-17)
    * variance ellipses

References:
    Fang, Y.-C., Weingartner, T. J., Potter, R. A., Winsor, P. R. and Statscewich, H., 2015. Quality assessment of
        HF radar-derived surface currents using optimal interpolation. J. Atmos. Oceanic Technol., 32, 282-296.
    Kundu, P. K., 1976. Ekman veering observed near the ocean bottom. J. Phys. Oceanogr., 6, 238-242.
    Willmott, C. J., 1981. On the validation of models. Phys. Geogr., 2, 184-194.
"""
import numpy as np
import pandas as pd
import xarray as xr
from pyproj import Geod

from hfradarpy.calc import lonlat2km


def drifter_velocity(lon, lat, time, ellps='GRS80'):
    """
    Velocity of a drifter from its positions, by central differences (Nwankwo et al. 2026, eq. 1).

    The velocity at position i is the displacement from position i-1 to position i+1, measured on the
    ellipsoid, divided by the time between them. The first and last positions get NaN.

    INPUT:
        lon, lat: positions of one drifter (degrees), in time order
        time: times of the positions (anything pandas.to_datetime accepts)
        ellps: ellipsoid for the distances and bearings (default 'GRS80', as in Nwankwo et al. 2026)

    OUTPUT:
        u, v: eastward and northward velocity (m/s), same length as the input
    """
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    time = pd.to_datetime(np.asarray(time))
    seconds = np.asarray((time - time[0]) / pd.Timedelta(seconds=1), dtype=float)
    u = np.full(lon.size, np.nan)
    v = np.full(lon.size, np.nan)
    if lon.size < 3:
        return u, v

    az, _, dist = Geod(ellps=ellps).inv(lon[:-2], lat[:-2], lon[2:], lat[2:])
    dt = seconds[2:] - seconds[:-2]
    with np.errstate(divide='ignore', invalid='ignore'):
        u[1:-1] = dist * np.sin(np.deg2rad(az)) / dt
        v[1:-1] = dist * np.cos(np.deg2rad(az)) / dt
    return u, v


def harmonize(u, v, time, window='3h', max_speed=None):
    """
    Prepare an in situ current time series for comparison with HF radar (Nwankwo et al. 2026, section 3.5).

    Velocities faster than max_speed are removed first, then a centered moving average of length window is
    applied, matching the 3-hour averaging of the radials.

    INPUT:
        u, v: velocity components
        time: times of the samples
        window: length of the centered moving average (pandas offset string, default '3h'); None skips it
        max_speed: remove samples faster than this (same units as u and v); None keeps every sample

    OUTPUT:
        pandas.DataFrame indexed by time with columns u and v
    """
    df = pd.DataFrame({'u': np.asarray(u, dtype=float), 'v': np.asarray(v, dtype=float)},
                      index=pd.to_datetime(np.asarray(time)))
    df = df.sort_index()
    if max_speed is not None:
        df[np.hypot(df['u'], df['v']) > max_speed] = np.nan
    if window is not None:
        df = df.rolling(window, center=True, min_periods=1).mean()
    return df


def sample_totals(ds, lon, lat, time, max_distance_km=None, max_dt=None, u='VELU', v='VELV'):
    """
    HF radar velocities at the positions and times of in situ measurements (nearest grid point and time step).

    INPUT:
        ds: gridded totals (xarray.Dataset with time, lat and lon coordinates), e.g. Total.to_xarray(model='gridded')
            concatenated over time
        lon, lat, time: positions and times of the in situ samples
        max_distance_km: samples farther than this from the nearest grid point get NaN (default: no limit)
        max_dt: samples farther than this from the nearest time step get NaN (pandas.Timedelta or offset string;
                default: no limit)
        u, v: names of the velocity variables in ds (default 'VELU', 'VELV')

    OUTPUT:
        u, v: numpy arrays of the HF radar velocities at the samples
    """
    lon = np.atleast_1d(np.asarray(lon, dtype=float))
    lat = np.atleast_1d(np.asarray(lat, dtype=float))
    time = pd.to_datetime(np.atleast_1d(np.asarray(time)))

    # Pointwise selection: one (lon, lat, time) per sample
    indexers = {'lon': xr.DataArray(lon, dims='sample'), 'lat': xr.DataArray(lat, dims='sample'),
                'time': xr.DataArray(time.values, dims='sample')}
    picked = ds[[u, v]].sel(indexers, method='nearest')
    # Drop any other dimension (e.g. depth) by taking its first level
    picked = picked.isel({d: 0 for d in picked.dims if d != 'sample'})
    hu = picked[u].values.astype(float)
    hv = picked[v].values.astype(float)

    bad = np.zeros(lon.size, dtype=bool)
    if max_distance_km is not None:
        dx, dy = lonlat2km(lon, lat, picked['lon'].values, picked['lat'].values)
        bad |= np.hypot(dx, dy) > max_distance_km
    if max_dt is not None:
        bad |= np.abs(pd.to_datetime(picked['time'].values) - time) > pd.Timedelta(max_dt)
    hu[bad] = np.nan
    hv[bad] = np.nan
    return hu, hv


def _paired(*arrays):
    """Arrays as floats, keeping only the samples that are finite in all of them."""
    arrays = [np.asarray(a, dtype=float).ravel() for a in arrays]
    ok = np.logical_and.reduce([np.isfinite(a) for a in arrays])
    return [a[ok] for a in arrays]


def complex_correlation(u1, v1, u2, v2, demean=False):
    """
    Complex correlation between two vector time series (Kundu 1976).

    With w = u + iv, r = <conj(w1) w2> / (<|w1|^2>^1/2 <|w2|^2>^1/2), so
        Re(r) = <u1 u2 + v1 v2> / norm,  Im(r) = <u1 v2 - v1 u2> / norm,
    and the veering angle is the argument of r: atan(<u1 v2 - v1 u2> / <u1 u2 + v1 v2>) (Kundu 1976;
    Nwankwo et al. 2023, eq. B1; Fang et al. 2015, eq. 18). A positive angle means series 2 is rotated
    counterclockwise from series 1. Samples missing in either series are ignored.

    INPUT:
        u1, v1: reference series (e.g. buoy or drifter)
        u2, v2: compared series (e.g. HF radar)
        demean: remove each series' mean before correlating (default False, as in the references)

    OUTPUT:
        r: magnitude of the complex correlation (0-1)
        angle: veering angle in degrees
    """
    u1, v1, u2, v2 = _paired(u1, v1, u2, v2)
    if u1.size == 0:
        return np.nan, np.nan
    w1, w2 = u1 + 1j * v1, u2 + 1j * v2
    if demean:
        w1, w2 = w1 - w1.mean(), w2 - w2.mean()
    rho = np.mean(np.conj(w1) * w2) / np.sqrt(np.mean(np.abs(w1) ** 2) * np.mean(np.abs(w2) ** 2))
    return float(np.abs(rho)), float(np.rad2deg(np.angle(rho)))


def skill_score(u_ref, v_ref, u_est, v_est, mean_of='estimate'):
    """
    Skill of estimated currents against reference currents (Fang et al. 2015, eqs. 15-17).

        Skill = (SU + SV) / 2
        SU = 1 - sum |u_ref - u_est|^2 / sum (|u_ref - mean| + |u_est - mean|)^2

    and likewise SV for v. A skill of 1 is perfect agreement and 0 no agreement. Fang et al. (2015) use
    the mean of the estimated series (their subscript o) in the denominator; mean_of='reference' gives the
    original index of agreement of Willmott (1981). Samples missing in any series are ignored.

    INPUT:
        u_ref, v_ref: reference currents (e.g. buoy or drifter; "m" in Fang et al. 2015)
        u_est, v_est: estimated currents (e.g. HF radar; "o" in Fang et al. 2015)
        mean_of: 'estimate' (default, Fang et al. 2015) or 'reference' (Willmott 1981)

    OUTPUT:
        skill, su, sv
    """
    if mean_of not in ('estimate', 'reference'):
        raise ValueError("mean_of must be 'estimate' or 'reference'")
    u_ref, v_ref, u_est, v_est = _paired(u_ref, v_ref, u_est, v_est)
    if u_ref.size == 0:
        return np.nan, np.nan, np.nan

    def component(ref, est):
        mean = (est if mean_of == 'estimate' else ref).mean()
        denominator = np.sum((np.abs(ref - mean) + np.abs(est - mean)) ** 2)
        if denominator == 0:
            return 1.0 if np.allclose(ref, est) else 0.0
        return 1 - np.sum((ref - est) ** 2) / denominator

    su, sv = component(u_ref, u_est), component(v_ref, v_est)
    return (su + sv) / 2, su, sv


def variance_ellipse(u, v):
    """
    Variance (principal axis) ellipse of a vector time series.

    INPUT:
        u, v: velocity components; missing samples are ignored

    OUTPUT:
        dict with
            u_mean, v_mean: mean velocity (the ellipse center)
            major, minor: semi-axes, the standard deviations along the principal axes
            angle: direction of the major axis in degrees counterclockwise from east, between -90 and 90
            n: number of samples used
    """
    u, v = _paired(u, v)
    out = dict(u_mean=np.nan, v_mean=np.nan, major=np.nan, minor=np.nan, angle=np.nan, n=int(u.size))
    if u.size < 2:
        return out
    cov = np.cov(np.vstack((u, v)), ddof=0)
    eigval, eigvec = np.linalg.eigh(cov)
    major_vec = eigvec[:, 1]
    angle = np.rad2deg(np.arctan2(major_vec[1], major_vec[0]))
    if angle > 90:
        angle -= 180
    elif angle <= -90:
        angle += 180
    out.update(u_mean=u.mean(), v_mean=v.mean(), major=np.sqrt(max(eigval[1], 0)),
               minor=np.sqrt(max(eigval[0], 0)), angle=angle)
    return out
