Analysis: gap filling and validation
====================================

HFRadarPy includes focused tools for filling gaps in regular total-vector
grids and comparing HFR currents with drifters or fixed instruments.

Gap filling
-----------

Gap-filling functions expect an xarray dataset with ``time``, ``lat``, and
``lon`` dimensions. A singleton ``depth`` dimension is supported. The default
variables are the native total fields ``VELU`` and ``VELV``.

Nearest-neighbour averaging
~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    from hfradarpy.gapfill import fill_nearest_neighbor

    filled = fill_nearest_neighbor(
        totals,
        variables=("VELU", "VELV"),
        min_neighbors=3,
        passes=1,
        neighborhood="box",
    )

The two components are filled together. A neighbour counts only when all
requested variables are finite there.

``neighborhood="box"`` uses the 26 adjacent cells in a 3 x 3 x 3
time/latitude/longitude cube. ``"cross"`` uses only the six cells one step
along an axis. Additional passes allow newly filled cells to support later
fills and therefore expand the affected area.

Gauss-Markov interpolation
~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    from hfradarpy.gapfill import fill_gauss_markov

    filled = fill_gauss_markov(
        totals,
        length_scale_km={"VELU": 3.6, "VELV": 4.72},
        time_scale_hours=4,
        max_distance_km=18,
        max_dt_hours=8,
        min_points=10,
        max_points=30,
    )

The algorithm estimates anomalies relative to each grid point's time mean.
A location that is never observed has no mean and remains missing. Keep
``min_points`` large enough for a stable local estimate and validate the
length and time scales for the region being studied.

Two-stage filling
~~~~~~~~~~~~~~~~~

::

    from hfradarpy.gapfill import fill_gaps

    filled = fill_gaps(
        totals,
        min_neighbors=3,
        passes=1,
        gauss_markov=True,
        max_distance_km=18,
        max_dt_hours=8,
    )

``fill_gaps`` applies nearest-neighbour averaging first and optionally applies
Gauss-Markov interpolation to the remaining gaps. It returns a copy and adds
``fill_flag``:

=====  =========================
Value  Meaning
=====  =========================
0      Observed
1      Nearest-neighbour average
2      Gauss-Markov interpolation
9      Missing
=====  =========================

Keep this flag with derived products; filled values should not be presented
as observations.

Drifter velocity
----------------

Convert an ordered drifter position series to eastward and northward velocity::

    from hfradarpy.validation import drifter_velocity

    u_drifter, v_drifter = drifter_velocity(lon, lat, time)

Central differences use positions ``i-1`` and ``i+1`` on the configured
ellipsoid. The first and last result are NaN. Output units are metres per
second.

Harmonize in-situ data
----------------------

Remove unrealistic speeds and apply a centred moving average::

    from hfradarpy.validation import harmonize

    observed = harmonize(
        u_drifter,
        v_drifter,
        time,
        window="3h",
        max_speed=2.5,
    )

``max_speed`` uses the input velocity units. Pass ``window=None`` to skip
averaging.

Sample the HFR grid
-------------------

Collocate the nearest grid point and time step to each in-situ observation::

    from hfradarpy.validation import sample_totals

    u_hfr, v_hfr = sample_totals(
        totals,
        lon=drifter_lon,
        lat=drifter_lat,
        time=drifter_time,
        max_distance_km=6,
        max_dt="30min",
        u="VELU",
        v="VELV",
    )

Selections beyond either limit become NaN. Set ``u`` and ``v`` to ``"u"``
and ``"v"`` for datasets using enhanced variable names.

Validation statistics
---------------------

::

    from hfradarpy.validation import (
        complex_correlation,
        skill_score,
        variance_ellipse,
    )

    correlation, veering = complex_correlation(
        u_observed, v_observed, u_hfr, v_hfr
    )
    skill, skill_u, skill_v = skill_score(
        u_observed, v_observed, u_hfr, v_hfr
    )
    ellipse = variance_ellipse(u_hfr, v_hfr)

``complex_correlation`` returns correlation magnitude and veering angle in
degrees. A positive angle means the second series is rotated counterclockwise
from the first. ``demean=True`` compares anomalies.

``skill_score`` returns the mean component score and its u/v components. Its
default denominator uses the estimate mean following Fang et al. (2015); use
``mean_of="reference"`` for the original Willmott form.

``variance_ellipse`` returns mean velocity, major and minor standard-deviation
axes, major-axis angle counterclockwise from east, and the number of paired
samples.

Unit consistency
----------------

These comparison statistics do not convert units. Convert drifter velocities
from m/s when totals are in cm/s, or convert totals to m/s::

    u_drifter_cm_s = 100 * u_drifter
    v_drifter_cm_s = 100 * v_drifter

Apply the same convention to speed filters and reported metrics.

