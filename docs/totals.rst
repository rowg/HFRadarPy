Total-vector data
=================

HFRadarPy can parse existing total-vector files with
:class:`hfradarpy.totals.Total` or create totals from radial measurements with
:func:`hfradarpy.totals.combineRadials`.

Read an existing total
----------------------

::

    from hfradarpy.totals import Total

    total = Total("hourly_total.tuv")
    print(total.time)
    print(total.data[["LOND", "LATD", "VELU", "VELV", "GDOP"]].head())

The parser accepts CODAR LLUV totals and WERA CUR totals. As with radial
objects, ``data`` is a DataFrame and ``metadata`` preserves source header
fields.

Convert to a regular xarray grid::

    ds = total.to_xarray(model="gridded")
    print(ds["VELU"].dims)  # time, depth, lat, lon

Supply a two-column longitude/latitude grid file or a GeoSeries to retain
unobserved cells in a prescribed grid::

    ds = total.to_xarray(grid="examples/data/grid_files/totals.grid.txt")

``Total`` currently exposes conversion to xarray rather than a dedicated
NetCDF writer. Use xarray to write the result::

    ds.to_netcdf("output/total.nc")

Create an empty total on a grid
-------------------------------

::

    import geopandas as gpd
    from shapely.geometry import Point
    from hfradarpy.totals import Total

    grid = gpd.GeoSeries(
        [Point(-75.0, 38.0), Point(-74.5, 38.0)],
        crs="EPSG:4326",
    )
    total = Total(grid=grid)

This initializes the standard total columns and coordinate metadata. It is
also how ``combineRadials`` creates its output object.

Combine radials
---------------

Inputs are a DataFrame indexed by site code with a ``Radial`` column, a
GeoSeries of grid points, distances in metres, and the target timestamp::

    import datetime as dt
    import geopandas as gpd
    import pandas as pd
    from shapely.geometry import Point

    from hfradarpy.radials import Radial
    from hfradarpy.totals import combineRadials

    radials = pd.DataFrame({
        "Radial": pd.Series({
            "SITE1": Radial("RDLm_SITE1_2025_01_01_0000.ruv"),
            "SITE2": Radial("RDLm_SITE2_2025_01_01_0000.ruv"),
        })
    })
    grid = gpd.GeoSeries(
        [Point(-74.5, 39.0), Point(-74.4, 39.0)],
        crs="EPSG:4326",
    )

    total, warning = combineRadials(
        radials,
        grid,
        sRad=10_000,
        gRes=6_000,
        tStp=dt.datetime(2025, 1, 1),
        method="uwls",
        engine="vectorized",
        minContrSites=2,
        minContrRads=3,
        useLandMask=False,
        dropPoints=False,
    )
    if warning:
        raise RuntimeError(warning)

The function always returns ``(total, warning)``. A non-empty warning means
the requested combination was not completed as expected.

Combination methods
-------------------

``method="wls"``
   Weighted least squares. The radial ``ETMP`` field supplies the standard
   deviation used for inverse-variance weights.

``method="uwls"``
   Unweighted least squares. Use when valid, comparable uncertainty estimates
   are unavailable.

``method="oi"``
   Optimal interpolation. Configure model variance ``mdlvar``, error variance
   ``errvar``, east-west and north-south decorrelation lengths ``sx`` and
   ``sy``, and ``oi_option`` (``"exponential"`` or ``"gaussian"``).

For OI, verify the units expected by your processing configuration. The
example workflows use velocity variances in ``(cm/s)^2`` and spatial scales
in kilometres.

Engines
-------

``engine="vectorized"`` uses a spatial tree and array operations and is the
recommended engine for new work. It also:

* enforces minimum site and radial counts for every method;
* permits several time steps from the same site when ``tempthreshold`` is set;
* creates ``NSIT``, ``ROR``, and ``MSFR`` site-contribution diagnostics;
* supports ``count_all_points=True``; and
* supports alternate ``oi_search_factor`` values.

``engine="loop"`` retains the original grid-point implementation. It cannot
accept duplicate site indices, custom OI search factors, or
``count_all_points=True``.

Spatial selection
-----------------

For least squares, ``sRad`` is the radial search radius in metres.
``distance="geodesic"`` measures on the Total CRS ellipsoid;
``distance="flat"`` follows the HFRProgs flat-earth selection rule.

For OI, candidates satisfy an ellipse based on ``sx`` and ``sy`` and bounded
by ``oi_search_factor`` (default 2). ``sRad`` is still recorded as total
metadata but is not the OI candidate rule.

Time-window combinations
------------------------

The vectorized engine can combine neighbouring radial times::

    total, warning = combineRadials(
        radial_rows,
        grid,
        10_000,
        6_000,
        target_time,
        method="oi",
        engine="vectorized",
        tempthreshold=pd.Timedelta("1h"),
        mdlvar=625,
        errvar=100,
        sx=10,
        sy=10,
    )

The DataFrame may repeat a site code. All rows from that site within the time
window count as one contributing site.

Output fields
-------------

The main fields include eastward and northward velocity (``VELU``, ``VELV``),
speed and heading (``VELO``, ``HEAD``), uncertainty/covariance terms
(``UQAL``, ``VQAL``, ``CQAL``), geometry (``GDOP``), and contributing radial
count (``NRAD``). Vectorized output also includes the site diagnostics above.

Set ``dropPoints=False`` to preserve the full requested grid. With
``count_all_points=True``, ``NRAD`` is populated even at cells where minimum
contribution rules prevented a total vector.

Quality control and export
--------------------------

Apply total QC before conversion::

    total.initialize_qc()
    total.qc_qartod_maximum_velocity(max_speed=150, high_speed=100)
    total.qc_qartod_gdop(max_GDOP=2)
    total.qc_qartod_primary_flag()

    ds = total.to_xarray(grid=grid)
    ds.to_netcdf("output/qc_total.nc")

See :doc:`quality_control` for all total-vector tests.

