Radial data
===========

The :class:`hfradarpy.radials.Radial` class is the main entry point for CODAR
and WERA radial current files.

Loading files safely
--------------------

::

    from hfradarpy.radials import Radial

    radial = Radial("RDLi_SEAB_2019_01_01_0000.ruv")
    if radial._iscorrupt or not radial.is_valid():
        raise ValueError("The radial could not be parsed")

By default common invalid values such as 999 and 1080 are replaced with NaN.
Set ``replace_invalid=False`` to preserve literal source values. Other
constructor options are:

``mask_over_land=True``
   Remove land vectors during construction. For most workflows it is clearer
   to load first and call ``mask_over_land`` explicitly.

``empty_radial=True``
   Parse the file's structure and metadata, then empty its data tables.

``vflip=True``
   Reverse native CODAR radial velocity so positive means away from the
   radar. See :doc:`concepts` before using this option.

The source table
----------------

``radial.data`` is a DataFrame. Common columns are:

=======  ===============================================================
Column   Meaning
=======  ===============================================================
LOND     Longitude in decimal degrees
LATD     Latitude in decimal degrees
VELO     Signed radial velocity
VELU     Eastward component of radial velocity
VELV     Northward component of radial velocity
BEAR     Bearing clockwise from true north
HEAD     Vector heading
RNGE     Range from the site in kilometres
VFLG     Manufacturer vector/location flag, when present
ETMP     Temporal quality or variability field, when present
ESPC     Spatial quality or variability field, when present
=======  ===============================================================

File variants do not all provide the same fields. Test for a column before
using it::

    if "ETMP" in radial.data:
        good = radial.data[radial.data["ETMP"] < 15]

Metadata remains close to the source header::

    site = radial.metadata.get("Site")
    origin = radial.metadata.get("Origin")
    pattern = radial.metadata.get("PatternType")

Land masking
------------

Return a Boolean water mask without changing the radial::

    over_water = radial.mask_over_land(subset=False, res="low")
    water_data = radial.data.loc[over_water]

Or remove land rows in place::

    radial.mask_over_land(subset=True, res="high")

``res="high"`` uses the 10 m Natural Earth data bundled with the package;
other values select the smaller 110 m dataset. The resolution names describe
the source map scale, not a distance in metres.

xarray models
-------------

Tabular output keeps one row per source vector::

    ds = radial.to_xarray(model="tabular", enhance=True)

Gridded output maps rows to a regular radial grid::

    ds = radial.to_xarray(
        model="gridded",
        range_minmax=(5, 100),
        bearing=(0, 1),
        enhance=True,
    )

For gridded output, ``range_minmax`` sets minimum and maximum range in
kilometres. ``bearing`` specifies antenna bearing and angular resolution in
degrees. When omitted, HFRadarPy infers both grids from metadata and data.

Concatenating files
-------------------

::

    from pathlib import Path
    from hfradarpy.radials import concat

    files = sorted(Path("radials/SEAB").glob("*.ruv"))
    ds = concat(
        files,
        method="gridded",
        range_minmax=(5, 100),
        enhance=True,
        parallel=False,
    )

Files with empty primary tables are skipped and the result is sorted by time.
Every retained file must be compatible with xarray concatenation. For a
mixed collection, group files by site and product type before concatenating.

Writing NetCDF
--------------

::

    radial.to_netcdf(
        "output/radial.nc",
        model="gridded",
        enhance=True,
        user_attributes={
            "title": "Quality-controlled SEAB radial currents",
            "institution": "Example observing network",
        },
    )

When ``user_attributes`` is supplied, HFRadarPy adds required time-coverage
and geospatial attributes. ``prepend_extension=True`` changes the suffix to
``.tabular.nc`` or ``.gridded.nc``. Existing NetCDF destinations are replaced.

Writing CTF radial files
------------------------

After editing the DataFrame or applying QC::

    radial.to_ruv("output/RDLi_SEAB_2019_01_01_0000_qc.ruv")

Validation is enabled by default. If the destination is the original source
and ``overwrite=False``, a modified suffix is used to protect the input. A
round trip preserves the data table, key metadata, and supported QC columns.

Batch-processing pattern
------------------------

::

    from pathlib import Path
    from hfradarpy.radials import Radial

    source = Path("radials/SEAB")
    destination = Path("netcdf/SEAB")

    for path in sorted(source.glob("*.ruv")):
        radial = Radial(path)
        if radial._iscorrupt or not radial.is_valid():
            continue
        radial.to_netcdf(
            destination / path.with_suffix(".nc").name,
            model="gridded",
        )

In operational code, log the exception and source path instead of silently
discarding failed files.

See also
--------

* :doc:`quality_control`
* :doc:`plotting`
* :class:`hfradarpy.radials.Radial`
* :func:`hfradarpy.radials.concat`

