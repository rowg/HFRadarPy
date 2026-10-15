Core concepts
=============

Native objects and analysis objects
-----------------------------------

HFRadarPy deliberately supports two views of the same product.

The parsed domain object (``Radial``, ``Total``, or ``Waves``) preserves the
manufacturer-oriented table and header. Its ``data`` attribute is a pandas
``DataFrame`` and its ``metadata`` attribute contains source header fields.
Use this view for quality control, header-aware processing, and CTF round
trips.

The xarray view is created by ``to_xarray``. It is better for labelled array
operations, multi-file analysis, plotting, and NetCDF output. Conversion
returns a dataset; it does not replace ``data``.

Native and enhanced names
-------------------------

With ``enhance=False``, xarray variables retain manufacturer abbreviations.
With ``enhance=True``, known radial fields are renamed, including:

================  =============================
Native field      Enhanced radial name
================  =============================
``VELO``          ``velocity``
``VELU``          ``u``
``VELV``          ``v``
``LOND``          ``lon``
``LATD``          ``lat``
``BEAR``          ``bearing``
``RNGE``          ``range``
``VFLG``          ``vector_flag``
``ETMP``          ``temporal_quality``
``ESPC``          ``spatial_quality``
================  =============================

Only names that exist in the source are renamed. Print the dataset or inspect
``dataset.data_vars`` instead of assuming every format contains every field.

Tabular and gridded radial models
---------------------------------

``model="tabular"``
   Preserves one observation per source row. This is the safest choice for
   round-trip checks, irregular products, and DataFrame-like filtering.

``model="gridded"``
   Places radial observations onto ``(time, range, bearing)`` coordinates and
   provides two-dimensional longitude and latitude coordinates. Missing
   range/bearing cells become NaN. This model is convenient for concatenation
   and map-based analysis.

A total-vector gridded dataset normally uses ``(time, depth, lat, lon)``.
Functions in :mod:`hfradarpy.gapfill` and :mod:`hfradarpy.validation` accept a
singleton ``depth`` dimension where relevant.

Units and signs
---------------

Native CODAR radial velocity fields are centimetres per second. Their radial
sign convention is positive toward the radar. ``Radial(..., vflip=True)``
reverses ``MINV``, ``MAXV``, and ``VELO`` in the parsed table so positive is
away from the radar.

Do not infer units from a variable name. Inspect ``DataArray.attrs["units"]``
after enhancement or consult the input format when working with native tables.
The validation functions preserve the units supplied to them, except
``drifter_velocity``, which returns metres per second.

Mutation and copies
-------------------

QC methods and ``mask_over_land(subset=True)`` mutate the domain object's
``data`` table. ``to_xarray`` and ``to_netcdf`` work from a deep copy so the
conversion does not reshape the original table. Gap-filling functions return
a modified copy of the input xarray dataset.

Missing data
------------

Parsers replace common manufacturer fill values with NaN by default. Pass
``replace_invalid=False`` only when the literal source values are required.
NetCDF writers encode missing floating-point values as ``-999`` on disk; an
xarray read normally decodes the fill value back to NaN.

Time
----

Product timestamps are interpreted from file headers or standard LLUV file
names. Operational workflows should use UTC consistently. The parsed CTF
``TimeZone`` value remains available in ``metadata`` for provenance.
