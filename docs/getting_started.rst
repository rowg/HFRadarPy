Getting started
===============

This walkthrough uses a CODAR radial file, but the same object model applies
throughout HFRadarPy: parse a file into a domain object, work with its pandas
table and metadata, then convert to xarray or write a product.

The repository includes sample files under ``examples/data``. From the
repository root::

    from pathlib import Path

    from hfradarpy.radials import Radial

    path = Path("examples/data/radials/ruv/SEAB/")
    radial = Radial(path / "RDLi_SEAB_2019_01_01_0100.ruv")

Inspect a parsed file
---------------------

The most useful attributes are:

``data``
   The primary data table as a ``pandas.DataFrame``. Column names initially
   follow the source format (for example ``VELO``, ``LOND``, and ``LATD``).

``metadata``
   An ordered mapping of CTF header fields. Values generally remain strings
   because the exact source header must be preserved for round trips.

``time``
   The product timestamp as a Python ``datetime``.

``file_name`` and ``full_file``
   The source name and path.

``is_wera`` and ``_iscorrupt``
   Parser state describing the manufacturer and whether parsing failed.

For example::

    print(radial)
    print(radial.time)
    print(radial.metadata["Origin"])
    print(radial.data.columns.tolist())
    print(radial.data[["LOND", "LATD", "VELO", "BEAR", "RNGE"]].head())

Check ``radial.is_valid()`` before starting a batch conversion. Empty or
corrupt inputs can otherwise fail later when a required coordinate or header
field is accessed.

Convert to xarray
-----------------

Use the tabular model when every source row should remain an observation::

    tabular = radial.to_xarray(model="tabular", enhance=True)
    print(tabular.sizes)

Use the gridded model for range-by-bearing analysis::

    gridded = radial.to_xarray(model="gridded", enhance=True)
    print(gridded["velocity"].dims)
    print(gridded["velocity"].attrs["units"])

``enhance=True`` replaces known manufacturer abbreviations with descriptive
names and adds metadata. Use ``enhance=False`` when downstream code requires
the native column names. See :doc:`concepts` for the implications.

Quality control
---------------

Quality-control methods update ``radial.data`` in place. Initialize the QC
metadata once, run the individual tests, then calculate the aggregate primary
flag last::

    radial.initialize_qc()
    radial.qc_qartod_syntax()
    radial.qc_qartod_maximum_velocity(max_speed=250, high_speed=150)
    radial.qc_qartod_valid_location(use_mask=True, res="low")
    radial.qc_qartod_spatial_median()
    radial.qc_qartod_primary_flag()

    print(radial.data[["Q201", "Q202", "Q203", "Q205", "PRIM"]].head())

Flag values follow the IOOS convention: 1 pass, 2 not evaluated, 3 suspect,
4 fail, and 9 missing. :doc:`quality_control` documents every test and its
output column.

Write output
------------

Write compressed NetCDF::

    radial.to_netcdf(
        "examples/output/SEAB_2019_01_01_0100.nc",
        model="gridded",
        enhance=True,
    )

Or preserve the CTF representation, including QC columns::

    radial.to_ruv(
        "examples/output/RDLi_SEAB_2019_01_01_0100_qc.ruv",
        validate=True,
    )

HFRadarPy creates parent directories for these exports. NetCDF export replaces
an existing destination. CTF export protects the original source name unless
``overwrite=True`` is passed.

Work with a time series
-----------------------

``hfradarpy.radials.concat`` accepts paths, ``Radial`` objects, or a mixture::

    from hfradarpy.radials import concat

    files = sorted(path.glob("*.ruv"))
    series = concat(files, method="gridded", enhance=True, parallel=False)
    print(series.time.min().item(), series.time.max().item())

Use ``parallel=True`` for a large collection after validating a small serial
sample. Parsing uses all available CPU cores and may raise memory use.

Next steps
----------

* :doc:`radials` explains selection, concatenation, units, and round trips.
* :doc:`quality_control` gives a production-ready QC sequence.
* :doc:`totals` covers total files and radial combination.
* :doc:`waves_diagnostics` covers the other supported CODAR products.
* :doc:`analysis` covers gap filling and in-situ validation.

