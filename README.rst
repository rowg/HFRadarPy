=========
HFRadarPy
=========

HFRadarPy is a Python toolbox for reading, quality controlling, transforming,
combining, and exporting high-frequency radar (HFR) ocean-surface-current
data. It understands CODAR tabular (CTF/LLUV) and WERA products and exposes
their data as familiar pandas and xarray objects.

HFRadarPy is useful when you need to:

* inspect CODAR radial, total, wave, diagnostic, or antenna-pattern files;
* read WERA radial and total files;
* run IOOS QARTOD quality-control tests and retain their flags;
* convert radial, total, and wave products to compressed NetCDF;
* concatenate a time series of files into an xarray dataset;
* combine radial measurements into total vectors with least squares or
  optimal interpolation; or
* fill gaps and validate total-vector products against in-situ observations.

Quick start
-----------

Install the latest release from conda-forge::

    mamba install -c conda-forge hfradarpy

Load a CODAR radial, inspect its native table, run QARTOD tests, and create a
CF-oriented NetCDF file::

    from pathlib import Path

    from hfradarpy.radials import Radial

    source = Path("RDLi_SEAB_2019_01_01_0100.ruv")
    radial = Radial(source)

    print(radial.time)
    print(radial.metadata["Site"])
    print(radial.data[["LOND", "LATD", "VELO"]].head())

    radial.initialize_qc()
    radial.qc_qartod_syntax()
    radial.qc_qartod_maximum_velocity(max_speed=250, high_speed=150)
    radial.qc_qartod_spatial_median()
    radial.qc_qartod_primary_flag()

    radial.to_netcdf("output/SEAB_2019_01_01_0100.nc", model="gridded")

For in-memory analysis, convert without writing a file::

    dataset = radial.to_xarray(model="gridded", enhance=True)
    speed = abs(dataset["velocity"])

The ``Radial`` object keeps the manufacturer's table in ``radial.data`` as a
``pandas.DataFrame``. ``to_xarray`` returns a new ``xarray.Dataset`` and does
not modify that table.

Documentation
-------------

The `online user guide <https://hfradarpy.readthedocs.io/>`_ includes
practical documentation for radial data, quality control, totals, waves,
plotting, gap filling, validation, the CLI, and the complete public API.

The repository also contains `runnable examples
<https://github.com/rowg/HFRadarPy/tree/master/examples>`_ and sample input
files for CODAR and WERA products.

Supported data
--------------

.. list-table::
   :header-rows: 1
   :widths: 25 45 30

   * - Product
     - Common inputs
     - Primary interface
   * - Radial currents
     - CODAR ``.ruv``; WERA radial
     - ``Radial``
   * - Total currents
     - CODAR ``.tuv``; WERA ``.cur``
     - ``Total``
   * - Wave products
     - CODAR ``.wls``
     - ``Waves``
   * - Diagnostics
     - CODAR ``.rdt/.hdt/.xdt``
     - ``Diags``
   * - Antenna patterns
     - CODAR ``.loop``
     - ``Loop``

Velocity convention
-------------------

Native CODAR radial velocities are stored in centimetres per second and are
positive toward the radar. HFRadarPy preserves those source values in
``Radial.data`` unless ``vflip=True`` is requested. Enhanced xarray output
adds descriptive names and CF-style metadata; always inspect each variable's
``units`` attribute before combining it with another data source.

Project links
-------------

* Source: https://github.com/rowg/HFRadarPy
* Issues: https://github.com/rowg/HFRadarPy/issues
* Documentation: https://hfradarpy.readthedocs.io/
* License: MIT

Citation
--------

If HFRadarPy contributes to published work, cite the software using the
metadata in ``CITATION.cff`` and record the version used::

    import hfradarpy
    print(hfradarpy.__version__)
