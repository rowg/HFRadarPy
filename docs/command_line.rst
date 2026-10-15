Command-line interface
======================

Installing the package creates the ``hfradarpy`` command. Show its commands::

    hfradarpy --help

The command prints the installed HFRadarPy version before running a
subcommand.

Concatenate radial files
------------------------

``concat-radials`` finds ``*.ruv`` files in one directory, optionally selects
a time interval from standard filenames, creates a gridded enhanced dataset,
and writes NetCDF::

    hfradarpy concat-radials \
      --path_data /data/radials/SEAB \
      --path_save /data/netcdf \
      --time_start "2019-01-01 00:00:00" \
      --time_end "2019-01-02 00:00:00" \
      --fname seab \
      --parallel

The output filename includes ``_radials_`` and the actual first and last
times. Use ``hfradarpy concat-radials --help`` for the option list.

Extract radial time series
--------------------------

``extract-timeseries`` concatenates gridded radials and selects the nearest
bearing/range cell to one or more longitude/latitude points::

    hfradarpy extract-timeseries \
      --path_data /data/radials/SEAB \
      --path_save /data/timeseries \
      --point=-74.25,39.10 \
      --point=-74.00,38.80 \
      --time_start "2019-01-01 00:00:00" \
      --time_end "2019-01-02 00:00:00" \
      --type netcdf \
      --fname seab

``--point`` is repeatable and must be ``longitude,latitude``. The selection is
nearest-neighbour in bearing and range from the receiving site; it is not a
spatial interpolation. Output type is ``csv`` or ``netcdf``.

Plot radial files
-----------------

Plot one radial file::

    hfradarpy plotruv \
      --radial_file /data/radials/SEAB/RDLi_SEAB_2019_01_01_0000.ruv \
      --path_save /data/figures \
      --speed_display color \
      --limits -100 100

Or plot every matching file in a directory and time interval::

    hfradarpy plotruv \
      --path_data /data/radials/SEAB \
      --path_save /data/figures \
      --prefix RDLi \
      --time_start "2019-01-01 00:00:00" \
      --time_end "2019-01-02 00:00:00" \
      --plotflag PRIM

When plotting a directory, ``--prefix`` is required because it forms the file
pattern. ``--plotflag`` highlights vectors whose named QC flag is greater than
1. See :doc:`plotting` for the Python plotting interface.

Operational notes
-----------------

* Times should be supplied in UTC.
* Input files must use the standard ``YYYY_MM_DD_HHMM`` filename timestamp.
* A directory with no matching files cannot be concatenated.
* ``--parallel`` uses all CPUs and may increase peak memory.
* The CLI does not run quality control; use the Python API for configurable
  QC workflows.

The repository's scripts under ``examples`` expose additional purpose-built
``argparse`` interfaces. They are examples rather than installed CLI
commands; run each script with ``--help``.
