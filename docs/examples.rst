Examples and recipes
====================

The ``examples`` directory contains runnable programs, sample input files,
and notebooks. Run scripts from the repository root after installing the
package in the active environment.

Core conversions
----------------

``radials_to_netcdf.py``
   Convert the included radial files to tabular and gridded NetCDF.

``wave_to_netcdf.py``
   Apply a wave-height range check and convert CODAR wave files to NetCDF.

``plot_timeseries_from_radials.ipynb``
   Explore a radial time series interactively.

``radial_class_exploration.ipynb``
   Inspect the ``Radial`` object, metadata, and tabular fields.

``plot_radials_and_totals.ipynb``
   Map radial and total-vector products.

Total-vector workflows
----------------------

``radials_to_totals_oi.py``
   End-to-end optimal interpolation with a configurable time window,
   total-vector QC, optional daily averaging, and optional gap filling. Run::

       python examples/radials_to_totals_oi.py --help

``radials_to_totals_benchmark.py``
   Compare loop and vectorized engines on archived radial data.

``synthetic_totals_benchmark.py``
   Generate a known current field in memory, simulate radial observations,
   and assess WLS/UWLS/OI accuracy and performance.

``realtime_totals.py``
   A larger operational example with configuration files, input change
   detection, QC, NetCDF encoding, and optional MongoDB bookkeeping. Treat it
   as an application built with HFRadarPy, not as part of the stable library
   API.

Validation study
----------------

``maracoos_study.py`` and ``maracoos_study_analysis.py`` demonstrate a
checkpointed, reproducible comparison of optimal-interpolation settings and
validation metrics. ``maracoos_totals_study.ipynb`` provides the notebook
version.

Choose an example
-----------------

===============================  =============================================
Goal                             Start with
===============================  =============================================
Read one radial                  :doc:`getting_started`
Convert a radial directory       ``examples/radials_to_netcdf.py``
QC a radial                      :doc:`quality_control`
Build hourly total vectors       ``examples/radials_to_totals_oi.py``
Compare engine performance       ``examples/synthetic_totals_benchmark.py``
Fill a gridded time series       :doc:`analysis`
Validate against a drifter       :doc:`analysis`
Build an operational processor   ``examples/realtime_totals.py``
===============================  =============================================

Example data are intentionally small enough for tests and demonstrations.
They are not reference observations for scientific interpretation.

