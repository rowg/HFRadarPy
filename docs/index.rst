HFRadarPy documentation
=======================

HFRadarPy turns manufacturer HFR files into analysis-ready pandas and xarray
objects. It covers the path from raw radial measurements through QARTOD
quality control, total-vector production, NetCDF export, plotting,
gap-filling, and validation.

If you are new to the package, begin with :doc:`getting_started`. If you
already know the workflow you need, use the task guides below. The
:doc:`api/index` is the exhaustive lookup reference.

Choose a workflow
-----------------

:doc:`radials`
   Parse CODAR and WERA radial products, understand the data model, and
   produce xarray or NetCDF output.

:doc:`quality_control`
   Apply IOOS QARTOD tests, interpret flags, create a primary flag, and write
   QC results back to CTF or NetCDF.

:doc:`totals`
   Load total files or combine radial measurements with WLS, UWLS, and
   optimal interpolation.

:doc:`analysis`
   Fill gaps, collocate observations, compute validation statistics, and make
   maps.

.. toctree::
   :maxdepth: 2
   :caption: Start here

   installation
   getting_started
   concepts

.. toctree::
   :maxdepth: 2
   :caption: User guides

   radials
   quality_control
   totals
   waves_diagnostics
   analysis
   plotting
   command_line
   examples

.. toctree::
   :maxdepth: 2
   :caption: Reference

   api/index
   contributing
   authors
   history

Getting help
------------

Search these docs first, then open an issue on the `HFRadarPy issue tracker
<https://github.com/rowg/HFRadarPy/issues>`_. Include the package version,
input file type, a minimal example, and the full traceback. Do not attach
sensitive operational data to a public issue.

Indices
-------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
