
Installation
============

HFRadarPy supports Python 3.10 and newer. Geospatial packages such as Cartopy,
GeoPandas, GEOS, and PROJ have compiled dependencies, so conda-forge is the
most reliable installation route.

Conda-forge (recommended)
-------------------------

Create an isolated environment and install the release::

    mamba create -n hfradarpy -c conda-forge python=3.11 hfradarpy
    mamba activate hfradarpy

``conda`` may be used in place of ``mamba``. Verify the installation::

    python -c "import hfradarpy; print(hfradarpy.__version__)"

PyPI
----

If the required geospatial libraries are already available on your system::

    python -m pip install hfradarpy

Use a virtual environment; avoid installing into the operating system's
Python installation.

Install from source
-------------------

For development or unreleased changes::

    git clone https://github.com/rowg/HFRadarPy.git
    cd HFRadarPy
    mamba env create -f environment_dev.yml
    mamba activate hfradarpy
    python -m pip install -e .

Run the tests from the repository root::

    pytest

Build the documentation
-----------------------

Install the documentation dependencies, then ask Sphinx to treat warnings as
errors::

    python -m pip install -r docs/requirements.txt
    python -m sphinx -W --keep-going -b html docs docs/_build/html

Open ``docs/_build/html/index.html`` in a browser after a successful build.

Troubleshooting
---------------

``ImportError`` for Cartopy, GeoPandas, GEOS, or PROJ
   Install HFRadarPy from conda-forge in a clean environment. Mixing packages
   from several conda channels or mixing system and conda geospatial
   libraries commonly causes binary-compatibility failures.

Natural Earth data or land masking fails
   ``mask_over_land`` and valid-location tests depend on the land-mask assets
   installed with the package and on GeoPandas. Confirm that ``hfradarpy`` was
   installed as a package rather than copying an individual Python file.

NetCDF writing fails
   HFRadarPy writes with the ``netcdf4`` engine. Confirm that ``netCDF4`` is
   importable and that the destination directory is writable.

Large concatenations use too much memory
   Start with a short time range, use gridded subsets via ``range_minmax`` and
   ``bearing``, or process files in batches. ``parallel=True`` speeds parsing
   but can increase peak memory use.
