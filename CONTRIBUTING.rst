============
Contributing
============

Contributions to HFRadarPy are welcome. Bug reports, focused feature
proposals, documentation improvements, new file-format fixtures, and tested
code changes all help the project.

Before opening an issue
-----------------------

Search the `issue tracker <https://github.com/rowg/HFRadarPy/issues>`_ for an
existing report. For a bug, include:

* the HFRadarPy and Python versions;
* operating system and installation method;
* manufacturer and product type (for example CODAR radial or WERA total);
* a minimal code example;
* the complete traceback; and
* a small, shareable input file when licensing and operational policy permit.

Do not publish credentials, private endpoints, or restricted observational
data. A feature proposal should explain the use case, expected interface, and
how success can be tested.

Development setup
-----------------

Fork the repository, then clone your fork::

    git clone https://github.com/YOUR-USER/HFRadarPy.git
    cd HFRadarPy

Create the development environment and install the checkout in editable
mode::

    mamba env create -f environment_dev.yml
    mamba activate hfradarpy
    python -m pip install -e .

Create a focused branch::

    git switch -c fix/short-description

Running tests
-------------

Run the full suite from the repository root::

    pytest

During development, run the narrowest relevant file first::

    pytest tests/test_radials.py
    pytest tests/test_qc.py -k spatial_median

Changes to parsers or exporters should normally include a round-trip test.
Algorithm changes should cover normal input, missing values, boundary
conditions, and invalid configuration. Keep new fixtures as small as possible.

Documentation
-------------

User-visible changes require documentation. Put task-oriented explanation in
the appropriate guide under ``docs/`` and keep public docstrings accurate.
Examples must be runnable and should state units, coordinate conventions,
mutation behavior, and return types when those are not obvious.

Build the documentation with warnings treated as errors::

    python -m pip install -r docs/requirements.txt
    python -m sphinx -W --keep-going -b html docs docs/_build/html

The documentation is organized as:

* ``getting_started.rst`` for a first successful workflow;
* task guides such as ``radials.rst`` and ``totals.rst``;
* ``concepts.rst`` for conventions shared across workflows; and
* ``api/`` for generated symbol reference.

Pull requests
-------------

Keep each pull request scoped to one problem. A strong pull request:

* explains the user-visible problem and the chosen behavior;
* includes tests that fail before and pass after the change;
* updates documentation and ``HISTORY.rst`` when appropriate;
* avoids unrelated formatting or generated-file changes; and
* passes the test and strict documentation builds.

If a behavior changes scientific output, document the previous and new
results, units, tolerances, and validation method. Preserve source-file
provenance in any comparison.

Code style
----------

Follow the surrounding code, use descriptive names, and keep public APIs
backward compatible unless a breaking change is intentional and documented.
Prefer ``pathlib.Path`` for new path handling. Public functions and methods
need docstrings that describe parameters, return values, units, mutation, and
raised exceptions.

Releases
--------

Releases are performed by maintainers. Update the history entry before
tagging, run the full tests and documentation build, and verify package
metadata and installation from the built artifact.
