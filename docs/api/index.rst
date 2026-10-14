API reference
=============

The API reference is organized by task. Public classes and functions are
listed with their live signatures and docstrings. Start with a user guide if
you need a workflow rather than a symbol lookup.

.. toctree::
   :maxdepth: 2

   radials
   totals
   waves
   gapfill
   validation
   plotting
   io
   utilities
   diagnostics_patterns

Stability
---------

Names beginning with an underscore are internal. The low-level structures in
:mod:`hfradarpy.hfrnet` primarily support totals computation and may be less
stable than the high-level ``Radial``, ``Total``, and ``Waves`` interfaces.

