Totals API
==========

High-level interface
--------------------

.. currentmodule:: hfradarpy.totals

.. autosummary::
   :nosignatures:

   Total
   combineRadials
   combineRadials_hfrnet

.. automodule:: hfradarpy.totals
   :members: Total, combineRadials, combineRadials_hfrnet
   :show-inheritance:

Combination primitives
----------------------

These functions are useful for algorithm-level work. Most applications
should use :func:`combineRadials`.

.. autofunction:: totalLeastSquare
.. autofunction:: totalUnWeightedLeastSquare
.. autofunction:: totalOI
.. autofunction:: radBinsInSearchRadius
.. autofunction:: makeTotalVector
.. autofunction:: makeTotalVector_uwls
.. autofunction:: makeTotalVector_oi
.. autofunction:: makeTotalVectors_ls_vectorized
.. autofunction:: makeTotalVectors_oi_vectorized

HFRNet-compatible structures
----------------------------

.. automodule:: hfradarpy.hfrnet
   :members:
   :show-inheritance:

