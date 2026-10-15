Quality control
===============

HFRadarPy implements IOOS QARTOD-style tests for radial and total currents,
plus European HFR Node (EHN) tests. This guide focuses on the QARTOD workflow.

Flag meanings
-------------

=====  =================
Flag   Meaning
=====  =================
1      Pass
2      Not evaluated
3      Suspect
4      Fail
9      Missing data
=====  =================

Tests add columns to ``data`` and descriptions to ``metadata["QCTest"]``.
They do not automatically delete suspect or failed observations.

Radial tests
------------

Call ``initialize_qc`` before the first test. A typical sequence is::

    from hfradarpy.radials import Radial

    current = Radial("RDLi_SEAB_2019_01_01_0100.ruv", replace_invalid=False)
    previous = "RDLi_SEAB_2019_01_01_0000.ruv"

    current.initialize_qc()
    current.qc_qartod_syntax()                              # Q201
    current.qc_qartod_maximum_velocity(250, 150)            # Q202
    current.qc_qartod_valid_location(use_mask=True)         # Q203
    current.qc_qartod_radial_count(150, 300)                # Q204
    current.qc_qartod_spatial_median()                      # Q205
    current.qc_qartod_temporal_gradient(previous)           # Q206
    current.qc_qartod_avg_radial_bearing(180)               # Q207
    current.qc_qartod_stuck_value(resolution=0.01, N=3)     # Q209
    current.qc_qartod_primary_flag()                        # PRIM

=================================  ======  ========================================
Method                             Column  Purpose
=================================  ======  ========================================
``qc_qartod_syntax``               Q201    File structure and required fields
``qc_qartod_maximum_velocity``     Q202    Gross range and suspect-speed thresholds
``qc_qartod_valid_location``       Q203    Vendor flags, land mask, angular segments
``qc_qartod_radial_count``         Q204    File-wide count thresholds
``qc_qartod_spatial_median``       Q205    Difference from neighbouring velocities
``qc_qartod_temporal_gradient``    Q206    Change from the previous radial
``qc_qartod_avg_radial_bearing``   Q207    Mean bearing relative to a reference
``qc_qartod_stuck_value``          Q209    Repeated values over prior time steps
``qc_operator_test``               Q211    Operator-wide or segment override
``qc_qartod_primary_flag``         PRIM    Worst flag across selected tests
=================================  ======  ========================================

The maximum-velocity and temporal-gradient thresholds use the same units as
the native ``VELO`` column—normally cm/s for CODAR radial files.

Valid-location segments
~~~~~~~~~~~~~~~~~~~~~~~

Use ``angseg`` to fail known obstructed bearing/range sectors in addition to
the vendor flag and optional land mask::

    blocked = [{
        "range_min": 10,
        "range_max": 30,
        "bearingTrue_start": 350,
        "bearingTrue_stop": 20,
    }]
    current.qc_qartod_valid_location(use_mask=False, angseg=blocked)

Segments can wrap through north, as in this example.

Operator flags
~~~~~~~~~~~~~~

Flag a whole file::

    current.qc_operator_test(suspect_all=True)

Or a range/bearing segment::

    current.qc_operator_test(flag_segment=[{
        "range_min": 20,
        "range_max": 40,
        "bearingTrue_start": 90,
        "bearingTrue_stop": 120,
        "flag_value": 4,
    }])

Primary flag
~~~~~~~~~~~~

Run the primary flag last. By default it includes all recorded tests and
ignores not-evaluated results when finding the worst flag. Restrict it to a
named subset when required by a network policy::

    current.qc_qartod_primary_flag(
        include=["qc_qartod_syntax", "qc_qartod_maximum_velocity"]
    )

Re-running tests
~~~~~~~~~~~~~~~~

A test will not overwrite its existing column by default. Pass
``recalculate=True`` to remove the old result and metadata before recomputing::

    current.qc_qartod_maximum_velocity(
        max_speed=200,
        high_speed=120,
        recalculate=True,
    )

Total-vector tests
------------------

::

    total.initialize_qc()
    total.qc_qartod_gdop(max_GDOP=2)                       # Q302
    total.qc_qartod_maximum_velocity(250, 150)             # Q303
    total.qc_qartod_valid_location(use_mask=True)          # Q305
    total.qc_qartod_u_uncertainty(max_Uerr=0.6)            # Q306
    total.qc_qartod_v_uncertainty(max_Verr=0.6)            # Q307
    total.qc_site_contribution(max_ror=4, max_share=0.75)  # SCON
    total.qc_qartod_primary_flag()

Site-contribution QC requires the ``NSIT``, ``ROR``, and ``MSFR`` diagnostics
produced by ``combineRadials(..., engine="vectorized")``. ``qc_threshold`` can
create a project-specific maximum-threshold flag from any numeric column.

Keep flags with the data
------------------------

Filtering only at the last step retains provenance::

    passing = current.data[current.data["PRIM"] == 1]
    usable = current.data[current.data["PRIM"].isin([1, 3])]

Write the complete flagged product before producing a filtered derivative.
Both ``to_ruv`` and ``to_netcdf`` retain supported QC columns.

Production checklist
--------------------

* Store thresholds and land-mask resolution with the processing configuration.
* Keep input filenames and package version with the output provenance.
* Run individual tests before ``PRIM``.
* Use ``recalculate=True`` only when replacing an intentional prior result.
* Validate threshold units against the native file.
* Distinguish ``not evaluated`` from ``pass`` in downstream analysis.
