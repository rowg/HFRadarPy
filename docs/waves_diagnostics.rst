Waves, diagnostics, and antenna patterns
========================================

Wave files
----------

Load a CODAR ``.wls`` product with :class:`hfradarpy.waves.Waves`::

    from hfradarpy.waves import Waves

    waves = Waves("WVLM_SEAB_2019_03_01_0000.wls")
    print(waves.data.head())

The source layout determines the xarray dimensions. An averaged product uses
``time``; a range-resolved product uses ``time`` and ``DIST``.

Flag implausible wave heights without removing rows::

    waves.flag_wave_heights(minimum=0.2, maximum=5, remove=False)
    print(waves.data["mwht_flag"].value_counts())

The flag is 1 inside the inclusive range and 4 outside. To remove bad rows
instead, pass ``remove=True``.

Convert one file or a series::

    from hfradarpy.waves import concat

    ds = waves.to_xarray(enhance=True)
    series = concat(sorted(wave_directory.glob("*.wls")), enhance=True)

Enhanced wave names include ``wave_height``, ``wave_period``,
``wave_bearing``, ``wind_bearing``, and ``distance_from_origin`` where the
corresponding native fields exist.

Write compressed NetCDF::

    waves.to_netcdf("output/waves.nc", enhance=True)

With ``prepend_extension=True``, the filename indicates ``averaged`` or
``ranged`` output.

Diagnostic files
----------------

:class:`hfradarpy.diagnostics.Diags` reads CODAR receiver, radial, and other
diagnostic CTF files such as ``.rdt``, ``.hdt``, and ``.xdt``::

    from hfradarpy.diagnostics import Diags, concat

    diagnostic = Diags("STAT_SEAB_2018_01_01.rdt")
    print(diagnostic.data.index)  # DatetimeIndex

    history = concat(sorted(diagnostic_directory.glob("*.rdt")))

``concat`` returns a pandas DataFrame sorted according to the input file list.
Write a modified CTF product with ``diagnostic.to_file(path)``.

Antenna-pattern files
---------------------

:class:`hfradarpy.patterns.Loop` parses CODAR ``.loop`` files::

    from hfradarpy.patterns import Loop, concat

    loop = Loop("LOOP_VIEW_2025_10_23_2114.loop")
    print(loop.data.head())
    combined = concat(sorted(loop_directory.glob("*.loop")))

The combined result is a pandas DataFrame. ``Loop.to_file`` writes a CTF
round trip and protects the input filename unless ``overwrite=True``.

Batch-processing advice
-----------------------

* Group files by product and site before concatenation.
* Check ``_iscorrupt`` and ``is_valid()`` on each parsed file.
* Keep native columns when a receiving system expects manufacturer names.
* Use enhanced xarray names for new analytical code.

