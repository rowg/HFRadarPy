Plotting
========

HFRadarPy provides plotting functions for enhanced xarray/NetCDF products and
a direct plotter for CODAR radial files. Cartopy and Matplotlib are required.

Plot an enhanced radial dataset
-------------------------------

::

    from hfradarpy.plot.plot_nc import plot_radials

    plot_radials(
        "radial.gridded.nc",
        plot_type="velocity",
        output_file="figures/radial.png",
        extent=[-75, -72, 37, 40],
        velocity_min=-100,
        velocity_max=100,
        sub=2,
        title="SEAB radial currents",
    )

The dataset must contain enhanced names ``u``, ``v``, ``velocity``, ``lon``,
and ``lat``. Available ``plot_type`` values are:

``"velocity"``
   Colour vectors by signed radial velocity.

``"motion"``
   Distinguish toward/away motion with two colours.

``"qc_pass_fail"``
   Colour by ``primary_flag_qc``.

``prim_filter=True`` keeps only rows whose enhanced primary flag equals 1.

Plot a total-vector dataset
---------------------------

::

    from hfradarpy.plot.plot_nc import plot_totals

    plotting_ds = total_ds.rename({"VELU": "u", "VELV": "v"})
    plot_totals(
        plotting_ds,
        output_file="figures/total.png",
        extent=[-76, -69, 35, 42],
        sub=2,
        title="Surface currents",
    )

``plot_totals`` expects one-dimensional ``lon`` and ``lat`` coordinates and
variables named ``u`` and ``v``. HFRadarPy's ``Total.to_xarray`` retains
native ``VELU`` and ``VELV`` names, so rename them as shown.

Directly plot a radial file
---------------------------

::

    from hfradarpy.plot.plot_ruv import plot_ruv

    plot_ruv(
        "RDLi_SEAB_2019_01_01_0000.ruv",
        save_path="figures",
        speed_display="color",
        redblue=True,
        scale=50,
        vlims=(-100, 100),
    )

Use ``plotflag="PRIM"`` to highlight failed or suspect values when that
column is present.

Markers and map extent
----------------------

The NetCDF plotters accept marker specifications as
``[longitude, latitude, matplotlib_kwargs]``::

    markers = [
        [-74.6, 38.5, {"marker": "o", "markersize": 7, "color": "red"}],
    ]

Pass explicit ``extent`` bounds for reproducible maps. Natural Earth coast,
border, lake, and state/province layers are added automatically.

Interactive use
---------------

When ``output_file`` is omitted, the shared plotting function returns a
Matplotlib-like object. In scripts and batch jobs, always set an output path
and close any figures you create outside HFRadarPy.

