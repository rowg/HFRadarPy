import unittest
from pathlib import Path

import numpy as np
import numpy.testing as npt
import pandas as pd
import pytest
import sys
import xarray as xr

from hfradarpy.radials import Radial
from hfradarpy.radials import concat as concatenate_radials

data_path = (Path(__file__).parent.with_name("examples") / "data").resolve()
output_path = (Path(__file__).parent.with_name("examples") / "output").resolve()

ROUND_TRIP_FILES = [
    "SEAB/RDLi_SEAB_2019_01_01_0000.ruv",
    "WERA/RDL_csw_2019_10_24_162300.ruv",
    "WERA/RDL_UMiami_STF_2019_06_01_0000.hfrweralluv1.0",
]


class TestRoundTrip:
    """Values written to .ruv and netCDF must match the radial they came from."""

    def test_ruv_round_trip(self, tmp_path):
        rad = Radial(data_path / "radials" / "ruv" / "SEAB" / "RDLi_SEAB_2019_01_01_0000.ruv")
        out_file = tmp_path / "RDLi_SEAB_2019_01_01_0000.ruv"
        rad.to_ruv(out_file, validate=False)

        copy = Radial(out_file)
        pd.testing.assert_frame_equal(copy.data, rad.data)
        for key in ["Site", "TimeStamp", "Origin", "PatternType", "AngularResolution", "RangeResolutionKMeters"]:
            assert copy.metadata[key] == rad.metadata[key]

    def test_ruv_round_trip_keeps_qc_columns(self, tmp_path):
        rad = Radial(data_path / "radials" / "ruv" / "SEAB" / "RDLi_SEAB_2019_01_01_0100.ruv", replace_invalid=False)
        rad.initialize_qc()
        rad.qc_qartod_syntax()
        rad.qc_qartod_maximum_velocity(max_speed=30, high_speed=20)
        rad.qc_qartod_spatial_median()
        rad.qc_qartod_primary_flag()
        out_file = tmp_path / "RDLi_SEAB_2019_01_01_0100.ruv"
        rad.to_ruv(out_file, validate=False)

        copy = Radial(out_file, replace_invalid=False)
        for col in ["Q201", "Q202", "Q205", "PRIM"]:
            npt.assert_array_equal(copy.data[col].to_numpy(), rad.data[col].to_numpy())
        assert set(rad.data["Q202"]) == {1, 3, 4}

    @pytest.mark.parametrize("name", ROUND_TRIP_FILES)
    def test_tabular_netcdf_values_match_source(self, tmp_path, name):
        rad = Radial(data_path / "radials" / "ruv" / name)
        nc_file = tmp_path / "radial.nc"
        rad.to_netcdf(str(nc_file), model="tabular")

        with xr.open_dataset(nc_file) as ds:
            assert ds.sizes["time"] == len(rad.data)
            for var, col in [("velocity", "VELO"), ("u", "VELU"), ("v", "VELV"), ("lon", "LOND"),
                             ("lat", "LATD"), ("bearing", "BEAR"), ("range", "RNGE")]:
                npt.assert_allclose(ds[var].values, rad.data[col].to_numpy(dtype=float), rtol=0, atol=1e-6,
                                    err_msg=var)
            assert ds["velocity"].attrs["units"] == "cm s-1"

    @pytest.mark.parametrize("name", ROUND_TRIP_FILES)
    def test_gridded_netcdf_values_match_source(self, tmp_path, name):
        """Each radial lands in the (range, bearing) cell containing its RNGE and BEAR; every other cell is empty.
        (In the hfrweralluv file RNGE is the true distance of the point, so it falls within a cell, not on its centre.)"""
        rad = Radial(data_path / "radials" / "ruv" / name)
        nc_file = tmp_path / "radial.nc"
        rad.to_netcdf(str(nc_file), model="gridded")

        with xr.open_dataset(nc_file) as ds:
            assert int(np.isfinite(ds["velocity"].values).sum()) == len(rad.data)
            rng = xr.DataArray(rad.data["RNGE"].to_numpy(), dims="row")
            bear = xr.DataArray(rad.data["BEAR"].to_numpy(), dims="row")
            cells = ds.isel(time=0).sel(range=rng, bearing=bear, method="nearest")
            for coord, col in [("range", "RNGE"), ("bearing", "BEAR")]:
                half_cell = np.diff(ds[coord].values.astype(float)).min() / 2
                assert (np.abs(cells[coord].values - rad.data[col].to_numpy()) <= half_cell + 1e-4).all(), coord
            for var, col in [("velocity", "VELO"), ("u", "VELU"), ("v", "VELV")]:
                npt.assert_allclose(cells[var].values, rad.data[col].to_numpy(dtype=float), rtol=0, atol=1e-4,
                                    err_msg=var)

            # Cell lon/lat are cell centres: exact when the radials sit on them, otherwise within a cell (~2 km)
            on_centres = np.allclose(cells["range"].values, rad.data["RNGE"], atol=1e-4)
            for var, col in [("lon", "LOND"), ("lat", "LATD")]:
                npt.assert_allclose(cells[var].values, rad.data[col].to_numpy(dtype=float), rtol=0,
                                    atol=1e-4 if on_centres else 0.025, err_msg=var)


def test_codar_radial_to_tabular_netcdf():
    radial_file = data_path / "radials" / "ruv" / "SEAB" / "RDLi_SEAB_2019_01_01_0000.ruv"
    nc_file = output_path / "radials" / "nc" / "tabular" / "SEAB" / "RDLi_SEAB_2019_01_01_0000.nc"

    # Converts the underlying .data (natively a pandas DataFrame)
    # to an xarray object when `create_netcdf` is called.
    # This automatically 'enhances' the netCDF file
    # with better variable names and attributes.
    rad1 = Radial(radial_file)
    rad1.to_netcdf(str(nc_file), model="tabular", enhance=True)

    # Convert it to an xarray Dataset with no variable
    # or attribute enhancements
    xds2 = rad1.to_xarray("tabular", enhance=False)

    # Convert it to xarray Dataset with increased usability
    # by changing variables names, adding attributes,
    # and decoding the CF standards like scale_factor
    xds3 = rad1.to_xarray("tabular", enhance=True)

    with xr.open_dataset(nc_file) as xds1:
        # The two enhanced files should be identical
        assert xds1.identical(xds3)

        # Enhanced and non-enhanced files should not
        # be equal
        assert not xds1.identical(xds2)


def test_codar_radial_to_gridded_netcdf():
    radial_file = data_path / "radials" / "ruv" / "SEAB" / "RDLi_SEAB_2019_01_01_0000.ruv"
    nc_file = output_path / "radials" / "nc" / "gridded" / "SEAB" / "RDLi_SEAB_2019_01_01_0000.nc"

    # Converts the underlying .data (natively a pandas DataFrame)
    # to an xarray object when `create_netcdf` is called.
    # This automatically 'enhances' the netCDF file
    # with better variable names and attributes.
    rad1 = Radial(radial_file)
    rad1.to_netcdf(str(nc_file), model="gridded")

    # Convert it to an xarray Dataset with no variable
    # or attribte enhancements
    xds2 = rad1.to_xarray("gridded", enhance=False)

    # Convert it to xarray Dataset with increased usability
    # by changing variables names, adding attributes,
    # and decoding the CF standards like scale_factor
    xds3 = rad1.to_xarray("gridded", enhance=True)

    with xr.open_dataset(nc_file) as xds1:
        # The two enhanced files should be identical
        assert xds1.identical(xds3)

        # Enhanced and non-enhanced files should not
        # be equal
        assert not xds1.identical(xds2)


def test_codar_mask():
    radial_file = data_path / "radials" / "ruv" / "SEAB" / "RDLi_SEAB_2019_01_01_0000.ruv"
    rad1 = Radial(radial_file, replace_invalid=False)
    # Total points before masking
    assert len(rad1.data) == 745
    rad1.mask_over_land(subset=True, res='low')
    # Make sure we subset the land points
    assert len(rad1.data) == 592


def test_codar_qc():
    radial_file = data_path / "radials" / "ruv" / "SEAB" / "RDLi_SEAB_2019_01_01_0100.ruv"
    radial_file_previous = data_path / "radials" / "ruv" / "SEAB" / "RDLi_SEAB_2019_01_01_0000.ruv"
    rad1 = Radial(radial_file, replace_invalid=False)
    rad1.initialize_qc()
    # assert len(rad1.data) == 733
    # rad1.mask_over_land(subset=True)
    rad1.qc_qartod_syntax() #QC201
    rad1.qc_qartod_radial_count()
    rad1.qc_qartod_valid_location()
    rad1.qc_qartod_maximum_velocity()
    rad1.qc_qartod_spatial_median()
    rad1.qc_qartod_temporal_gradient(radial_file_previous)
    rad1.qc_qartod_avg_radial_bearing(reference_bearing=180)
    rad1.qc_qartod_primary_flag()
    rad1.qc_qartod_stuck_value()
    # assert len(rad1.data) == 587
    assert "Q201" in rad1.data
    assert "Q202" in rad1.data
    assert "Q203" in rad1.data
    assert "Q204" in rad1.data
    assert "Q205" in rad1.data  # temporal gradient test
    assert "Q206" in rad1.data
    assert "Q207" in rad1.data
    assert "Q209" in rad1.data
    assert "PRIM" in rad1.data


def test_wera_radial_to_tabular_netcdf():
    radial_file = data_path / "radials" / "ruv" / "WERA" / "RDL_csw_2019_10_24_162300.ruv"
    nc_file = output_path / "radials" / "nc" / "tabular" / "WERA" / "RDL_csw_2019_10_24_162300.nc"

    # Converts the underlying .data (natively a pandas DataFrame)
    # to an xarray object when `create_netcdf` is called.
    # This automatically 'enhances' the netCDF file
    # with better variable names and attributes.
    rad1 = Radial(radial_file)
    rad1.to_netcdf(str(nc_file), model="tabular")

    # Convert it to an xarray Dataset with no variable
    # or attribute enhancements
    xds2 = rad1.to_xarray("tabular", enhance=False)

    # Convert it to xarray Dataset with increased usability
    # by changing variables names, adding attributes,
    # and decoding the CF standards like scale_factor
    xds3 = rad1.to_xarray("tabular", enhance=True)

    with xr.open_dataset(nc_file) as xds1:
        # The two enhanced files should be identical
        assert xds1.identical(xds3)

        # Enhanced and non-enhanced files should not
        # be equal
        assert not xds1.identical(xds2)


def test_wera_radial_to_gridded_netcdf():
    radial_file = data_path / "radials" / "ruv" / "WERA" / "RDL_csw_2019_10_24_162300.ruv"
    nc_file = output_path / "radials" / "nc" / "gridded" / "WERA" / "RDL_csw_2019_10_24_162300.nc"

    # Converts the underlying .data (natively a pandas DataFrame)
    # to an xarray object when `create_netcdf` is called.
    # This automatically 'enhances' the netCDF file
    # with better variable names and attributes.
    rad1 = Radial(radial_file)
    rad1.to_netcdf(str(nc_file), model="gridded")

    # Convert it to an xarray Dataset with no variable
    # or attribute enhancements
    xds2 = rad1.to_xarray("gridded", enhance=False)

    # Convert it to xarray Dataset with increased usability
    # by changing variables names, adding attributes,
    # and decoding the CF standards like scale_factor
    xds3 = rad1.to_xarray("gridded", enhance=True)

    with xr.open_dataset(nc_file) as xds1:
        # The two enhanced files should be identical
        assert xds1.identical(xds3)

        # Enhanced and non-enhanced files should not
        # be equal
        assert not xds1.identical(xds2)


def test_wera_mask():
    radial_file = data_path / "radials" / "ruv" / "WERA" / "RDL_csw_2019_10_24_162300.ruv"
    rad1 = Radial(radial_file, replace_invalid=False)
    # Total points before masking
    assert len(rad1.data) == 6327
    rad1.mask_over_land(subset=True, res='low')
    # Make sure we subset the land points
    assert len(rad1.data) == 5745


def test_wera_qc():
    radial_file = data_path / "radials" / "ruv" / "WERA" / "RDL_csw_2019_10_24_162300.ruv"
    rad1 = Radial(radial_file, replace_invalid=False)
    rad1.initialize_qc()
    rad1.mask_over_land(subset=True,res='low')
    rad1.qc_qartod_radial_count()
    rad1.qc_qartod_valid_location(use_mask=True, res='low')
    rad1.qc_qartod_maximum_velocity()
    rad1.qc_qartod_spatial_median()
    rad1.qc_qartod_avg_radial_bearing(reference_bearing=180)
    rad1.qc_qartod_primary_flag()
    assert "Q204" in rad1.data
    assert "Q203" in rad1.data
    assert "Q202" in rad1.data
    assert "Q205" in rad1.data
    assert "Q207" in rad1.data
    assert "PRIM" in rad1.data


def test_wera_raw_to_quality_gridded_nc():
    radial_file = data_path / "radials" / "ruv" / "WERA" / "RDL_csw_2019_10_24_162300.ruv"
    nc_file = output_path / "radials" / "qc" / "nc" / "gridded" / "WERA" / "RDL_csw_2019_10_24_162300.nc"
    rad1 = Radial(radial_file, replace_invalid=False)
    rad1.initialize_qc()
    # rad1.mask_over_land(subset=True)
    rad1.qc_qartod_radial_count()
    rad1.qc_qartod_valid_location()
    rad1.qc_qartod_maximum_velocity()
    rad1.qc_qartod_spatial_median()
    rad1.to_netcdf(str(nc_file), model="gridded")

    xds2 = rad1.to_xarray("gridded", enhance=True)

    with xr.open_dataset(nc_file) as xds1:
        assert len(xds1.QCTest) == 4  # all QC tests should have run
        # The two enhanced files should be identical
        assert xds1.identical(xds2)


def test_wera_raw_to_quality_tabular_nc():
    radial_file = data_path / "radials" / "ruv" / "WERA" / "RDL_csw_2019_10_24_162300.ruv"
    nc_file = output_path / "radials" / "qc" / "nc" / "tabular" / "WERA" / "RDL_csw_2019_10_24_162300.nc"
    rad1 = Radial(radial_file, replace_invalid=False)
    rad1.initialize_qc()
    # rad1.mask_over_land()
    rad1.qc_qartod_radial_count()
    rad1.qc_qartod_valid_location()
    rad1.qc_qartod_maximum_velocity()
    rad1.qc_qartod_spatial_median()
    rad1.to_netcdf(str(nc_file), model="tabular")

    xds2 = rad1.to_xarray("tabular", enhance=True)

    with xr.open_dataset(nc_file) as xds1:
        assert len(xds1.QCTest) == 4  # all QC tests should have run
        # The two enhanced files should be identical
        assert xds1.identical(xds2)


def test_miami_radial_gridded_nc():
    radial_file = data_path / "radials" / "ruv" / "WERA" / "RDL_UMiami_STF_2019_06_01_0000.hfrweralluv1.0"
    nc_file = output_path / "radials" / "nc" / "gridded" / "WERA" / "RDL_UMiami_STF_2019_06_01_0000.nc"

    # Converts the underlying .data (natively a pandas DataFrame)
    # to an xarray object when `create_netcdf` is called.
    # This automatically 'enhances' the netCDF file
    # with better variable names and attributes.
    rad1 = Radial(radial_file)
    rad1.to_netcdf(str(nc_file), model="gridded")

    # Convert it to an xarray Dataset with no variable
    # or attribute enhancements
    xds2 = rad1.to_xarray("gridded", enhance=False)

    # Convert it to xarray Dataset with increased usability
    # by changing variables names, adding attributes,
    # and decoding the CF standards like scale_factor
    xds3 = rad1.to_xarray("gridded", enhance=True)

    with xr.open_dataset(nc_file) as xds1:
        # The two enhanced files should be identical
        assert xds1.identical(xds3)

        # Enhanced and non-enhanced files should not
        # be equal
        assert not xds1.identical(xds2)


def test_miami_radial_tabular_nc():
    radial_file = data_path / "radials" / "ruv" / "WERA" / "RDL_UMiami_STF_2019_06_01_0000.hfrweralluv1.0"
    nc_file = output_path / "radials" / "nc" / "tabular" / "WERA" / "RDL_UMiami_STF_2019_06_01_0000.nc"

    # Converts the underlying .data (natively a pandas DataFrame)
    # to an xarray object when `create_netcdf` is called.
    # This automatically 'enhances' the netCDF file
    # with better variable names and attributes.
    rad1 = Radial(radial_file)
    rad1.to_netcdf(str(nc_file), model="tabular")
    # rad1.export(str(nc_file), file_type='netcdf', model='tabular')

    # Convert it to an xarray Dataset with no variable
    # or attribute enhancements
    xds2 = rad1.to_xarray("tabular", enhance=False)
    # xds2 = rad1.to_xarray_tabular(enhance=False)

    # Convert it to xarray Dataset with increased usability
    # by changing variables names, adding attributes,
    # and decoding the CF standards like scale_factor
    xds3 = rad1.to_xarray("tabular", enhance=True)
    # xds3 = rad1.to_xarray_tabular(enhance=True)

    with xr.open_dataset(nc_file) as xds1:
        # The two enhanced files should be identical
        assert xds1.identical(xds3)

        # Enhanced and non-enhanced files should not
        # be equal
        assert not xds1.identical(xds2)


class TestCombineRadials(unittest.TestCase):
    def setUp(self):
        self.file_paths = list((data_path / "radials" / "ruv" / "SEAB").glob("*.ruv"))

        self.radial_files = [str(r) for r in self.file_paths]

        self.radial_objects = [Radial(str(r)) for r in self.radial_files]

        # Select even indexed file_paths and odd indexed radial objects
        # into one array of mixed content types for concating
        self.radial_mixed = self.radial_files[::2] + self.radial_objects[1:][::2]

    def test_concat_radial_objects(self):
        combined = concatenate_radials(self.radial_objects)
        assert combined.time.size == len(self.file_paths)
        # Make sure the dataset was sorted by time
        assert np.array_equal(combined.time.values, np.sort(combined.time.values))

    def test_concat_radial_files(self):
        combined = concatenate_radials(self.radial_files)
        assert combined.time.size == len(self.file_paths)
        # Make sure the dataset was sorted by time
        assert np.array_equal(combined.time.values, np.sort(combined.time.values))

    def test_concat_mixed_radials(self):
        combined = concatenate_radials(self.radial_mixed)
        assert combined.time.size == len(self.file_paths)
        # Make sure the dataset was sorted by time
        assert np.array_equal(combined.time.values, np.sort(combined.time.values))

    def test_concat_mixed_radials_enhance(self):
        # Select even indexed file_paths and odd indexed radial objects
        # into one array of mixed content types for concating
        combined = concatenate_radials(self.radial_mixed, enhance=True)
        assert combined.time.size == len(self.file_paths)
        # Make sure the dataset was sorted by time
        assert np.array_equal(combined.time.values, np.sort(combined.time.values))
