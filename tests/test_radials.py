import unittest
from pathlib import Path

import numpy as np
import pytest
import sys
import xarray as xr

from hfradarpy.radials import Radial
from hfradarpy.radials import concat as concatenate_radials

data_path = (Path(__file__).parent.with_name("examples") / "data").resolve()
output_path = (Path(__file__).parent.with_name("examples") / "output").resolve()


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
