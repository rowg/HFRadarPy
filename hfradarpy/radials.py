import datetime as dt
from dateutil.relativedelta import relativedelta
import geopandas as gpd
import logging
import numpy as np
import pandas as pd
from pyproj import Geod, CRS
from collections import OrderedDict
import re
import os
from pathlib import Path
from shapely.geometry import Point
import copy
import xarray as xr
import netCDF4
from hfradarpy.common import fileParser
from hfradarpy.common import timestamp_from_lluv_filename as get_time
from hfradarpy.calc import reckon, dms2dd, createLonLatGridFromBB, createLonLatGridFromBBwera, \
    createLonLatGridFromTopLeftPointWera
from hfradarpy.io.nc import make_encoding
from joblib import Parallel, delayed
import multiprocessing
import json
import warnings
try:
    from mpl_toolkits.basemap import Basemap
except Exception as err:
    pass
import matplotlib.pyplot as plt
from matplotlib import colors
try:
    import cartopy
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import cartopy.io.img_tiles as cimgt
except Exception as err:
    pass
from scipy.spatial import ConvexHull
import geopy.distance

logger = logging.getLogger(__name__)


def velocityMedianInDistLimits(cell, radData, distLim, g):
    """
    This function evaluates the median of all radial velocities contained in radData
    lying within a distance of distLim km from the origin grid cell.
    The native CRS of the Radial is used for distance calculations.

    INPUT:
        cell: Series containing longitude and latitude of the origin grid cell
        radData: DataFrame containing radial data
        distLim: range limit in km for selecting velocities for median calculation
        g: Geod object according to the Radial CRS

    OUTPUT:
        median: median of the selected velocities.
    """
    # Convert grid cell Series and radial bins DataFrame to numpy arrays
    cell = cell.to_numpy()
    radLon = radData['LOND'].to_numpy()
    radLat = radData['LATD'].to_numpy()
    # Evaluate distances between the origin grid cell and radial bins
    az12, az21, cellToRadDist = g.inv(len(radLon) * [cell[0]], len(radLat) * [cell[1]], radLon, radLat)

    # Remove the origin from the radial data (i.e. distance = 0)
    radData = radData.drop(radData[cellToRadDist == 0].index)
    # Remove the origin from the distance array (i.e. distance = 0)
    cellToRadDist = cellToRadDist[cellToRadDist != 0]

    # Figure out which radial bins are within the range limit from the origin cell
    distSelectionIndices = np.where(cellToRadDist < distLim * 1000)[0].tolist()

    # Evaluate the median of the selected velocities
    # median = np.median(radData.iloc[distSelectionIndices]['VELO'])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        median = np.nanmedian(radData.iloc[distSelectionIndices]['VELO'])

    return median


###### COMMENT FROM LORENZO: ######
######I WOULD REMOVE THIS FUNCTION FROM radials.py FILE AND PUT IT IN A HIGHER LEVEL SCRIPT - ANYWAY, IT SEEMS TO BE UNUSED IN THE TOOLBOX ######
def qc_radial_file(radial_file, qc_values=None, export=None, save_path=None, clean=False, clean_path=None):
    """
    Main function to parse and qc radial files.

    Setting clean to True will create two separate quality controlled radial files. Must set clean_path.
    The first radial file with containing flag metadata will be saved to save_path. This file contains data along with flags.
    The second radial file with data that failed qc removed will be saved to clean_path. This file does not contain any flags.


    Args:
        radial_file (str or Path):
            Path to radial file or a Radial object
        qc_values (dict, optional):
            Dictionary containing thresholds for each QC test. Defaults to None.
        export (str, optional):
            None or 'radial' or 'netcdf-tabular' or 'netcdf-multidimensional'. Defaults to None.
        save_path (str or Path, optional):
            Path to save quality controlled radial file. Defaults to None.
        clean (bool, optional):
            Remove any row of data where the primary flag equals 4 (Failure flag). Defaults to False.
        clean_path (str or Path, optional):
            Path to save quality controlled radial file with data that fails qc removed. Defaults to None.

    Returns:
        Radial object: A quality controlled radial file.
    """
    qc_values = qc_values or dict(
        qc_qartod_avg_radial_bearing=dict(reference_bearing=151, warning_threshold=15, failure_threshold=30),
        qc_qartod_radial_count=dict(min_count=75.0, low_count=225.0),
        qc_qartod_maximum_velocity=dict(max_speed=300.0, high_speed=100.0),
        qc_qartod_spatial_median=dict(
            radial_smed_range_cell_limit=2.1, radial_smed_angular_limit=10, radial_smed_current_difference=30
        ),
        qc_qartod_temporal_gradient=dict(gradient_temp_fail=32, gradient_temp_warn=25),
        qc_qartod_primary_flag=dict(
            include=[
                "qc_qartod_syntax",
                "qc_qartod_valid_location",
                "qc_qartod_radial_count",
                "qc_qartod_maximum_velocity",
                "qc_qartod_spatial_median",
            ]
        ),
    )

    if not isinstance(radial_file, Radial):
        r = Radial(radial_file)
    else:
        r = radial_file

    if r.is_valid():
        if clean:
            rclean = copy.deepcopy(r)
        t0 = r.time - dt.timedelta(hours=1)
        previous_radial = "{}_{}{}".format(
            "_".join(r.file_name.split("_")[:2]), t0.strftime("%Y_%m_%d_%H%M"), os.path.splitext(r.file_name)[1]
        )
        previous_full_file = os.path.join(os.path.dirname(r.full_file), previous_radial)
        qc_keys = qc_values.keys()

        # run high frequency radar qartod tests on open radial file
        r.initialize_qc()
        r.qc_qartod_syntax()

        if "qc_qartod_maximum_velocity" in qc_keys:
            r.qc_qartod_maximum_velocity(**qc_values["qc_qartod_maximum_velocity"])

        # run valid location test whether or not there is an option specified in qc_values
        if "qc_qartod_valid_location" in qc_keys:
            r.qc_qartod_valid_location(
                **qc_values["qc_qartod_valid_location"])  # allows the use_mask option to be passed to the test
        else:
            r.qc_qartod_valid_location()

        if "qc_qartod_radial_count" in qc_keys:
            r.qc_qartod_radial_count(**qc_values["qc_qartod_radial_count"])

        if "qc_qartod_spatial_median" in qc_keys:
            r.qc_qartod_spatial_median(**qc_values["qc_qartod_spatial_median"])

        if "qc_qartod_temporal_gradient" in qc_keys:
            r.qc_qartod_temporal_gradient(previous_full_file, **qc_values["qc_qartod_temporal_gradient"])

        if "qc_qartod_stuck_value" in qc_keys:
            r.qc_qartod_stuck_value(**qc_values["qc_qartod_stuck_value"])
            # r.qc_qartod_stuck_value_v2(**qc_values['qc_qartod_stuck_value'])

        if "qc_qartod_avg_radial_bearing" in qc_keys:
            r.qc_qartod_avg_radial_bearing(**qc_values["qc_qartod_avg_radial_bearing"])

        # --------------------------------------------------------------------------
        # Tests that have not been included in the QARTOD manual
        if "qc_qartod_stuck_value_version_2" in qc_keys:
            r.qc_qartod_stuck_value_version_2(**qc_values["qc_qartod_stuck_value"])
            # r.qc_qartod_stuck_value_v2(**qc_values['qc_qartod_stuck_value'])
        # --------------------------------------------------------------------------

        # Primary flag test is performed last
        if "qc_qartod_primary_flag" in qc_keys:
            r.qc_qartod_primary_flag(**qc_values["qc_qartod_primary_flag"])

        if clean:
            d = rclean.data
            dqc = r.data
            if "PRIM" in r.data:
                rt = d[dqc["PRIM"] != 4]
                rclean.data = rt

                for key in rclean._tables.keys():
                    table = rclean._tables[key]
                    if "LLUV" in table["TableType"]:
                        rclean._tables[1]["TableRows"] = rt.shape[0]
                # else:
                #   warning that it didn't update number of table rows
            # else:
            # warning of failure to update file, the original will be exported

        # Export radial file to either a radial or netcdf
        if export:
            if export == "radial":
                r.to_ruv(os.path.join(save_path, r.file_name))
            elif export == "netcdf-tabular":
                r.to_netcdf(os.path.join(save_path, r.file_name), "tabular")
            elif export == "netcdf-gridded":
                r.to_netcdf(os.path.join(save_path, r.file_name), "gridded")

            if clean:
                if export == "radial":
                    rclean.to_ruv(os.path.join(clean_path, rclean.file_name))
                elif export == "netcdf-tabular":
                    rclean.to_netcdf(os.path.join(clean_path, rclean.file_name), "tabular", prepend_extension=True)
                elif export == "netcdf-gridded":
                    rclean.to_netcdf(os.path.join(clean_path, rclean.file_name), "gridded", prepend_extension=True)
        else:
            return r


def concat(rlist, range_minmax=None, bearing=None,
           method="gridded", enhance=False, parallel=False):
    """
    This function takes a list of Radial objects or radial file paths and
    combines them along the time dimension using xarrays built-in concatenation
    routines.

    Args:
        rlist (list):
            list of radial files or Radial objects that you want to concatenate
        method (str, optional):
            'gridded' or 'tabular'. Defaults to 'gridded'
        enhance (bool, optional):
            Change manufacturer variables to cf compliant variable names. Add CF attributes and metadata. Defaults to False.
        parallel (bool, optional): Enables parallel processing. Defaults to False.

    Returns:
        xarray dataset: radials concatenated into an xarray dataset
    """

    def load_radials(radial, method, enhance=False):
        if not isinstance(radial, Radial):
            radial = Radial(radial)
            ds = radial.to_xarray(method, enhance=enhance, range_minmax=range_minmax, bearing=bearing)
        return (radial.file_name, ds)

    if parallel:
        num_cores = multiprocessing.cpu_count()
        radials = Parallel(n_jobs=num_cores)(
            delayed(load_radials)(radial=r, method=method, enhance=enhance) for r in rlist
        )
        radial_dict = {radial: ds for (radial, ds) in radials}
    else:
        radial_dict = {}
        for radial in rlist:
            if not isinstance(radial, Radial):
                radial = Radial(radial)
            if radial.data.shape[0] > 0:
                radial_dict[radial.file_name] = radial.to_xarray(method, enhance=enhance, range_minmax=range_minmax,
                                                                 bearing=bearing)

    ds = xr.concat(radial_dict.values(), "time")
    return ds.sortby("time")


class Radial(fileParser):
    """
    Radial Subclass.

    This class should be used when loading a CODAR radial (.ruv) file. This class utilizes the generic LLUV class from
    ~/hfradarpy/common.py in order to load CODAR Radial files
    """

    def __init__(self, fname, replace_invalid=True, mask_over_land=False, empty_radial=False, vflip=False):
        """
        Initialize the radial object

        Args:
            fname (str or Path): Full file path to the radial to be loaded
            replace_invalid (bool, optional): Replace invalid/dummy manufacturer fill values with NaN. Defaults to True.
            empty_radial (bool, optional): Returns an empty Radial object. Defaults to False.
            mask_over_land (bool, optional): Applies masking over land
            vflip (bool, optional): Flip sign of velocities [MINV, MAXV, VELO] so + is away and - is towards the radar. Defaults to False.
        """
        logging.info("Loading radial file: {}".format(fname))
        super().__init__(fname)

        if self._iscorrupt:
            return

        self.data = pd.DataFrame()
        self.velocity_sign = "+: towards"

        for key in self._tables.keys():
            table = self._tables[key]
            if 'LLUV' in table['TableType']:
                self.data = table['data']
            elif 'CRAD' in table['TableType']:
                self.crad_data = table['data']
            elif 'rads' in table['TableType']:
                self.diagnostics_radial = table['data']
                self.diagnostics_radial['datetime'] = self.diagnostics_radial[
                    ['TYRS', 'TMON', 'TDAY', 'THRS', 'TMIN', 'TSEC']].apply(lambda s: dt.datetime(*s), axis=1)
            elif 'rcvr' in table['TableType']:
                self.diagnostics_hardware = table['data']
                self.diagnostics_hardware['datetime'] = self.diagnostics_hardware[
                    ['TYRS', 'TMON', 'TDAY', 'THRS', 'TMIN', 'TSEC']].apply(lambda s: dt.datetime(*s), axis=1)
            elif 'RINF' in table['TableType']:
                self.range_information = table['data']
            elif 'MRGS' in table['TableType']:
                self.merge_information = table['data']

        if 'Site' in self.metadata.keys():
            self.metadata['Site'] = re.sub(r'[\W_]+', '', self.metadata['Site'])

        if not self.data.empty:
            if replace_invalid:
                self.replace_invalid_values()

            if mask_over_land:
                self.mask_over_land()

            if empty_radial:
                self.empty_radial()

        if vflip:
            # CODAR standard: + velocities towards radar, - velocities away from radar
            # Science standard: - velocities towards radar, + velocities away from radar
            # Flip sign so positive velocities are away from the radar
            self.flip_velocities()

    def __repr__(self):
        """
        Represent class's object as a string
        """
        return "<Radial: {}>".format(self.file_name)

    def flip_velocities(self):
        # CODAR standard: + velocities towards radar, - velocities away from radar
        # Science standard: - velocities towards radar, + velocities away from radar
        # Flip sign so positive velocities are away from the radar as per cf conventions
        logging.info(f"Flipping the sign of radial velocities")
        flips = ["MINV", "MAXV", "VELO"]
        for f in flips:
            if f in self.data:
                self.data[f] = -self.data[f]
        self.flipped_velocites = True

        if self.velocity_sign == "+: towards":
            self.velocity_sign = "-: away"
        elif self.velocity_sign == "-: away":
            self.velocity_sign = "+: towards"

    def empty_radial(self):
        """
        Create an empty Radial object.
        """

        self.file_path = ''
        self.file_name = ''
        self.full_file = ''
        self.metadata = ''
        self._iscorrupt = False
        self.time = []

        for key in self._tables.keys():
            table = self._tables[key]
            self._tables[key]['TableRows'] = '0'
            if 'LLUV' in table['TableType']:
                self.data.drop(self.data.index[:], inplace=True)
                self._tables[key]['data'] = self.data
            elif 'CRAD' in table['TableType']:
                self.crad_data.drop(self.crad_data.index[:], inplace=True)
                self._tables[key]['data'] = self.crad_data
            elif 'rads' in table['TableType']:
                self.diagnostics_radial.drop(self.diagnostics_radial.index[:], inplace=True)
                self._tables[key]['data'] = self.diagnostics_radial
            elif 'rcvr' in table['TableType']:
                self.diagnostics_hardware.drop(self.diagnostics_hardware.index[:], inplace=True)
                self._tables[key]['data'] = self.diagnostics_hardware
            elif 'RINF' in table['TableType']:
                self.range_information.drop(self.range_information.index[:], inplace=True)
                self._tables[key]['data'] = self.range_information

    def mask_over_land(self, subset=False, res='high'):
        """
        This function masks the radial vectors lying on land.
        Radial vector coordinates are checked against a reference file containing information
        about which locations are over land or in an unmeasurable area (for example, behind an
        island or point of land).
        The Natural Earth public domain maps are used as reference.
        If "res" option is set to "high", the map with 10 m resolution is used, otherwise the map with 110 m resolution is used.
        The EPSG:4326 CRS is used for distance calculations.
        If "subset" option is set to True, the radial vectors lying on land are removed.

        INPUT:
            subset: option enabling the removal of radial vectors on land (if set to True)
            res: resolution of the www.naturalearthdata.com dataset used to perform the masking; None or 'low' or 'high'. Defaults to 'high'.

        OUTPUT:
            waterIndex: list containing the indices of radial vectors lying on water.
        """
        # Load the reference file (GeoPandas "naturalearth_lowres")
        mask_dir = Path(__file__).with_name(".hfradarpy")
        if (res == 'high'):
            maskfile = os.path.join(mask_dir, 'ne_10m_admin_0_countries.shp')
        else:
            maskfile = os.path.join(mask_dir, 'ne_110m_admin_0_countries.shp')
        land = gpd.read_file(maskfile)

        # Build the GeoDataFrame containing total points
        geodata = gpd.GeoDataFrame(
            self.data[['LOND', 'LATD']],
            crs="EPSG:4326",
            geometry=[
                Point(xy) for xy in zip(self.data.LOND.values, self.data.LATD.values)
            ]
        )
        # Join the GeoDataFrame containing total points with GeoDataFrame containing leasing areas
        geodata = gpd.sjoin(geodata.to_crs(4326), land.to_crs(4326), how="left", predicate="intersects")

        # All data in the continent column that lies over water should be nan.
        waterIndex = geodata['CONTINENT'].isna()

        if subset:
            # Subset the data to water only
            self.data = self.data.loc[waterIndex].reset_index()
        else:
            return waterIndex

    def plot_Basemap(self, lon_min=None, lon_max=None, lat_min=None, lat_max=None, shade=False, show=True):
        """
        This function plots the current radial velocity field (i.e. VELU and VELV components) on a
        Cartesian grid. The grid is defined either from the input values or from the Radial object
        metadata. If no input is passed and no metadata related to the bounding box are present, the
        grid is defined from data content (i.e. LOND and LATD values).
        If 'shade' is False (default), a quiver plot with color and magnitude of the vectors proportional to
        current velocity is produced. If 'shade' is True, a quiver plot with uniform vetor lenghts is produced,
        superimposed to a pseudo-color map representing velocity magnitude.

        INPUT:
            lon_min: minimum longitude value in decimal degrees (if None it is taken from Total metadata)
            lon_max: maximum longitude value in decimal degrees (if None it is taken from Total metadata)
            lat_min: minimum latitude value in decimal degrees (if None it is taken from Total metadata)
            lat_max: maximum latitude value in decimal degrees (if None it is taken from Total metadata)
            shade: boolean for enabling/disabling shade plot (default False)
            show: boolean for enabling/disabling plot visualization (default True)

        OUTPUT:

        """
        # Initialize figure
        fig = plt.figure(figsize=(24, 16), tight_layout={'pad': 0})

        # Get the bounding box limits
        if not lon_min:
            if 'BBminLongitude' in self.metadata:
                lon_min = float(self.metadata['BBminLongitude'].split()[0])
            else:
                lon_min = self.data.LOND.min() - 1

        if not lon_max:
            if 'BBmaxLongitude' in self.metadata:
                lon_max = float(self.metadata['BBmaxLongitude'].split()[0])
            else:
                lon_max = self.data.LOND.max() + 1

        if not lat_min:
            if 'BBminLatitude' in self.metadata:
                lat_min = float(self.metadata['BBminLatitude'].split()[0])
            else:
                lat_min = self.data.LATD.min() - 1

        if not lat_max:
            if 'BBmaxLatitude' in self.metadata:
                lat_max = float(self.metadata['BBmaxLatitude'].split()[0])
            else:
                lat_max = self.data.LATD.max() + 1

                # Evaluate longitude and latitude of the center of the map
        lon_C = (lon_max + lon_min) / 2
        lat_C = (lat_max + lat_min) / 2

        # Set the background map
        m = Basemap(llcrnrlon=lon_min, llcrnrlat=lat_min, urcrnrlon=lon_max, urcrnrlat=lat_max, lon_0=lon_C,
                    lat_0=lat_C,
                    resolution='i', ellps='WGS84', projection='tmerc')
        m.drawcoastlines()
        # m.fillcontinents(color='#cc9955', lake_color='white')
        m.fillcontinents()
        m.drawparallels(np.arange(lat_min, lat_max))
        m.drawmeridians(np.arange(lon_min, lon_max))
        m.drawmapboundary(fill_color='white')

        # m.bluemarble()

        # Get station coordinates and codes
        if not self.is_wera:
            siteLon = float(self.metadata['Origin'].split()[1])
            siteLat = float(self.metadata['Origin'].split()[0])
            siteCode = self.metadata['Site']
        else:
            if 'Longitude(dd)OfTheCenterOfTheReceiveArray' in self.metadata.keys():
                siteLon = float(self.metadata['Longitude(dd)OfTheCenterOfTheReceiveArray'][:-1])
                siteLat = float(self.metadata['Latitude(dd)OfTheCenterOfTheReceiveArray'][:-1])
                if self.metadata['Latitude(dd)OfTheCenterOfTheReceiveArray'][-1] == 'S':
                    siteLat = -siteLat
                if self.metadata['Longitude(dd)OfTheCenterOfTheReceiveArray'][-1] == 'W':
                    siteLon = -siteLon
            elif 'Longitude(deg-min-sec)OfTheCenterOfTheReceiveArray' in self.metadata.keys():
                siteLon = dms2dd(
                    list(map(int, self.metadata['Longitude(deg-min-sec)OfTheCenterOfTheReceiveArray'][:-2].split('-'))))
                siteLat = dms2dd(
                    list(map(int, self.metadata['Latitude(deg-min-sec)OfTheCenterOfTheReceiveArray'][:-2].split('-'))))
                if self.metadata['Latitude(deg-min-sec)OfTheCenterOfTheReceiveArray'][-1] == 'S':
                    siteLat = -siteLat
                if self.metadata['Longitude(deg-min-sec)OfTheCenterOfTheReceiveArray'][-1] == 'W':
                    siteLon = -siteLon
            siteCode = self.metadata['StationName']

        # Compute the native map projection coordinates for the stations
        xS, yS = m(siteLon, siteLat)

        # Plot radial stations
        m.plot(xS, yS, 'rD')
        plt.text(xS, yS, siteCode, fontdict={'fontsize': 22, 'fontweight': 'bold'})

        # Plot velocity field
        if shade:
            self.to_xarray_multidimensional()

            if self.is_wera:
                # Create grid from longitude and latitude
                [longitudes, latitudes] = np.meshgrid(self.xdr['LONGITUDE'].data, self.xdr['LATITUDE'].data)

                # Compute the native map projection coordinates for the pseudo-color cells
                X, Y = m(longitudes, latitudes)
            else:
                # Compute the native map projection coordinates for the pseudo-color cells
                X, Y = m(self.xdr['LONGITUDE'].data, self.xdr['LATITUDE'].data)

            # Create velocity variable in the shape of the grid
            V = abs(self.xdr['VELO'][0, 0, :, :].to_numpy())
            # V = V[:-1,:-1]

            # Make the pseudo-color plot
            warnings.simplefilter("ignore", category=UserWarning)
            c = m.pcolormesh(X, Y, V, shading='nearest', cmap=plt.cm.jet, vmin=0, vmax=1)

            # Compute the native map projection coordinates for the vectors
            x, y = m(self.data.LOND, self.data.LATD)

            # Create the velocity component variables
            if self.is_wera:
                u = self.data.VELU
                v = self.data.VELV
            else:
                u = self.data.VELU / 100  # CODAR velocities are in cm/s
                v = self.data.VELV / 100  # CODAR velocities are in cm/s

            # Make the quiver plot
            m.quiver(x, y, u, v, width=0.001, headwidth=4, headlength=4, headaxislength=4)

            warnings.simplefilter("default", category=UserWarning)

        else:
            # Compute the native map projection coordinates for the vectors
            x, y = m(self.data.LOND, self.data.LATD)

            # Create the velocity variables
            if self.is_wera:
                u = self.data.VELU
                v = self.data.VELV
                vel = abs(self.data.VELO)
            else:
                u = self.data.VELU / 100  # CODAR velocities are in cm/s
                v = self.data.VELV / 100  # CODAR velocities are in cm/s
                vel = abs(self.data.VELO) / 100  # CODAR velocities are in cm/s

            # Make the quiver plot
            m.quiver(x, y, u, v, vel, cmap=plt.cm.jet, width=0.001, headwidth=4, headlength=4, headaxislength=4)

        # Add colorbar
        cbar = plt.colorbar()
        cbar.set_label('m/s', fontsize='x-large')

        # Add title
        plt.title(self.file_name + ' radial velocity field', fontdict={'fontsize': 30, 'fontweight': 'bold'})

        if show:
            plt.show()

        return fig
    
    def plot(self, lon_min=None, lon_max=None, lat_min=None, lat_max=None, shade=False, show=True):
        """
        This function plots the current radial velocity field (i.e. VELU and VELV components) on a 
        Cartesian grid. The grid is defined either from the input values or from the Radial object
        metadata. If no input is passed and no metadata related to the bounding box are present, the
        grid is defined from data content (i.e. LOND and LATD values).
        If 'shade' is False (default), a quiver plot with color and magnitude of the vectors proportional to
        current velocity is produced. If 'shade' is True, a quiver plot with uniform vetor lenghts is produced,
        superimposed to a pseudo-color map representing velocity magnitude.
        
        INPUT:
            lon_min: minimum longitude value in decimal degrees (if None it is taken from Total metadata)
            lon_max: maximum longitude value in decimal degrees (if None it is taken from Total metadata)
            lat_min: minimum latitude value in decimal degrees (if None it is taken from Total metadata)
            lat_max: maximum latitude value in decimal degrees (if None it is taken from Total metadata)
            shade: boolean for enabling/disabling shade plot (default False)
            show: boolean for enabling/disabling plot visualization (default True)
            
        OUTPUT:

        """
        # Get the bounding box limits
        if not lon_min:
            if 'BBminLongitude' in self.metadata:
                lon_min = float(self.metadata['BBminLongitude'].split()[0])
            else:
                lon_min = self.data.LOND.min() - 1
                
        if not lon_max:
            if 'BBmaxLongitude' in self.metadata:
                lon_max = float(self.metadata['BBmaxLongitude'].split()[0])
            else:
                lon_max = self.data.LOND.max() + 1
                
        if not lat_min:
            if 'BBminLatitude' in self.metadata:
                lat_min = float(self.metadata['BBminLatitude'].split()[0])
            else:
                lat_min = self.data.LATD.min() - 1
                
        if not lat_max:
            if 'BBmaxLatitude' in self.metadata:
                lat_max = float(self.metadata['BBmaxLatitude'].split()[0])
            else:
                lat_max = self.data.LATD.max() + 1   
                
        # Get station coordinates and codes
        if not self.is_wera:
            siteLon = float(self.metadata['Origin'].split()[1])
            siteLat = float(self.metadata['Origin'].split()[0])
            siteCode = self.metadata['Site']
        else:
            if 'Longitude(dd)OfTheCenterOfTheReceiveArray' in self.metadata.keys():
                siteLon = float(self.metadata['Longitude(dd)OfTheCenterOfTheReceiveArray'][:-1])
                siteLat = float(self.metadata['Latitude(dd)OfTheCenterOfTheReceiveArray'][:-1])
                if self.metadata['Latitude(dd)OfTheCenterOfTheReceiveArray'][-1] == 'S':
                    siteLat = -siteLat
                if self.metadata['Longitude(dd)OfTheCenterOfTheReceiveArray'][-1] == 'W':
                    siteLon = -siteLon
            elif 'Longitude(deg-min-sec)OfTheCenterOfTheReceiveArray' in self.metadata.keys():
                siteLon = dms2dd(list(map(int,self.metadata['Longitude(deg-min-sec)OfTheCenterOfTheReceiveArray'][:-2].split('-'))))
                siteLat = dms2dd(list(map(int,self.metadata['Latitude(deg-min-sec)OfTheCenterOfTheReceiveArray'][:-2].split('-'))))
                if self.metadata['Latitude(deg-min-sec)OfTheCenterOfTheReceiveArray'][-1] == 'S':
                    siteLat = -siteLat
                if self.metadata['Longitude(deg-min-sec)OfTheCenterOfTheReceiveArray'][-1] == 'W':
                    siteLon = -siteLon
            siteCode = self.metadata['StationName']
                
        # Initialize the figure     
        cmap = 'jet'                                                                                                        # set the colorbar
        norm = colors.Normalize(vmin=0, vmax=1)                                                                             # set colorbar limits
        extent = [lon_min, lon_max,lat_min, lat_max]                                                                        # set the map extent
        fig = plt.figure(num=None, figsize=(24, 24), dpi=100, facecolor='w', edgecolor='k')
        ax = plt.axes(projection=ccrs.Mercator())                                                                           # set the map projection
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--') # add grid lines
        ax.add_feature(cfeature.LAND)                                                                                       # add land
        ax.add_feature(cfeature.OCEAN)                                                                                      # add ocean
        ax.add_feature(cfeature.COASTLINE)                                                                                  # add coastline
        ax.set_extent(extent)                                                                                               # apply the map extent
        
        # Plot radial stations
        plt.plot(siteLon, siteLat, color='r', markeredgecolor='k', marker='o', markersize=10, transform=ccrs.Geodetic())
        ax.text(siteLon, siteLat, siteCode, transform=ccrs.Geodetic(), fontdict={'fontsize': 22, 'fontweight' : 'bold'})        
        
        # Scale the velocity component variables
        if self.is_wera:
            u = self.data.VELU
            v = self.data.VELV
            vel = abs(self.data.VELO)
        else:
            u = self.data.VELU / 100                # CODAR velocities are in cm/s
            v = self.data.VELV / 100                # CODAR velocities are in cm/s
            vel = abs(self.data.VELO) / 100         # CODAR velocities are in cm/s     
        
        # Plot velocity field
        if shade:
            self.to_xarray_multidimensional()
            
            # Create grid from longitude and latitude
            if self.is_wera:
                # Create grid from longitude and latitude
                [X, Y] = np.meshgrid(self.xdr['LONGITUDE'].data, self.xdr['LATITUDE'].data)
            else:
                X, Y = self.xdr['LONGITUDE'].data, self.xdr['LATITUDE'].data
            
            # Create velocity variable in the shape of the grid
            M = abs(self.xdr['VELO'][0,0,:,:].to_numpy())
            # V = V[:-1,:-1]            
            
            # Make the pseudo-color plot
            warnings.simplefilter("ignore", category=UserWarning)
            ax.pcolormesh(X, Y, M, transform=ccrs.PlateCarree(), shading='nearest', cmap=cmap, vmin=0, vmax=1)
            
            # Get the longitude, latitude and velocity components
            x,y,U,V = self.data['LOND'].values, self.data['LATD'].values, u.values, v.values
            
            # Evaluate the velocity magnitude
            m = (U ** 2 + V ** 2) ** 0.5
            
            # Normalize velocity components
            Un, Vn = U/m, V/m
            
            # Make the quiver plot
            Q = ax.quiver(x, y, Un, Vn, transform=cartopy.crs.PlateCarree(), width=0.001, headwidth=4, headlength=4, headaxislength=4, cmap=cmap, norm=norm, scale=20)
            
            warnings.simplefilter("default", category=UserWarning)
            
        else:
            # Get the longitude, latitude and velocity components
            x,y,U,V = self.data['LOND'].values, self.data['LATD'].values, u.values, v.values
            
            # Make the quiver plot
            Q = ax.quiver(x, y, U, V, vel, transform=cartopy.crs.PlateCarree(), cmap=cmap, norm=norm, scale=5)
            # Add the reference arrow
            ax.quiverkey(Q, 0.1, 0.9, 0.5, r'$0.5 m/s$',fontproperties={'size': 12,'weight': 'bold'})
            
        # Add colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap,norm=norm)
        plt.colorbar(sm,ax=ax, orientation='vertical', pad=0.1).ax.set_xlabel('velocity (m/s)', labelpad=10, )

        
        # Add title
        plt.title(self.file_name + ' radial velocity field', fontdict={'fontsize': 30, 'fontweight' : 'bold'}, pad=25)

                
        if show:
            plt.show()
        
        return fig

    def to_xarray(self, model="tabular", enhance=False, range_minmax=None, bearing=None):
        """
        Helper function

        Args:
            model (str, optional):
                Create a 'tabular' (time) or 'gridded' (time, range, bearing) xarray dataset. Defaults to 'tabular'.
            enhance (bool, optional):
                Change manufacturer variables to meaningful variable names. Add attributes and other metadata. Defaults to True
            range_minmax (list or tuple, optional):
                For gridded dataset only, Range minimum (km) and range maximum (km) of the radial.
                Defaults to None which automatically infers the range based off of the radial data.
            bearing (list or tuple, optional):
                For gridded dataset only, Antenna bearing (degrees) and angular resolution (degrees) of the radial grid.
                Defaults to None which automatically infers the bearing based off of the radial data.
        """
        # Make a copy so that the original radial object is not altered by this function
        rcopy = copy.deepcopy(self)

        if model == "tabular":
            ds = rcopy._to_xarray_tabular(enhance)
        elif model == "gridded":
            ds = rcopy._to_xarray_gridded(range_minmax=range_minmax, bearing=bearing, enhance=enhance)
        else:
            raise ValueError("Please enter a valid data model type. Must be a string 'tabular' or 'gridded'")
        return ds

    def _to_xarray_gridded(self, range_minmax=None, bearing=None, enhance=True):
        """
        Convert radial file to an xarray dataset on a radial grid.
        This dataset has dimensions of time,range, and bearing.
        Add list range_min and range_max if you want to make the file only include range cells between the two.

        Adapted from MATLAB code by Mark Otero
        http://cordc.ucsd.edu/projects/mapping/documents/HFRNet_Radial_NetCDF.pdf

        Args:
            range_minmax (list or tuple, optional): Range minimum (km) and range maximum (km) of the radial.
                Defaults to None which automatically infers the range based off of the radial data.
            bearing (list or tuple, optional): Antenna bearing (degrees) and angular resolution (degrees) of the radial grid.
                Defaults to None which automatically infers the bearing based off of the radial data.
            enhance (bool, optional):
                Change manufacturer variables to cf compliant variable names. Add CF attributes and metadata. Defaults to True.

        Returns:
            xarray.Dataset: radial file converted to an xarray dataset with dimensions of time, range, and bearing.
        """
        logging.info("Converting radial matrix to multidimensional dataset")

        # CF Standard: T, Z, Y, X
        coords = ("time", "range", "bearing")

        # Intitialize empty xarray dataset
        ds = xr.Dataset()

        if range_minmax:
            # range_include = True, because it is not a NoneType
            if isinstance(range_minmax, (list, tuple)):
                # range_include must be a tuple or a list that has data inside
                if isinstance(range_minmax, (tuple)):
                    range_min = np.float32(range_minmax[0])
                    range_max = np.float32(range_minmax[1])
            else:
                # range_minmax is not a tuple or a list, so error out and return to the original function
                raise TypeError(
                    "range_minmax must be a list/tuple or None. See function documentation"
                )
        else:
            # range_minmax = False, because range_minmax is a NoneType.
            logging.info("Inferring min/max ranges based off of radial data")
            range_min = self.data.RNGE.min()
            range_max = self.data.RNGE.max()

        if isinstance(range_minmax, (list)):
            range_dim = range_minmax
        else:
            # Get the range resolution of the radar site from the metadata
            range_step = np.double(self.metadata["RangeResolutionKMeters"].split()[0])

            # Create an array from range_min to range_max with range_step between each range
            range_dim = np.arange(range_min, np.round(range_max + range_step, 4), range_step, dtype=np.float32)

        # round floats to 4 decimal places due to inaccuracy of np.arange
        range_dim = [round(x, 4) for x in range_dim]

        # If bearing has user input and
        if bearing:
            # bearing = True, because it is not a NoneType
            if isinstance(bearing, (list, tuple)):
                # bearing must be a tuple or a list contains data
                antenna_bearing = np.float32(bearing[0])
                angular_resolution = np.float32(bearing[1])
            else:
                # bearing is not a tuple or a list, so error out and return to the original function
                raise TypeError(
                    "bearing must be a list/tuple or None. See function documentation"
                )
        else:
            # bearing is a NoneType
            try:
                # Infer it from the radial metadata
                angular_resolution = np.double(self.metadata["AngularResolution"].split()[0])
                antenna_bearing = np.double(self.metadata["AntennaBearing"].split()[0])
            except KeyError:
                logging.info(
                    "AngularResolution and/or AntennaBearing not find in radial file metadata. "
                    "Creating a generic 360 degree grid with 1 degree bin resolution."
                )
                # If AngularResolution and/or AntennaBearing aren't in the radial file,
                # let's create a generic 360 degree grid with 1 degree resolution.
                angular_resolution = 1
                antenna_bearing = 0

        if isinstance(bearing, (list)):
            bearing_dim = bearing
        else:
            # Calculate bearings with the AntennaBearing as the min
            # and AntennaBearing+360 as the max
            # with angular_resolution as the spacing between each bearing.
            #  We use np.mod(bearings, 360) to convert back to 0 to 360
            bearing_dim = sorted(
                np.mod(np.arange(antenna_bearing, antenna_bearing + 360, angular_resolution, dtype=np.short), 360)
            )

        # create radial grid from bearing and range
        [bearings, ranges] = np.meshgrid(bearing_dim, range_dim)

        # calculate lat/lons from origin, bearing, and ranges
        latlon = [float(x) for x in re.findall(r"[-+]?\d*\.\d+|\d+", self.metadata["Origin"])]
        lons = np.full_like(bearings, latlon[1], dtype=float)
        lats = np.full_like(bearings, latlon[0], dtype=float)
        lond, latd = reckon(lons, lats, bearings, ranges)

        # create dictionary containing variables from dataframe in the shape of radial grid
        d = {key: np.tile(np.nan, bearings.shape) for key in self.data.keys()}

        # find grid indices from radial grid (bearing, ranges)
        range_map_idx = np.tile(np.nan, self.data["RNGE"].shape)
        bearing_map_idx = np.tile(np.nan, self.data["BEAR"].shape)

        for i, line in enumerate(self.data["RNGE"]):
            range_map_idx[i] = np.argmin(np.abs(range_dim - self.data.RNGE[i]))
            bearing_map_idx[i] = np.argmin(np.abs(bearing_dim - self.data.BEAR[i]))

        for k, v in d.items():
            v[range_map_idx.astype(int), bearing_map_idx.astype(int)] = self.data[k]
            d[k] = v

        # Add extra dimension for time
        d = {k: np.expand_dims(np.float32(v), axis=0) for (k, v) in d.items()}

        # Add coordinate variables to dataset
        timestamp = dt.datetime(*[int(s) for s in self.metadata["TimeStamp"].split()])
        ds.coords["bearing"] = bearing_dim
        ds.coords["range"] = range_dim
        ds.coords["time"] = pd.date_range(timestamp, periods=1)
        ds.coords["lon"] = (("range", "bearing"), lond.round(4))
        ds.coords["lat"] = (("range", "bearing"), latd.round(4))

        # Add all variables to dataset
        for k, v in d.items():
            ds[k] = (coords, v)

        # Check if calculated longitudes and latitudes align with given longitudes and latitudes
        # plt.plot(ds.lon, ds.lat, 'bo', ds.LOND.squeeze(), ds.LATD.squeeze(), 'rx')

        # Drop extraneous variables
        ds = ds.drop_vars(["LOND", "LATD", "BEAR", "RNGE"])

        # Assign header data to global attributes
        # convert QCTest dictionary to a list for NetCDF conversion
        if 'QCTest' in self.metadata:
            self.metadata["QCTest"] = list(self.metadata["QCTest"].values())
        ds = ds.assign_attrs(self.metadata)

        if enhance is True:
            ds = self.enhance_xarray(ds)
            ds = xr.decode_cf(ds)

        return ds

    def _to_xarray_tabular(self, enhance=False):
        """
        Convert radial file to an xarray dataset on a radial grid.
        This dataset has one dimension, time.

        Args:
            enhance (bool, optional): Rename variables to something meaningful and add cf attributes. Defaults to False.

        Returns:
            xarray.Dataset: Enhanced dataset with renamed variables and useful attributes.
        """
        logging.info("Converting radial matrix to tabular dataset")

        # Clean radial header
        # self.clean_header()

        # get timestamp from radial metadata
        timestamp = dt.datetime(*[int(s) for s in self.metadata["TimeStamp"].split()])

        self.data["time"] = timestamp
        self.data.set_index("time", inplace=True)

        # Intitialize xarray dataset
        ds = self.data.to_xarray()

        # Check if calculated longitudes and latitudes align with given longitudes and latitudes
        # plt.plot(ds.lon, ds.lat, 'bo', ds.LOND.squeeze(), ds.LATD.squeeze(), 'rx')

        # Assign header data to global attributes
        # convert QCTest dictionary to a list for NetCDF conversion
        if 'QCTest' in self.metadata:
            self.metadata["QCTest"] = list(self.metadata["QCTest"].values())
        ds = ds.assign_attrs(self.metadata)

        if enhance is True:
            ds = self.enhance_xarray(ds)
            ds = xr.decode_cf(ds)

        return ds

    def to_xarray_multidimensional(self, range_min=None, range_max=None, enhance=False):
        """
        This function creates a dictionary of xarray DataArrays containing the variables
        of the Radial object multidimensionally expanded along the coordinate axes (T,Z,Y,X).
        The spatial coordinate axes are chosen based on the file type: (RNGE,BEAR) for
        Codar radials and (LOND,LATD) for WERA radials. The coordinate limits and steps
        are taken from radial metadata. Only RNGE limits can be specified by the user.
        Some refinements are performed on Codar data in order to comply CF convention
        for positive velocities and to have velocities in m/s.
        The generated dictionary is attached to the Radial object, named as xdr.

        INPUT:
            range_min: minimum range value in km (if None it is taken from Radial metadata)
            range_max: maximum range value in km (if None it is taken from Radial metadata)

        OUTPUT:
        """
        # Initialize empty dictionary
        xdr = OrderedDict()

        # process Codar radial
        if not self.is_wera:
            # Check range limits
            if range_min is None:
                if 'RangeMin' in self.metadata:
                    range_min = float(self.metadata['RangeMin'].split()[0])
                else:
                    range_min = self.data.RNGE.min()
            if range_max is None:
                if 'RangeMax' in self.metadata:
                    range_max = float(self.metadata['RangeMax'].split()[0])
                else:
                    range_max = self.data.RNGE.max()
            # Get range step
            if 'RangeResolutionKMeters' in self.metadata:
                range_step = float(self.metadata['RangeResolutionKMeters'].split()[0])
            elif 'RangeResolutionMeters' in self.metadata:
                range_step = float(self.metadata['RangeResolutionMeters'].split()[0]) * 0.001
            else:
                range_step = float(1)
            # build range array
            range_dim = np.arange(
                range_min,
                range_max + range_step,
                range_step
            )

            # Check that the constraint on maximum number of range cells is satisfied
            range_dim = range_dim[range_dim <= range_max]

            # Get bearing step
            if 'AngularResolution' in self.metadata:
                bearing_step = float(self.metadata['AngularResolution'].split()[0])
            else:
                bearing_step = float(1)
            # build bearing array
            # bearing_dim = np.arange(1, 361, 1).astype(np.float)  # Complete 360 degree bearing coordinate allows for better aggregation
            bearing_dim_1 = np.sort(np.arange(np.min(self.data['BEAR']), -bearing_step, -bearing_step))
            bearing_dim_2 = np.sort(
                np.arange(np.min(self.data['BEAR']), np.max(self.data['BEAR']) + bearing_step, bearing_step))
            bearing_dim_3 = np.sort(np.arange(np.max(self.data['BEAR']), 360, bearing_step))
            bearing_dim = np.unique(np.concatenate((bearing_dim_1, bearing_dim_2, bearing_dim_3), axis=None))
            bearing_dim = bearing_dim[(bearing_dim >= 0) & (bearing_dim <= 360)]

            # Create radial grid from bearing and range
            [bearing, ranges] = np.meshgrid(bearing_dim, range_dim)

            # Get the ellipsoid of the radial coordinate reference system
            if 'GreatCircle' in self.metadata:
                radEllps = self.metadata['GreatCircle'].split()[0].replace('"', '')
            else:
                radEllps = 'WGS84'
            # Create Geod object with the retrieved ellipsoid
            g = Geod(ellps=radEllps)

            # Calculate lat/lons from origin, bearing, and ranges
            latlon = [float(x) for x in re.findall(r"[-+]?\d*\.\d+|\d+", self.metadata['Origin'])]
            bb = bearing.flatten()
            rr = ranges.flatten() * 1000  # distances to be expressed in meters for Geod
            lond, latd, backaz = g.fwd(len(bb) * [latlon[1]], len(bb) * [latlon[0]], bb, rr)
            lond = np.array(lond)
            latd = np.array(latd)
            lond = lond.reshape(bearing.shape)
            latd = latd.reshape(bearing.shape)

            # Find grid indices from radial grid (bearing, ranges)
            range_map_idx = np.tile(np.nan, self.data['RNGE'].shape)
            bearing_map_idx = np.tile(np.nan, self.data['BEAR'].shape)

            for i, line in enumerate(self.data['RNGE']):
                range_map_idx[i] = np.argmin(np.abs(range_dim - self.data.RNGE[i]))
                bearing_map_idx[i] = np.argmin(np.abs(bearing_dim - self.data.BEAR[i]))

            # Set X and Y coordinate mappings
            X_map_idx = bearing_map_idx  # BEAR is X axis
            Y_map_idx = range_map_idx  # RNGE is Y axis

            # Create dictionary containing variables from dataframe in the shape of radial grid
            d = {key: np.tile(np.nan, bearing.shape) for key in self.data.keys()}

        # process WERA radials
        else:
            # Get longitude limits and step
            if 'TopLeftLongitude' in self.metadata:
                topLeftLon = float(self.metadata['TopLeftLongitude'].split()[0])
            else:
                topLeftLon = float(0)
            if 'nx' in self.metadata:
                cellsLon = int(self.metadata['nx'].split()[0])
            else:
                cellsLon = 100

            # Get latitude limits and step
            if 'TopLeftLatitude' in self.metadata:
                topLeftLat = float(self.metadata['TopLeftLatitude'].split()[0])
            else:
                topLeftLat = float(90)
            if 'ny' in self.metadata:
                cellsLat = int(self.metadata['ny'].split()[0])
            else:
                cellsLat = 100

            # Get cell size in km
            if 'DGT' in self.metadata:
                cellSize = float(self.metadata['DGT'].split()[0])
            else:
                cellSize = float(2)

            # Generate grid coordinates
            gridGS = createLonLatGridFromTopLeftPointWera(topLeftLon, topLeftLat, cellSize, cellsLon, cellsLat)
            # extract longitudes and latitude from grid GeoSeries and insert them into numpy arrays
            lon_dim = np.unique(gridGS.x.to_numpy())
            lat_dim = np.unique(gridGS.y.to_numpy())
            # manage antimeridian crossing
            lon_dim = np.concatenate((lon_dim[lon_dim >= 0], lon_dim[lon_dim < 0]))

            # # Get the longitude and latitude values of the radial measurements
            # unqLon = np.sort(np.unique(self.data['LOND']))
            # unqLat = np.sort(np.unique(self.data['LATD']))

            # # Insert unqLon and unqLat values to replace the closest in lon_dim and lat_dim
            # replaceIndLon = abs(unqLon[None, :] - lon_dim[:, None]).argmin(axis=0).tolist()
            # replaceIndLat = abs(unqLat[None, :] - lat_dim[:, None]).argmin(axis=0).tolist()
            # lon_dim[replaceIndLon] = unqLon
            # lat_dim[replaceIndLat] = unqLat

            # Create radial grid from longitude and latitude
            [longitudes, latitudes] = np.meshgrid(lon_dim, lat_dim)

            # Find grid indices from lon/lat grid (longitudes, latitudes)
            lat_map_idx = np.tile(np.nan, self.data['LATD'].shape)
            lon_map_idx = np.tile(np.nan, self.data['LOND'].shape)

            for i, line in enumerate(self.data['LATD']):
                lat_map_idx[i] = np.argmin(np.abs(lat_dim - self.data.LATD[i]))
                lon_map_idx[i] = np.argmin(np.abs(lon_dim - self.data.LOND[i]))

            # set X and Y coordinate mappings
            X_map_idx = lon_map_idx  # LONGITUDE is X axis
            Y_map_idx = lat_map_idx  # LATITUDE is Y axis

            # create dictionary containing variables from dataframe in the shape of radial grid
            d = {key: np.tile(np.nan, longitudes.shape) for key in self.data.keys()}

            # Remap all variables
        for k, v in d.items():
            v[Y_map_idx.astype(int), X_map_idx.astype(int)] = self.data[k]
            d[k] = v

        # Add extra dimensions for time (T) and depth (Z) - CF Standard: T, Z, Y, X -> T=axis0, Z=axis1
        d = {k: np.expand_dims(np.float32(v), axis=(0, 1)) for (k, v) in d.items()}

        # Drop LOND, LATD, BEAR and RNGE variables (they are set as coordinates of the DataSet)
        if self.is_wera:
            d.pop('LOND')
            d.pop('LATD')
        else:
            d.pop('BEAR')
            d.pop('RNGE')
            d.pop('LOND')
            d.pop('LATD')

            # Refine Codar data
        if not self.is_wera:
            # Flip sign so positive velocities are away from the radar as per CF conventions (only for Codar radials)
            flips = ['MINV', 'MAXV', 'VELO']
            for f in flips:
                if f in d:
                    d[f] = -d[f]
            # Scale velocities to be in m/s (only for Codar radials)
            toMs = ['VELU', 'VELV', 'VELO', 'ESPC', 'ETMP', 'MINV', 'MAXV']
            for t in toMs:
                if t in d:
                    d[t] = d[t] * 0.01

        # Evaluate timestamp as number of days since 1950-01-01T00:00:00Z
        timeDelta = self.time - dt.datetime.strptime('1950-01-01T00:00:00Z', '%Y-%m-%dT%H:%M:%SZ')
        ncTime = timeDelta.days + timeDelta.seconds / (60 * 60 * 24)

        # Add all variables as xarray
        if not self.is_wera:
            for k, v in d.items():
                xdr[k] = xr.DataArray(v,
                                      dims={'TIME': v.shape[0], 'DEPTH': v.shape[1], 'RNGE': v.shape[2],
                                            'BEAR': v.shape[3]},
                                      coords={'TIME': [ncTime],
                                              'DEPTH': [0],
                                              'RNGE': range_dim,
                                              'BEAR': bearing_dim})
        else:
            for k, v in d.items():
                xdr[k] = xr.DataArray(v,
                                      dims={'TIME': v.shape[0], 'DEPTH': v.shape[1], 'LATITUDE': v.shape[2],
                                            'LONGITUDE': v.shape[3]},
                                      coords={'TIME': [ncTime],
                                              'DEPTH': [0],
                                              'LATITUDE': lat_dim,
                                              'LONGITUDE': lon_dim})

        # Add longitudes and latitudes evaluated from bearing/range grid to the dictionary (BEAR and RNGE are coordinates) for Codar data
        if not self.is_wera:
            xdr['LONGITUDE'] = xr.DataArray(np.float32(lond),
                                            dims={'RNGE': lond.shape[0], 'BEAR': lond.shape[1]},
                                            coords={'RNGE': range_dim,
                                                    'BEAR': bearing_dim})
            xdr['LATITUDE'] = xr.DataArray(np.float32(latd),
                                           dims={'RNGE': latd.shape[0], 'BEAR': latd.shape[1]},
                                           coords={'RNGE': range_dim,
                                                   'BEAR': bearing_dim})

            # Add DataArray for coordinate variables
        xdr['TIME'] = xr.DataArray(ncTime,
                                   dims={'TIME': len(pd.date_range(self.time, periods=1))},
                                   coords={'TIME': [ncTime]})
        xdr['DEPTH'] = xr.DataArray(0,
                                    dims={'DEPTH': 1},
                                    coords={'DEPTH': [0]})
        if not self.is_wera:
            xdr['RNGE'] = xr.DataArray(range_dim,
                                       dims={'RNGE': range_dim},
                                       coords={'RNGE': range_dim})
            xdr['BEAR'] = xr.DataArray(bearing_dim,
                                       dims={'BEAR': bearing_dim},
                                       coords={'BEAR': bearing_dim})
        else:
            xdr['LATITUDE'] = xr.DataArray(lat_dim,
                                           dims={'LATITUDE': lat_dim},
                                           coords={'LATITUDE': lat_dim})
            xdr['LONGITUDE'] = xr.DataArray(lon_dim,
                                            dims={'LONGITUDE': lon_dim},
                                            coords={'LONGITUDE': lon_dim})

            # Attach the dictionary to the Radial object
        self.xdr = xdr

        return

    @staticmethod
    def enhance_xarray(xds):
        """
        Change manufacturer variables to meaningful variable names.
        Add attributes and other metadata.

        Args:
            xds (xarray.Dataset): radial xarray Dataset to enhance

        Returns:
            xds (xarray.Dataset): enhanced xarray dataset
        """
        rename = dict(
            VELU="u",
            VELV="v",
            VFLG="vector_flag",
            ESPC="spatial_quality",
            ETMP="temporal_quality",
            MAXV="velocity_max",
            MINV="velocity_min",
            ERSC="spatial_count",
            ERTC="temporal_count",
            XDST="dist_east_from_origin",
            YDST="dist_north_from_origin",
            VELO="velocity",
            HEAD="heading",
            SPRC="range_cell",
            EACC="accuracy",  # WERA specific
            LOND="lon",
            LATD="lat",
            BEAR="bearing",
            RNGE="range",
        )

        rename_qc = dict()

        # rename variables to something meaningful if they existin
        # in the xarray dataset
        existing_renames = {k: v for k, v in rename.items() if k in xds}
        xds = xds.rename(existing_renames)

        # set time attribute
        xds["time"].attrs["standard_name"] = "time"

        # Set lon attributes
        xds["lon"].attrs["long_name"] = "Longitude"
        xds["lon"].attrs["standard_name"] = "longitude"
        xds["lon"].attrs["short_name"] = "lon"
        xds["lon"].attrs["units"] = "degrees_east"
        xds["lon"].attrs["axis"] = "X"
        xds["lon"].attrs["valid_min"] = np.float32(-180.0)
        xds["lon"].attrs["valid_max"] = np.float32(180.0)
        # xds['lon'] = xds['lon'].astype(np.float64)

        # Set lat attributes
        xds["lat"].attrs["long_name"] = "Latitude"
        xds["lat"].attrs["standard_name"] = "latitude"
        xds["lat"].attrs["short_name"] = "lat"
        xds["lat"].attrs["units"] = "degrees_north"
        xds["lat"].attrs["axis"] = "Y"
        xds["lat"].attrs["valid_min"] = np.float32(-90.0)
        xds["lat"].attrs["valid_max"] = np.float32(90.0)
        # xds['lat'] = xds['lat'].astype(np.float64)

        # Set u attributes
        xds["u"].attrs["long_name"] = "Eastward Surface Current (cm/s)"
        xds["u"].attrs["standard_name"] = "surface_eastward_sea_water_velocity"
        xds["u"].attrs["short_name"] = "u"
        xds["u"].attrs["units"] = "cm s-1"
        xds["u"].attrs["valid_min"] = np.float32(-300)
        xds["u"].attrs["valid_max"] = np.float32(300)
        xds["u"].attrs["coordinates"] = "lon lat"
        xds["u"].attrs["grid_mapping"] = "crs"
        # xds['u'] = xds['u'].astype(np.float64)

        # Set v attributes
        xds["v"].attrs["long_name"] = "Northward Surface Current (cm/s)"
        xds["v"].attrs["standard_name"] = "surface_northward_sea_water_velocity"
        xds["v"].attrs["short_name"] = "v"
        xds["v"].attrs["units"] = "cm s-1"
        xds["v"].attrs["valid_min"] = np.float32(-300)
        xds["v"].attrs["valid_max"] = np.float32(300)
        xds["v"].attrs["coordinates"] = "lon lat"
        xds["v"].attrs["grid_mapping"] = "crs"
        # xds['v'] = xds['v'].astype(np.float64)

        # Set bearing attributes
        xds["bearing"].attrs["long_name"] = "Bearing from origin (away from instrument)"
        xds["bearing"].attrs["short_name"] = "bearing"
        xds["bearing"].attrs["units"] = "degrees"
        xds["bearing"].attrs["valid_min"] = np.float32(0)
        xds["bearing"].attrs["valid_max"] = np.float32(360)
        xds["bearing"].attrs["grid_mapping"] = "crs"
        xds["bearing"].attrs["axis"] = "Y"
        # xds['bearing'] = xds['bearing'].astype(np.float64)

        # Set range attributes
        xds["range"].attrs["long_name"] = "Range from origin (away from instrument)"
        xds["range"].attrs["short_name"] = "range"
        xds["range"].attrs["units"] = "km"
        xds["range"].attrs["valid_min"] = np.float32(0)
        xds["range"].attrs["valid_max"] = np.float32(1000)
        xds["range"].attrs["grid_mapping"] = "crs"
        xds["range"].attrs["axis"] = "X"
        # xds['range'] = xds['range'].astype(np.float64)

        # velocity
        xds["velocity"].attrs["valid_range"] = [-1000, 1000]
        xds["velocity"].attrs["standard_name"] = "radial_sea_water_velocity_away_from_instrument"
        xds["velocity"].attrs["units"] = "cm s-1"
        xds["velocity"].attrs["coordinates"] = "lon lat"
        xds["velocity"].attrs["grid_mapping"] = "crs"
        # xds['velocity'] = xds['velocity'].astype(np.float64)

        # heading
        if "heading" in xds:
            xds["heading"].attrs["valid_range"] = [0, 3600]
            xds["heading"].attrs["standard_name"] = "direction_of_radial_vector_away_from_instrument"
            xds["heading"].attrs["units"] = "degrees"
            xds["heading"].attrs["coordinates"] = "lon lat"
            xds["heading"].attrs["scale_factor"] = 0.1
            xds["heading"].attrs["grid_mapping"] = "crs"
            # xds['heading'] = xds['heading'].astype(np.float64)

        # vector_flag
        if "vector_flag" in xds:
            xds["vector_flag"].attrs["long_name"] = "Vector Flag Masks"
            xds["vector_flag"].attrs["valid_range"] = [0, 2048]
            xds["vector_flag"].attrs["flag_masks"] = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
            xds["vector_flag"].attrs[
                "flag_meanings"
            ] = "grid_point_deleted grid_point_near_coast point_measurement no_radial_solution baseline_interpolation" \
                "exceeds_max_speed invalid_solution solution_beyond_valid_spatial_domain insufficient_angular_resolution" \
                "reserved reserved"
            xds["vector_flag"].attrs["coordinates"] = "lon lat"
            xds["vector_flag"].attrs["grid_mapping"] = "crs"
            # xds['vector_flag'] = xds['vector_flag'].astype(np.float64)

        # spatial_quality
        if "spatial_quality" in xds:
            xds["spatial_quality"].attrs["long_name"] = "Spatial Quality of radial sea water velocity"
            xds["spatial_quality"].attrs["units"] = "cm s-1"
            xds["spatial_quality"].attrs["coordinates"] = "lon lat"
            xds["spatial_quality"].attrs["grid_mapping"] = "crs"
            # xds['spatial_quality'] = xds['spatial_quality'].astype(np.float64)

        # temporal_quality
        if "temporal_quality" in xds:
            xds["temporal_quality"].attrs["long_name"] = "Temporal Quality of radial sea water velocity"
            xds["temporal_quality"].attrs["units"] = "cm s-1"
            xds["temporal_quality"].attrs["coordinates"] = "lon lat"
            xds["temporal_quality"].attrs["grid_mapping"] = "crs"
            # xds['temporal_quality'] = xds['temporal_quality'].astype(np.float64)

        # velocity_max
        if "velocity_max" in xds:
            xds["velocity_max"].attrs["long_name"] = "Maximum Velocity of sea water (away from instrument)"
            xds["velocity_max"].attrs["units"] = "cm s-1"
            xds["velocity_max"].attrs["coordinates"] = "lon lat"
            xds["velocity_max"].attrs["grid_mapping"] = "crs"
            # xds['velocity_max'] = xds['velocity_max'].astype(np.float64)

        # velocity_min
        if "velocity_min" in xds:
            xds["velocity_min"].attrs["long_name"] = "Minimum Velocity of sea water (away from instrument)"
            xds["velocity_min"].attrs["units"] = "cm s-1"
            xds["velocity_min"].attrs["coordinates"] = "lon lat"
            xds["velocity_min"].attrs["grid_mapping"] = "crs"
            # xds['velocity_min'] = xds['velocity_min'].astype(np.float64)

        # spatial_count
        if "spatial_count" in xds:
            xds["spatial_count"].attrs["long_name"] = "Spatial count of sea water velocity (away from instrument)"
            xds["spatial_count"].attrs["coordinates"] = "lon lat"
            xds["spatial_count"].attrs["grid_mapping"] = "crs"
            # xds['spatial_count'] = xds['spatial_count'].astype(np.float64)

        # temporal_count
        if "temporal_count" in xds:
            xds["temporal_count"].attrs["long_name"] = "Temporal count of sea water velocity (away from instrument)"
            xds["temporal_count"].attrs["coordinates"] = "lon lat"
            xds["temporal_count"].attrs["grid_mapping"] = "crs"
            # xds['temporal_count'] = xds['temporal_count'].astype(np.float64)

        # east_dist_from_origin
        if "dist_east_from_origin" in xds:
            xds["dist_east_from_origin"].attrs["long_name"] = "Eastward distance from instrument"
            xds["dist_east_from_origin"].attrs["units"] = "km"
            xds["dist_east_from_origin"].attrs["coordinates"] = "lon lat"
            xds["dist_east_from_origin"].attrs["grid_mapping"] = "crs"
            # xds['dist_east_from_origin'] = xds['dist_east_from_origin'].astype(np.float64)

        # north_dist_from_origin
        if "dist_north_from_origin" in xds:
            xds["dist_north_from_origin"].attrs["long_name"] = "Northward distance from instrument"
            xds["dist_north_from_origin"].attrs["units"] = "km"
            xds["dist_north_from_origin"].attrs["coordinates"] = "lon lat"
            xds["dist_north_from_origin"].attrs["grid_mapping"] = "crs"
            # xds['dist_north_from_origin'] = xds['dist_north_from_origin'].astype(np.float64)

        # range_cell
        if "range_cell" in xds:
            xds["range_cell"].attrs[
                "long_name"] = "Cross Spectra Range Cell  of sea water velocity (away from instrument)"
            xds["range_cell"].attrs["coordinates"] = "lon lat"
            xds["range_cell"].attrs["grid_mapping"] = "crs"
            # xds['range_cell'] = xds['range_cell'].astype(np.float64)

        # range_cell
        if "accuracy" in xds:
            xds["accuracy"].attrs["long_name"] = "Accuracy"
            xds["accuracy"].attrs["coordinates"] = "lon lat"
            xds["accuracy"].attrs["grid_mapping"] = "crs"
            xds["accuracy"].attrs["units"] = "cm s-1"
            # xds['accuracy'] = xds['accuracy'].astype(np.float64)

        # Q201
        if "Q201" in xds:
            xds["Q201"].attrs["long_name"] = "Syntax (QARTOD Test 201) Flag Masks"
            xds["Q201"].attrs["valid_range"] = [1, 9]
            xds["Q201"].attrs["flag_values"] = [1, 2, 3, 4, 5]
            xds["Q201"].attrs["flag_meanings"] = "pass not_evaluated suspect fail missing_data"
            xds["Q201"].attrs["coordinates"] = "lon lat"
            xds["Q201"].attrs["grid_mapping"] = "crs"
            # xds['Q201'] = xds['Q201'].astype(np.float64)
            rename_qc["Q201"] = "syntax_qc"

        # Q202
        if "Q202" in xds:
            xds["Q202"].attrs["long_name"] = "Maximum Velocity Threshold (QARTOD Test 202) Flag Masks"
            xds["Q202"].attrs["valid_range"] = [1, 9]
            xds["Q202"].attrs["flag_values"] = [1, 2, 3, 4, 5]
            xds["Q202"].attrs["flag_meanings"] = "pass not_evaluated suspect fail missing_data"
            xds["Q202"].attrs["coordinates"] = "lon lat"
            xds["Q202"].attrs["grid_mapping"] = "crs"
            # xds['Q202'] = xds['Q202'].astype(np.float64)
            rename_qc["Q202"] = "max_threshold_qc"

        # Q203
        if "Q203" in xds:
            xds["Q203"].attrs["long_name"] = "Valid Location (QARTOD Test 203) Flag Masks"
            xds["Q203"].attrs["valid_range"] = [1, 9]
            xds["Q203"].attrs["flag_values"] = [1, 2, 3, 4, 5]
            xds["Q203"].attrs["flag_meanings"] = "pass not_evaluated suspect fail missing_data"
            xds["Q203"].attrs["coordinates"] = "lon lat"
            xds["Q203"].attrs["grid_mapping"] = "crs"
            # xds['Q203'] = xds['Q203'].astype(np.float64)
            rename_qc["Q203"] = "valid_location_qc"

        # Q204
        if "Q204" in xds:
            xds["Q204"].attrs["long_name"] = "Radial Count (QARTOD Test 204) Flag Masks"
            xds["Q204"].attrs["valid_range"] = [1, 9]
            xds["Q204"].attrs["flag_values"] = [1, 2, 3, 4, 5]
            xds["Q204"].attrs["flag_meanings"] = "pass not_evaluated suspect fail missing_data"
            xds["Q204"].attrs["coordinates"] = "lon lat"
            xds["Q204"].attrs["grid_mapping"] = "crs"
            # xds['Q204'] = xds['Q204'].astype(np.float64)
            rename_qc["Q204"] = "radial_count_qc"

        # Q205
        if "Q205" in xds:
            xds["Q205"].attrs["long_name"] = "Spatial Median Filter (QARTOD Test 205) Flag Masks"
            xds["Q205"].attrs["valid_range"] = [1, 9]
            xds["Q205"].attrs["flag_values"] = [1, 2, 3, 4, 5]
            xds["Q205"].attrs["flag_meanings"] = "pass not_evaluated suspect fail missing_data"
            xds["Q205"].attrs["coordinates"] = "lon lat"
            xds["Q205"].attrs["grid_mapping"] = "crs"
            # xds['Q205'] = xds['Q205'].astype(np.float64)
            rename_qc["Q205"] = "spatial_median_filter_qc"

        # Q206
        if "Q206" in xds:
            xds["Q206"].attrs["long_name"] = "Temporal Gradient (QARTOD Test 206) Flag Masks"
            xds["Q206"].attrs["valid_range"] = [1, 9]
            xds["Q206"].attrs["flag_values"] = [1, 2, 3, 4, 5]
            xds["Q206"].attrs["flag_meanings"] = "pass not_evaluated suspect fail missing_data"
            xds["Q206"].attrs["coordinates"] = "lon lat"
            xds["Q206"].attrs["grid_mapping"] = "crs"
            # xds['Q206'] = xds['Q206'].astype(np.float64)
            rename_qc["Q206"] = "temporal_gradient_qc"

        # Q207
        if "Q207" in xds:
            xds["Q207"].attrs["long_name"] = "Average Radial Bearing (QARTOD Test 207) Flag Masks"
            xds["Q207"].attrs["valid_range"] = [1, 9]
            xds["Q207"].attrs["flag_values"] = [1, 2, 3, 4, 5]
            xds["Q207"].attrs["flag_meanings"] = "pass not_evaluated suspect fail missing_data"
            xds["Q207"].attrs["coordinates"] = "lon lat"
            xds["Q207"].attrs["grid_mapping"] = "crs"
            # xds['Q207'] = xds['Q207'].astype(np.float64)
            rename_qc["Q207"] = "average_radial_bearing_qc"

        # Q209
        if "Q209" in xds:
            xds["Q209"].attrs["long_name"] = "Radial Stuck Value (QARTOD Test 209) Flag Masks"
            xds["Q209"].attrs["valid_range"] = [1, 9]
            xds["Q209"].attrs["flag_values"] = [1, 2, 3, 4, 5]
            xds["Q209"].attrs["flag_meanings"] = "pass not_evaluated suspect fail missing_data"
            xds["Q209"].attrs["coordinates"] = "lon lat"
            xds["Q209"].attrs["grid_mapping"] = "crs"
            # xds['Q209'] = xds['Q209'].astype(np.float64)
            rename_qc["Q209"] = "radial_stuck_value_qc"

        # PRIM
        if "PRIM" in xds:
            xds["PRIM"].attrs["long_name"] = "Primary Flag Masks"
            xds["PRIM"].attrs["valid_range"] = [1, 9]
            xds["PRIM"].attrs["flag_values"] = [1, 2, 3, 4, 5]
            xds["PRIM"].attrs["flag_meanings"] = "pass not_evaluated suspect fail missing_data"
            xds["PRIM"].attrs["coordinates"] = "lon lat"
            xds["PRIM"].attrs["grid_mapping"] = "crs"
            # xds['PRIM'] = xds['PRIM'].astype(np.float64)
            rename_qc["PRIM"] = "primary_flag_qc"

        # QCOP
        if "QCOP" in xds:
            xds["QCOP"].attrs["long_name"] = "Operator Flag Masks"
            xds["QCOP"].attrs["valid_range"] = [1, 9]
            xds["QCOP"].attrs["flag_values"] = [1, 2, 3, 4, 5]
            xds["QCOP"].attrs["flag_meanings"] = "pass not_evaluated suspect fail missing_data"
            xds["QCOP"].attrs["coordinates"] = "lon lat"
            xds["QCOP"].attrs["grid_mapping"] = "crs"
            xds["QCOP"].attrs[
                "comment"
            ] = "Flag that is manually set by operator. Flag vectors that are not detected by QC tests but are clearly wrong."
            # xds['QCOP'] = xds['QCOP'].astype(np.float64)
            rename_qc["QCOP"] = "operator_flag_qc"

        # rename variables to something meaningful if they exist in the xarray dataset
        xds = xds.rename({k: v for k, v in rename_qc.items() if k in xds})
        # del xds.attrs['TimeStamp']

        return xds

    def clean_header(self, split_origin=False):
        """
        Clean up the radial header dictionary so that you can upload it to the HFR MySQL Database.

        Args:
            split_origin (bool, optional):
                Split the 'Origin' string for the receiver location into a list of floats [lon, lat]. Defaults to False.
        """
        keep = [
            "CTF",
            "FileType",
            "LLUVSpec",
            "UUID",
            "Manufacturer",
            "Site",
            "TimeStamp",
            "TimeZone",
            "TimeCoverage",
            "Origin",
            "GreatCircle",
            "GeodVersion",
            "LLUVTrustData",
            "RangeStart",
            "RangeEnd",
            "RangeResolutionKMeters",
            "AntennaBearing",
            "ReferenceBearing",
            "AngularResolution",
            "SpatialResolution",
            "PatternType",
            "PatternDate",
            "PatternResolution",
            "TransmitCenterFreqMHz",
            "DopplerResolutionHzPerBin",
            "FirstOrderMethod",
            "BraggSmoothingPoints",
            "CurrentVelocityLimits",
            "BraggHasSecondOrder",
            "RadialBraggPeakDropOff",
            "RadialBraggPeakNull",
            "RadialBraggNoiseThreshold",
            "PatternAmplitudeCorrections",
            "PatternPhaseCorrections",
            "PatternAmplitudeCalculations",
            "PatternPhaseCalculations",
            "RadialMusicParameters",
            "MergedCount",
            "RadialMinimumMergePoints",
            "FirstOrderCalc",
            "MergeMethod",
            "PatternMethod",
            "TransmitSweepRateHz",
            "TransmitBandwidthKHz",
            "SpectraRangeCells",
            "SpectraDopplerCells",
            "TableType",
            "TableColumns",
            "TableColumnTypes",
            "TableRows",
            "TableStart",
            "CurrentVelocityLimit",
        ]

        # TableColumnTypes
        key_list = list(self.metadata.keys())
        for key in key_list:
            if key not in keep:
                del self.metadata[key]

        for k, v in self.metadata.items():
            if "Site" in k:
                # WERA has lines like: '%Site: csw "CSW' and '%Site: gtn "gtn'
                # This should work for both CODAR and WERA files
                split_site = v.split(" ", 1)[0]
                self.metadata[k] = "".join(e for e in split_site if e.isalnum())
            elif k in ("TimeStamp", "PatternDate"):
                try:
                    t_list = [int(s) for s in v.split()]
                    self.metadata[k] = dt.datetime(*t_list)
                except ValueError:
                    # Can't parse a date, set to None
                    self.metadata[k] = None
            elif "TimeZone" in k:
                self.metadata[k] = v.split('"')[1]
            elif "TableColumnTypes" in k:
                self.metadata[k] = " ".join([x.strip() for x in v.strip().split(" ")])
            elif "Origin" in k:
                if split_origin:
                    self.metadata[k] = re.findall(r"[-+]?\d*\.\d+|\d+", v)
                else:
                    self.metadata[k] = v.lstrip()
            elif k in (
                    "RangeStart",
                    "RangeEnd",
                    "AntennaBearing",
                    "ReferenceBearing",
                    "AngularResolution",
                    "SpatialResolution",
                    "FirstOrderMethod",
                    "BraggSmoothingPoints",
                    "BraggHasSecondOrder",
                    "MergedCount",
                    "RadialMinimumMergePoints",
                    "FirstOrderCalc",
                    "SpectraRangeCells",
                    "SpectraDopplerCells",
                    "TableColumns",
                    "TableRows",
                    "PatternResolution",
                    "CurrentVelocityLimit",
                    "TimeCoverage",
            ):
                try:
                    self.metadata[k] = int(v)
                except ValueError:
                    temp = v.split(" ")[0]
                    try:
                        self.metadata[k] = int(temp)
                    except ValueError:
                        try:
                            self.metadata[k] = int(temp.split(".")[0])
                        except ValueError:
                            self.metadata[k] = None
            elif k in (
                    "RangeResolutionKMeters",
                    "CTF",
                    "TransmitCenterFreqMHz",
                    "DopplerResolutionHzPerBin",
                    "RadialBraggPeakDropOff",
                    "RadialBraggPeakNull",
                    "RadialBraggNoiseThreshold",
                    "TransmitSweepRateHz",
                    "TransmitBandwidthKHz",
            ):
                try:
                    self.metadata[k] = float(v)
                except ValueError:
                    try:
                        self.metadata[k] = float(v.split(" ")[0])
                    except ValueError:
                        self.metadata[k] = None
            else:
                continue

        required = ["Origin", "TransmitCenterFreqMHz"]
        present_keys = self.metadata.keys()
        for key in required:
            if key not in present_keys:
                self.metadata[key] = None

    def check_ehn_mandatory_variables(self):
        """
        This function checks if the Radial object contains all the mandatory data variables
        (i.e. not coordinate variables) required by the European standard data model developed in the framework of the
        EuroGOOS HFR Task Team.
        Missing variables are appended to the DataFrame containing data, filled with NaNs.

        INPUT:

        OUTPUT:
        """
        # Set mandatory variables based on the HFR manufacturer
        if self.is_wera:
            chkVars = ['VELU', 'VELV', 'VELO', 'HEAD', 'HCSS', 'EACC']
        else:
            chkVars = ['VELU', 'VELV', 'ESPC', 'ETMP', 'MAXV', 'MINV', 'ERSC', 'ERTC', 'XDST', 'YDST', 'VELO', 'HEAD',
                       'SPRC']

        # Check variables and add missing ones
        for vv in chkVars:
            if vv not in self.data.columns:
                self.data[vv] = np.nan

        return

    def apply_ehn_datamodel(self, network_data, station_data, version):
        """
        This function applies the European standard data model developed in the
        framework of the EuroGOOS HFR Task Team to the Radial object.
        The Radial object content is stored into an xarray Dataset built from the
        xarray DataArrays created by the Radial method to_xarray_multidimensional.
        Variable data types and data packing information are collected from
        "Data_Models/EHN/Radials/Radial_Data_Packing.json" file.
        Variable attribute schema is collected from
        "Data_Models/EHN/Radials/Radial_Variables.json" file.
        Global attribute schema is collected from
        "Data_Models/EHN/Global_Attributes.json" file.
        Global attributes are created starting from Radial object metadata and from
        DataFrames containing the information about HFR network and radial station
        read from the EU HFR NODE database.
        The generated xarray Dataset is attached to the Radial object, named as xds.

        INPUT:
            network_data: DataFrame containing the information of the network to which the radial site belongs
            station_data: DataFrame containing the information of the radial site that produced the radial
            version: version of the data model


        OUTPUT:
        """
        # Set the netCDF format
        ncFormat = 'NETCDF4_CLASSIC'

        # Expand Radial object variables along the coordinate axes
        self.to_xarray_multidimensional()

        # Set auxiliary coordinate sizes
        maxsiteSize = 1
        refmaxSize = 10
        maxinstSize = 1

        # Get data packing information per variable
        f = open('Data_Models/EHN/Radials/Radial_Data_Packing.json')
        dataPacking = json.loads(f.read())
        f.close()

        # Get variable attributes
        f = open('Data_Models/EHN/Radials/Radial_Variables.json')
        radVariables = json.loads(f.read())
        f.close()

        # Get global attributes
        f = open('Data_Models/EHN/Global_Attributes.json')
        globalAttributes = json.loads(f.read())
        f.close()

        # Rename velocity related variables
        self.xdr['DRVA'] = self.xdr.pop('HEAD')
        self.xdr['RDVA'] = self.xdr.pop('VELO')
        self.xdr['EWCT'] = self.xdr.pop('VELU')
        self.xdr['NSCT'] = self.xdr.pop('VELV')

        # Drop unnecessary DataArrays from the DataSet
        if not self.is_wera:
            self.xdr.pop('VFLG')
        toDrop = []
        for vv in self.xdr:
            if vv not in radVariables.keys():
                toDrop.append(vv)
        for rv in toDrop:
            self.xdr.pop(rv)

        # Add coordinate reference system to the dictionary
        self.xdr['crs'] = xr.DataArray(int(0), )

        # Add antenna related variables to the dictionary
        # Number of antennas
        self.xdr['NARX'] = xr.DataArray([np.expand_dims(station_data.iloc[0]['number_of_receive_antennas'], axis=0)],
                                        dims={'TIME': len(pd.date_range(self.time, periods=1)), 'MAXSITE': maxsiteSize})
        self.xdr['NATX'] = xr.DataArray([np.expand_dims(station_data.iloc[0]['number_of_transmit_antennas'], axis=0)],
                                        dims={'TIME': len(pd.date_range(self.time, periods=1)), 'MAXSITE': maxsiteSize})

        # Longitude and latitude of antennas
        self.xdr['SLTR'] = xr.DataArray([np.expand_dims(station_data.iloc[0]['site_lat'], axis=0)],
                                        dims={'TIME': len(pd.date_range(self.time, periods=1)), 'MAXSITE': maxsiteSize})
        self.xdr['SLNR'] = xr.DataArray([np.expand_dims(station_data.iloc[0]['site_lon'], axis=0)],
                                        dims={'TIME': len(pd.date_range(self.time, periods=1)), 'MAXSITE': maxsiteSize})
        self.xdr['SLTT'] = xr.DataArray([np.expand_dims(station_data.iloc[0]['site_lat'], axis=0)],
                                        dims={'TIME': len(pd.date_range(self.time, periods=1)), 'MAXSITE': maxsiteSize})
        self.xdr['SLNT'] = xr.DataArray([np.expand_dims(station_data.iloc[0]['site_lon'], axis=0)],
                                        dims={'TIME': len(pd.date_range(self.time, periods=1)), 'MAXSITE': maxsiteSize})

        # Codes of antennas
        antCode = ('%s' % station_data.iloc[0]['station_id']).encode()
        self.xdr['SCDR'] = xr.DataArray(np.array([[antCode]]),
                                        dims={'TIME': len(pd.date_range(self.time, periods=1)), 'MAXSITE': maxsiteSize})
        self.xdr['SCDR'].encoding['char_dim_name'] = 'STRING' + str(len(antCode))
        self.xdr['SCDT'] = xr.DataArray(np.array([[antCode]]),
                                        dims={'TIME': len(pd.date_range(self.time, periods=1)), 'MAXSITE': maxsiteSize})
        self.xdr['SCDT'].encoding['char_dim_name'] = 'STRING' + str(len(antCode))

        # Add SDN namespace variables to the dictionary
        siteCode = ('%s' % station_data.iloc[0]['network_id']).encode()
        self.xdr['SDN_CRUISE'] = xr.DataArray([siteCode], dims={'TIME': len(pd.date_range(self.time, periods=1))})
        self.xdr['SDN_CRUISE'].encoding['char_dim_name'] = 'STRING' + str(len(siteCode))
        platformCode = ('%s' % station_data.iloc[0]['network_id'] + '-' + station_data.iloc[0]['station_id']).encode()
        self.xdr['SDN_STATION'] = xr.DataArray([platformCode], dims={'TIME': len(pd.date_range(self.time, periods=1))})
        self.xdr['SDN_STATION'].encoding['char_dim_name'] = 'STRING' + str(len(platformCode))
        ID = ('%s' % platformCode.decode() + '_' + self.time.strftime('%Y-%m-%dT%H:%M:%SZ')).encode()
        self.xdr['SDN_LOCAL_CDI_ID'] = xr.DataArray([ID], dims={'TIME': len(pd.date_range(self.time, periods=1))})
        self.xdr['SDN_LOCAL_CDI_ID'].encoding['char_dim_name'] = 'STRING' + str(len(ID))
        self.xdr['SDN_EDMO_CODE'] = xr.DataArray([np.expand_dims(station_data.iloc[0]['EDMO_code'], axis=0)],
                                                 dims={'TIME': len(pd.date_range(self.time, periods=1)),
                                                       'MAXINST': maxinstSize})
        sdnRef = ('%s' % network_data.iloc[0]['metadata_page']).encode()
        self.xdr['SDN_REFERENCES'] = xr.DataArray([sdnRef], dims={'TIME': len(pd.date_range(self.time, periods=1))})
        self.xdr['SDN_REFERENCES'].encoding['char_dim_name'] = 'STRING' + str(len(sdnRef))
        sdnXlink = (
                    '%s' % '<sdn_reference xlink:href=\"' + sdnRef.decode() + '\" xlink:role=\"\" xlink:type=\"URL\"/>').encode()
        self.xdr['SDN_XLINK'] = xr.DataArray(np.array([[sdnXlink]]),
                                             dims={'TIME': len(pd.date_range(self.time, periods=1)),
                                                   'REFMAX': refmaxSize})
        self.xdr['SDN_XLINK'].encoding['char_dim_name'] = 'STRING' + str(len(sdnXlink))

        # Add spatial and temporal coordinate QC variables (set to good data due to the nature of HFR system)
        self.xdr['TIME_QC'] = xr.DataArray([1], dims={'TIME': len(pd.date_range(self.time, periods=1))})
        self.xdr['POSITION_QC'] = self.xdr['QCflag'] * 0 + 1
        self.xdr['DEPTH_QC'] = xr.DataArray([1], dims={'TIME': len(pd.date_range(self.time, periods=1))})

        # Create DataSet from DataArrays
        self.xds = xr.Dataset(self.xdr)

        # Add data variable attributes to the DataSet
        for vv in self.xds:
            self.xds[vv].attrs = radVariables[vv]

        # Update QC variable attribute "comment" for inserting test thresholds and attribute "flag_values" for assigning the right data type
        for qcv in self.metadata['QCTest']:
            if qcv in self.xds:
                self.xds[qcv].attrs['comment'] = self.xds[qcv].attrs['comment'] + ' ' + self.metadata['QCTest'][qcv]
                self.xds[qcv].attrs['flag_values'] = list(
                    np.int_(self.xds[qcv].attrs['flag_values']).astype(dataPacking[qcv]['dtype']))
        for qcv in ['TIME_QC', 'POSITION_QC', 'DEPTH_QC']:
            if qcv in self.xds:
                self.xds[qcv].attrs['flag_values'] = list(
                    np.int_(self.xds[qcv].attrs['flag_values']).astype(dataPacking[qcv]['dtype']))

        # Add coordinate variable attributes to the DataSet
        for cc in self.xds.coords:
            self.xds[cc].attrs = radVariables[cc]

        # Remove axis attributes from LONGITUDE and LATITUDE variables for polar geometry radials (i.e. Codar radials)
        if not self.is_wera:
            self.xds.LONGITUDE.attrs.pop('axis')
            self.xds.LATITUDE.attrs.pop('axis')

        # Evaluate measurement maximum depth
        vertMax = 3e8 / (8 * np.pi * station_data.iloc[0]['transmit_central_frequency'] * 1e6)

        # Evaluate time coverage start, end, resolution and duration
        timeCoverageStart = self.time - relativedelta(minutes=station_data.iloc[0]['temporal_resolution'] / 2)
        timeCoverageEnd = self.time + relativedelta(minutes=station_data.iloc[0]['temporal_resolution'] / 2)
        timeResRD = relativedelta(minutes=station_data.iloc[0]['temporal_resolution'])
        timeCoverageResolution = 'PT'
        if timeResRD.hours != 0:
            timeCoverageResolution += str(int(timeResRD.hours)) + 'H'
        if timeResRD.minutes != 0:
            timeCoverageResolution += str(int(timeResRD.minutes)) + 'M'
        if timeResRD.seconds != 0:
            timeCoverageResolution += str(int(timeResRD.seconds)) + 'S'

            # Fill global attributes
        globalAttributes['site_code'] = siteCode.decode()
        globalAttributes['platform_code'] = platformCode.decode()
        if station_data.iloc[0]['oceanops_ref'] != None:
            globalAttributes['oceanops_ref'] = station_data.iloc[0]['oceanops_ref']
            globalAttributes['wmo_platform_code'] = station_data.iloc[0]['wmo_code']
            globalAttributes['wigos_id'] = station_data.iloc[0]['wigos_id']
        globalAttributes['doa_estimation_method'] = station_data.iloc[0]['DoA_estimation_method']
        globalAttributes['calibration_type'] = station_data.iloc[0]['calibration_type']
        globalAttributes['last_calibration_date'] = station_data.iloc[0]['last_calibration_date'].strftime(
            '%Y-%m-%dT%H:%M:%SZ')
        globalAttributes['calibration_link'] = station_data.iloc[0]['calibration_link']
        # globalAttributes['title'] = network_data.iloc[0]['title']
        globalAttributes['title'] = 'Near Real Time Surface Ocean Radial Velocity by ' + globalAttributes[
            'platform_code']
        globalAttributes['summary'] = station_data.iloc[0]['summary']
        globalAttributes['institution'] = station_data.iloc[0]['institution_name']
        globalAttributes['institution_edmo_code'] = str(station_data.iloc[0]['EDMO_code'])
        globalAttributes['institution_references'] = station_data.iloc[0]['institution_website']
        globalAttributes['id'] = ID.decode()
        globalAttributes['project'] = network_data.iloc[0]['project']
        globalAttributes['comment'] = network_data.iloc[0]['comment']
        globalAttributes['network'] = network_data.iloc[0]['network_name']
        globalAttributes['data_type'] = globalAttributes['data_type'].replace('current data', 'radial current data')
        globalAttributes['geospatial_lat_min'] = str(network_data.iloc[0]['geospatial_lat_min'])
        globalAttributes['geospatial_lat_max'] = str(network_data.iloc[0]['geospatial_lat_max'])
        if not self.is_wera:
            globalAttributes['geospatial_lat_resolution'] = str(
                np.abs(np.mean(np.diff(self.xds.LATITUDE.values[:, 0]))))
        else:
            globalAttributes['geospatial_lat_resolution'] = self.metadata['DGT']
        globalAttributes['geospatial_lon_min'] = str(network_data.iloc[0]['geospatial_lon_min'])
        globalAttributes['geospatial_lon_max'] = str(network_data.iloc[0]['geospatial_lon_max'])
        if not self.is_wera:
            globalAttributes['geospatial_lon_resolution'] = str(
                np.abs(np.mean(np.diff(self.xds.LONGITUDE.values[:, 0]))))
        else:
            globalAttributes['geospatial_lon_resolution'] = self.metadata['DGT']
        globalAttributes['geospatial_vertical_max'] = str(vertMax)
        globalAttributes['geospatial_vertical_resolution'] = str(vertMax)
        globalAttributes['time_coverage_start'] = timeCoverageStart.strftime('%Y-%m-%dT%H:%M:%SZ')
        globalAttributes['time_coverage_end'] = timeCoverageEnd.strftime('%Y-%m-%dT%H:%M:%SZ')
        globalAttributes['time_coverage_resolution'] = timeCoverageResolution
        globalAttributes['time_coverage_duration'] = timeCoverageResolution
        globalAttributes['area'] = network_data.iloc[0]['area']
        globalAttributes['format_version'] = version
        globalAttributes['netcdf_format'] = ncFormat
        globalAttributes['citation'] += network_data.iloc[0]['citation_statement']
        globalAttributes['license'] = network_data.iloc[0]['license']
        globalAttributes['acknowledgment'] = network_data.iloc[0]['acknowledgment']
        globalAttributes['processing_level'] = '2B'
        globalAttributes['contributor_name'] = network_data.iloc[0]['contributor_name']
        globalAttributes['contributor_role'] = network_data.iloc[0]['contributor_role']
        globalAttributes['contributor_email'] = network_data.iloc[0]['contributor_email']
        globalAttributes['manufacturer'] = station_data.iloc[0]['manufacturer']
        globalAttributes['sensor_model'] = station_data.iloc[0]['manufacturer']
        globalAttributes['software_version'] = version

        creationDate = dt.datetime.utcnow()
        globalAttributes['metadata_date_stamp'] = creationDate.strftime('%Y-%m-%dT%H:%M:%SZ')
        globalAttributes['date_created'] = creationDate.strftime('%Y-%m-%dT%H:%M:%SZ')
        globalAttributes['date_modified'] = creationDate.strftime('%Y-%m-%dT%H:%M:%SZ')
        globalAttributes['history'] = 'Data measured at ' + self.time.strftime(
            '%Y-%m-%dT%H:%M:%SZ') + '. netCDF file created at ' \
                                      + creationDate.strftime('%Y-%m-%dT%H:%M:%SZ') + ' by the European HFR Node.'

        # Add global attributes to the DataSet
        self.xds.attrs = globalAttributes

        # Encode data types, data packing and _FillValue for the data variables of the DataSet
        for vv in self.xds:
            if vv in dataPacking:
                if 'dtype' in dataPacking[vv]:
                    self.xds[vv].encoding['dtype'] = dataPacking[vv]['dtype']
                if 'scale_factor' in dataPacking[vv]:
                    self.xds[vv].encoding['scale_factor'] = dataPacking[vv]['scale_factor']
                if 'add_offset' in dataPacking[vv]:
                    self.xds[vv].encoding['add_offset'] = dataPacking[vv]['add_offset']
                if 'fill_value' in dataPacking[vv]:
                    self.xds[vv].encoding['_FillValue'] = netCDF4.default_fillvals[
                        np.dtype(dataPacking[vv]['dtype']).kind + str(np.dtype(dataPacking[vv]['dtype']).itemsize)]
                else:
                    self.xds[vv].encoding['_FillValue'] = None

        # Update valid_min and valid_max variable attributes according to data packing
        for vv in self.xds:
            if 'valid_min' in radVariables[vv]:
                if ('scale_factor' in dataPacking[vv]) and ('add_offset' in dataPacking[vv]):
                    self.xds[vv].attrs['valid_min'] = np.float_(((radVariables[vv]['valid_min'] - dataPacking[vv][
                        'add_offset']) / dataPacking[vv]['scale_factor'])).astype(dataPacking[vv]['dtype'])
                else:
                    self.xds[vv].attrs['valid_min'] = np.float_(radVariables[vv]['valid_min']).astype(
                        dataPacking[vv]['dtype'])
            if 'valid_max' in radVariables[vv]:
                if ('scale_factor' in dataPacking[vv]) and ('add_offset' in dataPacking[vv]):
                    self.xds[vv].attrs['valid_max'] = np.float_(((radVariables[vv]['valid_max'] - dataPacking[vv][
                        'add_offset']) / dataPacking[vv]['scale_factor'])).astype(dataPacking[vv]['dtype'])
                else:
                    self.xds[vv].attrs['valid_max'] = np.float_(radVariables[vv]['valid_max']).astype(
                        dataPacking[vv]['dtype'])

        # Encode data types and avoid data packing, valid_min, valid_max and _FillValue for the coordinate variables of the DataSet
        for cc in self.xds.coords:
            if cc in dataPacking:
                if 'dtype' in dataPacking[cc]:
                    self.xds[cc].encoding['dtype'] = dataPacking[cc]['dtype']
                if 'valid_min' in radVariables[cc]:
                    del self.xds[cc].attrs['valid_min']
                if 'valid_max' in radVariables[cc]:
                    del self.xds[cc].attrs['valid_max']
                self.xds[cc].encoding['_FillValue'] = None

        return

    def to_netcdf(self, filename, model="tabular", prepend_extension=False, range_minmax=None, bearing=None,
                  enhance=True,
                  user_attributes=None):
        """
        Create a compressed netCDF4 (.nc) file from the Radial object

        Args:
            filename (str or Path): User defined filename of radial file you want to save
            model (str, optional): Create a 'tabular' (time) or 'gridded' (time, range, bearing) NetCDF. Defaults to 'tabular'
            prepend_extension (bool, optional):
                Prepend a descriptive term (bearing or gridded) to the .nc extension. Defaults to False.
            range_minmax (list or tuple, optional):
                For gridded NetCDF only, Range minimum (km) and range maximum (km) of the radial.
                Defaults to None which automatically infers the range based off of the radial data.
            bearing (list or tuple, optional):
                For gridded NetCDF only, Antenna bearing (degrees) and angular resolution (degrees) of the radial grid.
                Defaults to None which automatically infers the bearing based off of the radial data.
            enhance (bool, optional):
                Change manufacturer variables to meaningful variable names. Add attributes and other metadata. Defaults to True
            user_attributes (dictionary, optional): Dictionary containing metadata for the NetCDF. Defaults to None.
        """

        # Make a copy so that the original radial object is not altered by this function
        rcopy = copy.deepcopy(self)

        # Make sure filename is converted into a Path object
        filename = Path(filename)
        os.makedirs(filename.parent.resolve(), exist_ok=True)

        if not rcopy.is_valid():
            raise ValueError("Could not export ASCII data, the input file was invalid.")

        # If the filename does not have a .nc extension, we will add one.
        if ".nc" not in str(filename):
            filename = filename.with_suffix(".nc")

        # If the outputted file exists already, delete the existing file
        if os.path.isfile(filename):
            os.remove(filename)

        xds = rcopy.to_xarray(model, range_minmax=range_minmax, bearing=bearing, enhance=enhance)

        # if 'tabular' in netcdf_type:
        #     xds = self.to_xarray_tabular(enhance=enhance)
        # elif 'gridded' in netcdf_type:
        #     xds = self.to_xarray_gridded(range_minmax=range_minmax, bearing=bearing, enhance=enhance)

        # Check if dataset has distance_from_origin in coordinates. We will prepend the .nc extension
        # with the appropriate name depending on whether the radial file is tabular or gridded.
        if prepend_extension:
            if "bearing" in xds.coords:
                pre_ext = "gridded"
            else:
                pre_ext = "tabular"
            # Change the extension to reflect the type of wave file
            filename = filename.with_suffix(f".{pre_ext}.nc")

        # encoding = make_encoding(xds, comp_level=4, fillvalue=np.nan)

        # # Assign header data to global attributes
        # xds['site'] = self.metadata['Site'].strip('"').strip()
        # xds['site'] = xds['site'].assign_attrs(self.metadata)

        # Grab min and max time in dataset for entry into global attributes for cf compliance
        time_start = xds['time'].min().data
        time_end = xds['time'].max().data

        if user_attributes:
            from hfradarpy.io.nc import required_global_attributes
            global_attributes = required_global_attributes(user_attributes, time_start, time_end)

            global_attributes['geospatial_lat_min'] = np.double(xds.lat.min())
            global_attributes['geospatial_lat_max'] = np.double(xds.lat.max())
            global_attributes['geospatial_lon_min'] = np.double(xds.lon.min())
            global_attributes['geospatial_lon_max'] = np.double(xds.lon.max())

            logging.debug('{} - Assigning global attributes to dataset'.format(self.file_name))
            xds = xds.assign_attrs(global_attributes)

        # Encode and compress variables
        encoding = make_encoding(xds, comp_level=4, fillvalue=-999.0)

        if "gridded" in model:
            encoding["bearing"] = dict(zlib=False, _FillValue=None)
            encoding["range"] = dict(zlib=False, _FillValue=None)
        encoding["time"] = dict(zlib=False, _FillValue=None)

        xds.to_netcdf(filename, encoding=encoding, format="netCDF4", engine="netcdf4", unlimited_dims=["time"])


    def to_ruv(self, filename, validate=True, overwrite=False):
        """
        Create a CODAR Radial (.ruv) file from radial instance

        Args:
            filename (str or Path): User defined filename of radial file you want to save
            validate (boolean): If False, no validation check will be performed before creating the file.
            overwrite (bool): If True, an exported file can overwrite an existing file with the same name. Defaults to False.

        """
        # Make a copy so that the original radial object is not altered by this function
        rcopy = copy.deepcopy(self)

        # Make sure filename is converted into a Path object
        filename = Path(filename)

        if validate:
            if not rcopy.is_valid():
                raise ValueError("Could not export ASCII data, the input file was invalid.")

        # Ensure that the filename passed into the export function is not the same as the filename that we read in.
        # # We do not want to overwrite the original wave file by accident.
        if not overwrite:
            if rcopy.full_file == str(filename):
                suffix = f'.mod{filename.suffix}'
                filename = filename.with_suffix(suffix)

        if os.path.isfile(filename):
            os.remove(filename)

        os.makedirs(os.path.dirname(filename), exist_ok=True)
        with open(filename, "w") as f:
            # Write header
            for metadata_key, metadata_value in rcopy.metadata.items():
                if "ProcessedTimeStamp" in metadata_key:
                    break
                else:
                    # print(metadata_key)
                    f.write("%{}: {}\n".format(metadata_key, metadata_value))

            # Write data tables. Anything beyond the first table is commented out.
            for table in rcopy._tables.keys():

                if "datetime" in rcopy._tables[table]["data"].keys():
                    rcopy._tables[table]["data"] = rcopy._tables[table]["data"].drop(["datetime"], axis=1)

                for table_key, table_value in rcopy._tables[table].items():
                    if table_key != 'data':
                        if (table_key == 'TableType') & (table == 1):
                            if 'QCD' in rcopy.metadata:
                                for qcd_info in rcopy.metadata['QCD']:
                                    f.write('%{}\n'.format(qcd_info))
                            if 'QCTest' in rcopy.metadata:
                                f.write('%QCFileVersion: 2.0.0\n')
                                f.write(
                                    '%QCReference: Quality control reference: IOOS QARTOD HF Radar ver 2.0 June 2022\n')
                                f.write('%QCFlagDefinitions: 1=pass 2=not_evaluated 3=suspect 4=fail 9=missing_data\n')
                                f.write('%QCTestFormat: "test_name [qc_thresholds]: test_result"\n')

                                for test in rcopy.metadata["QCTest"].values():
                                    f.write("%QCTest: {}\n".format(test))
                            f.write("%{}: {}\n".format(table_key, table_value))
                        elif table_key == "TableColumns":
                            f.write("%TableColumns: {}\n".format(len(rcopy._tables[table]["data"].columns)))
                        elif table_key == "TableColumnTypes":
                            f.write("%TableColumnTypes: {}\n".format(
                                " ".join(rcopy._tables[table]["data"].columns.to_list())))
                        elif table_key == "TableStart":
                            f.write("%{}: {}\n".format(table_key, table_value))
                        elif table_key == "_TableHeader":
                            pass
                        else:
                            f.write("%{}: {}\n".format(table_key, table_value))

                if table == 1:
                    # Fill NaN with 999.000 which is the standard fill value for codar lluv files
                    rcopy.data = rcopy.data.fillna(999.000)

                    try:
                        rcopy.data["LOND"] = rcopy.data["LOND"].apply(lambda x: "{:.7f}".format(x))
                        rcopy.data["LATD"] = rcopy.data["LATD"].apply(lambda x: "{:.7f}".format(x))
                        if "ESPC" in rcopy.data.columns:
                            rcopy.data["ESPC"] = rcopy.data["ESPC"].apply(lambda x: "{:.3f}".format(x))
                        if "ETMP" in rcopy.data.columns:
                            rcopy.data["ETMP"] = rcopy.data["ETMP"].apply(lambda x: "{:.3f}".format(x))
                        rcopy.data["BEAR"] = rcopy.data["BEAR"].apply(lambda x: "{:.1f}".format(x))
                        rcopy.data["HEAD"] = rcopy.data["HEAD"].apply(lambda x: "{:.1f}".format(x))
                    except:
                        print("Unexpected error in formatting one of these columns: LOND LATD ESPC ETMP BEAR HEAD")

                    # Convert _TableHeader to a new dataframe and concatenate to dataframe containing radial data
                    # This allows for the output format to follow CODARS CTF specifications
                    # The below block of code adds the weird header and units format that codar uses in their files
                    row_df = pd.DataFrame([rcopy._tables[1]["_TableHeader"][1]],
                                          columns=rcopy._tables[1]["_TableHeader"][0])
                    rcopy.data.columns = rcopy._tables[1]["_TableHeader"][0]
                    rcopy.data = pd.concat([row_df, rcopy.data], ignore_index=True)
                    rcopy.data.insert(0, "%%", np.nan)  # Insert column at the beginning of dataframe of NaNs
                    rcopy.data.iloc[
                        0, rcopy.data.columns.get_loc("%%")] = "%%"  # make the first row in the first column a '%%'

                    # Output data table to string
                    # rcopy.data.to_string(f, index=False, justify='center', header=True, na_rep=' ')
                    rcopy.data.temp = re.sub(
                        " %%", "%%", rcopy.data.to_string(index=False, justify="right", header=True, na_rep=" ")
                    )
                    f.write(rcopy.data.temp)
                else:
                    table_alias = rcopy._tables[table]
                    table_alias["data"].insert(0, "%%", "%")
                    table_alias["data"] = table_alias["data"].fillna(999.000)
                    if table_alias["TableType"] == "rads rad1":
                        f.write(
                            "%%   Time       Calculated      Calculated     Corrected        Noise Floor     SignalToNoise   Diag  Valid  Dual  Radial RadsV Rads   Max   Vel    Vel    Bearing  Radial Spectra Time\n"
                        )
                        f.write(
                            "%% FromStart   Amp1    Amp2   Phase13 Phase23 Phase1 Phase2   NF1   NF2   NF3   SN1  SN2  SN3   Range Dopplr Angle Vector  per  Range Range  Max    Aver   Average   Type   Type   Year Mo Dy  Hr Mn  S\n"
                        )
                        f.write(
                            "%%  Seconds  (1/v^2) (1/v^2)   (deg)   (deg)   (deg)  (deg)  (dBm) (dBm) (dBm)  (dB) (dB) (dB)  Cell  Cells  Prcnt Count  Range Cells  (km) (cm/s) (cm/s) (deg CWN)\n"
                        )
                        table_alias["data"].to_string(f, index=False, justify="center", header=False)
                    elif table_alias["TableType"] == "rcvr rcv4":
                        f.write(
                            "%%   Minutes  Rcvr  Awg XmitTrip  AwgRun   Supply  +5VDC  -5VDC +12VDC XInt XAmp XForw XRefl  Xmit  X+Ampl  X+5V  X2Int X2Amp X2Forw X2Refl  Xmit2 X2+Amp  X2+5V  GpsRcv GpsDsp GpsSat GpsSat   PLL   HiRcvr  Humid Supply Extern Extern CompRunTime    Date                \n"
                        )
                        f.write(
                            "%%  FromStart degC degC  HexCode  Seconds   Volts  Volts  Volts  Volts degC degC Watts Watts  VSWR   Volts  Volts  degC  degC  Watts  Watts  VSWR   Volts  Volts   Mode   Mode   Lock  Unlock Unlock   degC     %    Amps  InputA InputB   Minutes      Year Mo Dy Hr Mn Sec\n"
                        )
                        table_alias["data"].to_string(f, index=False, justify="center", header=False)
                    elif table_alias["TableType"] == "rcvr rcv3":
                        f.write(
                            "%%   Minutes  Rcvr  Awg XmitTrip  AwgRun   Supply  +5VDC  -5VDC +12VDC XInt XAmp XForw XRefl  Xmit  X+Ampl  X+5V   GpsRcv GpsDsp GpsSat GpsSat   PLL   HiRcvr  Humid Supply Extern Extern CompRunTime    Date\n"
                        )
                        f.write(
                            "%%  FromStart degC degC  HexCode  Seconds   Volts  Volts  Volts  Volts degC degC Watts Watts  VSWR   Volts  Volts   Mode   Mode   Lock  Unlock Unlock   degC     %    Amps  InputA InputB   Minutes      Year Mo Dy Hr Mn Sec\n"
                        )
                        table_alias["data"].to_string(f, index=False, justify="center", header=False)
                    elif table_alias["TableType"] == "rcvr rcv2":
                        f.write(
                            "%%   LogTime Rcvr Awg3 XmtTrip Awg3Run Supply +5VDC -5VDC +12VDC XInt XAmp XForw XRefl X+Ampl X+5VDC GpsRcv GpsDsp  PLL   HiRcvr Humid Supply Extern Extern  CompRunTime   Year Mo Dy Hr Mn  S\n"
                        )
                        f.write(
                            "%%   Minutes   °C   °C HexCode Seconds Volts  Volts Volts  Volts   °C   °C Watts Watts  Volts  Volts  Mode   Mode  Unlock    °C     %   Amps  InputA InputB    Minutes                        \n"
                        )
                        table_alias["data"].to_string(f, index=False, justify="center", header=False)
                    else:
                        table_alias["data"].to_string(f, index=False, justify="center", header=True)

                if table > 1:
                    f.write("\n%TableEnd: {}\n".format(table))
                else:
                    f.write("\n%TableEnd: \n")
                f.write("%%\n")

            # Write footer containing processing information
            if 'ProcessedTimeStamp' in rcopy.metadata.keys():
                f.write("%ProcessedTimeStamp: {}\n".format(rcopy.metadata["ProcessedTimeStamp"]))
            for tool in rcopy.metadata["ProcessingTool"]:
                f.write("%ProcessingTool: {}\n".format(tool))
                # f.write('%{}: {}\n'.format(footer_key, footer_value))
            f.write("%End:")
            
    def file_type(self):
        """Return a string representing the type of file this is."""
        return 'radial'

    def initialize_qc(self):
        """
        Initialize dictionary entry for QC metadata.
        """
        # Initialize dictionary entry for QC metadta
        self.metadata['QCTest'] = {}

    # QARTOD QC TESTS

    def qc_qartod_avg_radial_bearing(self, reference_bearing, warning_threshold=15, failure_threshold=30):
        """
        Integrated Ocean Observing System (IOOS)
        Quality Assurance of Real-Time Oceanographic Data (QARTOD)
        Valid Location (Test 12)
        Check that the average radial bearing remains relatively constant (Roarty et al. 2012).

        It is expected that the average of all radial velocity bearings AVG_RAD_BEAR obtained during a sample
        interval (e.g., 1 hour) should be close to a reference bearing REF_RAD_BEAR and not vary beyond warning
        or failure thresholds.

        Args:
            reference_bearing (int or float): Reference bearing
            warning_threshold (int, optional): Warning Threshold. Defaults to 15.
            failure_threshold (int, optional): Failure Threshold. Defaults to 30.
        """
        test_str = "Q207"
        # Absolute value of the difference between the bearing mean and reference bearing
        absolute_difference = np.abs(self.data["BEAR"].mean() - reference_bearing)

        if absolute_difference >= failure_threshold:
            flag = 4
        elif (absolute_difference >= warning_threshold) & (absolute_difference < failure_threshold):
            flag = 3
        elif absolute_difference < warning_threshold:
            flag = 1

        self.data[test_str] = flag  # Assign the flags to the column
        self.metadata['QCTest'][
            test_str] = f"qc_qartod_avg_radial_bearing ({test_str}) - Test applies to entire file. Thresholds=" \
                        + "[ " + f"reference_bearing={reference_bearing} (degrees) " \
                        + f"warning={warning_threshold} (degrees) " \
                        + f"failure={failure_threshold} (degrees) " \
                        + f"]: See result in column {test_str} below"
        self.append_to_tableheader(test_str, "(flag)")

    def qc_qartod_valid_location(self, use_mask=True, res='high', angseg=None):
        """
        Integrated Ocean Observing System (IOOS)
        Quality Assurance of Real-Time Oceanographic Data (QARTOD)
        Valid Location (Test 8)
        Removes radial vectors placed over land or in other unmeasureable areas

        Radial vector coordinates are checked against a reference file containing information about which locations
        are over land or in an unmeasurable area (for example, behind an island or point of land). Radials in these
        areas will be flagged with a code (FLOC) in the radial file (+128 in CODAR radial files) and are not included
        in total vector calculations.

        Link: https://ioos.noaa.gov/ioos-in-action/manual-real-time-quality-control-high-frequency-radar-surface-current-data/

        Args:
            use_mask (bool, optional): Use mask_over_land function in addition to manufacturers flags. Defaults to False.
            res (string, optional): either 'low' or 'high', Defaults to low.
            angseg (array of dictionaries, optional): Used to define areas (segments) for invalid locations.
                Each array element is a dictionary defined by the following keys:
                range_min (int) :
                    Minimum range given in kilometers
                range_max (int) :
                    Maximum range given in kilometers
                bearingTrue_start (int) :
                    Start bearing of a segment (The angular segment will start at this value given in
                    degrees True and it will extend in a CLOCKWISE direction until it reaches the stop bearing.)
                bearingTrue_stop (int) :
                    Stop bearing of a segment (The angular segment begins at the start bearing and extends CLOCKWISE
                    until it reaches this value given in degrees True.)
       """

        test_str = "Q203"
        flag_column = "VFLG"
        applied_test_str = ''
        success = 0  # start with none of the tests running successfully

        self.data[test_str] = 1  # add new column of passing values

        if flag_column in self.data:
            try:
                self.data.loc[(self.data[
                                   flag_column] == 128), test_str] = 4  # set to 4 where SeaSonde AngSeg is flagged (manufacturer)
                applied_test_str += f"({flag_column}==128)"
                success = 1
            except:
                logger.warning(f"qc_qartod_valid_location with VFLG = 128 did not run successfully")

        if use_mask:
            try:
                self.data.loc[
                    ~self.mask_over_land(res=res), test_str] = 4  # set to 4 where land is flagged (mask_over_land)
                applied_test_str += f"(land mask {res} res)"
                success = 1
            except:
                logger.warning(f"qc_qartod_valid_location hfradarpy land mask did not run successfully")

        if angseg:
            try:
                for seg in angseg:
                    # find locations of radials in specified range and bearing limits and assign fail flags
                    if seg['bearingTrue_stop'] - seg['bearingTrue_start'] > 0:
                        try:
                            self.data.loc[
                                (self.data["RNGE"] >= seg['range_min']) & (
                                        self.data["RNGE"] <= seg['range_max']) & (
                                        self.data["BEAR"] >= seg['bearingTrue_start']) & (
                                        self.data["BEAR"] <= seg['bearingTrue_stop']), test_str
                            ] = 4  # set to 4 where hfradarpy angseg sections are flagged
                            success = 1
                        except:
                            logger.warning(
                                f"qc_qartod_valid_location hfradarpy angseg: one segment not applicable or did not run successfully")

                    else:
                        try:
                            self.data.loc[
                                (self.data["RNGE"] >= seg['range_min']) & (
                                        self.data["RNGE"] <= seg['range_max']) & (
                                        self.data["BEAR"] >= seg['bearingTrue_start']) & (
                                        self.data["BEAR"] <= 360), test_str
                            ] = 4  # set to 4 where hfradarpy angseg sections are flagged
                            self.data.loc[
                                (self.data["RNGE"] >= seg['range_min']) & (
                                        self.data["RNGE"] <= seg['range_max']) & (
                                        self.data["BEAR"] >= 0) & (
                                        self.data["BEAR"] <= seg['bearingTrue_stop']), test_str
                            ] = 4  # set to 4 where hfradarpy angseg sections are flagged
                            success = 1
                        except:
                            logger.warning(
                                f"qc_qartod_valid_location hfradarpy angseg: one segment not applicable or did not run successfully")
                applied_test_str += f"(angseg)"
            except:
                logger.warning(f"qc_qartod_valid_location hfradarpy angseg not applicable or did not run successfully")

        self.metadata["QCTest"][test_str] = f"qc_qartod_valid_location ({test_str}) - Test applies to each row. Thresholds=[{applied_test_str}]: " \
            + f"See results in column {test_str} below"

        if success == 0:
            self.data[test_str] = 2  # add column of "not evaluated" flags if none of the test methods were successful
            logger.warning(
                f"qc_qartod_valid_location did not run, no {flag_column} column, land mask and angseg either not used or not successfully applied")

        self.append_to_tableheader(test_str, "(flag)")


    def qc_qartod_radial_count(self, min_count=150, low_count=300):
        """
        Integrated Ocean Observing System (IOOS)
        Quality Assurance of Real-Time Oceanographic Data (QARTOD)
        Radial Count (Test 9)
        Rejects radials in files with low radial counts (poor radial map coverage).

        The number of radials (RCNT) in a radial file must be above a threshold value RCNT_MIN to pass the test and
        above a value RC_LOW to not be considered suspect. If the number of radials is below the minimum level,
        it indicates a problem with data collection. In this case, the file should be rejected and none of the
        radials used for total vector processing.

        Link: https://ioos.noaa.gov/ioos-in-action/manual-real-time-quality-control-high-frequency-radar-surface-current-data/

        Args:
            min_count (int, optional):
                Minimum radial count threshold (failure) below which the file should be rejected. Defaults to 150.
            low_count (int, optional):
                Low radial count threshold (warning) below which the file should be considered suspect. Defaults to 300.
        """
        test_str = "Q204"
        column_flag = "VFLG"

        # If a vector flag is supplied by the vendor, subset by that first
        if column_flag in self.data:
            num_radials = len(self.data[self.data[column_flag] != 128])
        else:
            num_radials = len(self.data)

        if num_radials < min_count:
            radial_count_flag = 4
        elif (num_radials >= min_count) and (num_radials <= low_count):
            radial_count_flag = 3
        elif num_radials > low_count:
            radial_count_flag = 1

        self.data[test_str] = radial_count_flag
        self.metadata['QCTest'][
            test_str] = f"qc_qartod_radial_count ({test_str}) - Test applies to entire file. Thresholds=" \
                        + "[ " + f"failure={min_count} (radials) " \
                        + f"warning_num={low_count} (radials) " \
                        + f"<valid_radials={num_radials}> " \
                        + f"]:  See results in column {test_str} below"
        self.append_to_tableheader(test_str, "(flag)")

    def qc_qartod_maximum_velocity(self, max_speed=250, high_speed=150):
        """
        Integrated Ocean Observing System (IOOS)
        Quality Assurance of Real-Time Oceanographic Data (QARTOD)
        Max Threshold (Test 7)
        Ensures that a radial current speed is not unrealistically high.

        The maximum radial speed threshold (RSPDMAX) represents the maximum reasonable surface radial velocity
        for the given domain.

        Link: https://ioos.noaa.gov/ioos-in-action/manual-real-time-quality-control-high-frequency-radar-surface-current-data/

        Args:
            max_speed (int, optional):
                Maximum Radial Speed (cm/s). Radials beyond this speed will be flagged a failure. Defaults to 250
            high_speed (int, optional):
                High Radial Speed (cm/s). Radials between high and max speed will be flagged suspect. Defaults to 150
        """
        test_str = "Q202"

        self.data["VELO"] = self.data["VELO"].astype(float)  # make sure VELO is a float

        # Add new column to dataframe for test, and set every row as passing, 1, flag
        self.data[test_str] = 1

        # velocity is less than radial_max_speed but greater than radial_high_speed, set that row as a warning, 3, flag
        self.data.loc[
            (self.data["VELO"].abs() < max_speed) & (self.data["VELO"].abs() > high_speed), test_str
        ] = 3

        # if velocity is greater than radial_max_speed, set that row as a fail, 4, flag
        self.data.loc[(self.data["VELO"].abs() > max_speed), test_str] = 4

        self.metadata['QCTest'][
            test_str] = f"qc_qartod_maximum_velocity ({test_str}) - Test applies to each row. Thresholds=" \
                        + "[ " + f"high_vel={str(high_speed)} (cm/s) " \
                        + f"max_vel={str(max_speed)} (cm/s) " \
                        + f"]: See results in column {test_str} below"
        self.append_to_tableheader(test_str, "(flag)")

    def qc_qartod_spatial_median(
            self, smed_range_cell_limit=2.1, smed_angular_limit=10, smed_current_difference=30
    ):
        """
        Integrated Ocean Observing System (IOOS)
        Quality Assurance of Real-Time Oceanographic Data (QARTOD)
        Spatial Median (Test 10)

        Ensures that the radial velocity is not too different from nearby radial velocities.
        RV is the radial velocity
        NV is a set of radial velocities for neighboring radial cells
        (cells within radius of 'radial_smed_range_cell_limit' * Range Step (km)
        and whose vector bearing (angle of arrival at site)
        is also within 'radial_smed_angular_limit' degrees of the source vector's bearing)
        Required to pass the test: |RV - median(NV)| <= radial_smed_current_difference
        Link: https://ioos.noaa.gov/ioos-in-action/manual-real-time-quality-control-high-frequency-radar-surface-current-data/

        Args:
            smed_range_cell_limit (float, optional):
                Multiple of range step which depends on the radar type. Defaults to 2.1.
            smed_angular_limit (int, optional):
                Limit for number of degrees from source radial's bearing (degrees). Defaults to 10.
            smed_current_difference (int, optional):
                Current difference (cm/s). Defaults to 30.
        """
        test_str = "Q205"

        self.data[test_str] = 1
        try:
            Bstep = [float(s) for s in re.findall(r"-?\d+\.?\d*", self.metadata["AngularResolution"])]
            Bstep = Bstep[0]
            # Bstep = int(min(np.diff(np.unique(self.data['BEAR']))))  #use as backup method if other fails?

            RLim = int(round(smed_range_cell_limit))  # if not an integer will cause an error later on
            BLim = int(smed_angular_limit / Bstep)  # if not an integer will cause an error later on

            # convert bearing into bearing cell numbers
            adj = np.mod(min(self.data["BEAR"]), Bstep)
            Bcell = ((self.data["BEAR"] - adj) / Bstep) - 1
            Bcell = Bcell.astype(int)
            # Btable = np.column_stack((self.data['BEAR'], Bcell))  #only for debugging

            if "SPRC" in self.data.columns:
                Rcell = self.data["SPRC"]
            else:
                ranges = self.data['RNGE'].unique()
                ranges_sorted = sorted(ranges)
                intervals = pd.Series(ranges_sorted).diff()
                range_interval = float(intervals.mode().iloc[0])
                Rcell = (self.data['RNGE'] / range_interval).round().astype(int)

            # Rtable = np.column_stack((self.data['RNGE'], Rcell))   #only for debugging

            # place velocities into a matrix with rows defined as bearing cell# and columns as range cell#
            BRvel = np.zeros((int(360 / Bstep), max(Rcell) + 1), dtype=int) + np.nan
            BRind = np.zeros((int(360 / Bstep), max(Rcell) + 1), dtype=int) + np.nan

            for xx in range(len(self.data["VELO"])):
                BRvel[Bcell[xx]][Rcell[xx]] = self.data["VELO"][xx]
                BRind[Bcell[xx]][Rcell[xx]] = xx  # keep track of indices so easier to return to original format

            # deal with 359 to 0 transition in bearing by
            # repeating first BLim rows at the bottom and last BLim rows at the top
            # also pad ranges with NaNs by adding extra columns on the left and right of the array
            # this keeps the indexing for selecting the neighbors from breaking

            BRtemp = np.append(np.append(BRvel[-BLim:], BRvel, axis=0), BRvel[:BLim], axis=0)
            rangepad = np.zeros((BRtemp.shape[0], RLim), dtype=int) + np.nan
            BRpad = np.append(np.append(rangepad, BRtemp, axis=1), rangepad, axis=1)

            # calculate median of neighbors (neighbors include the point itself)
            BRmed = BRpad + np.nan  # initialize with an array of NaN
            for rr in range(RLim, BRvel.shape[1] + RLim):
                for bb in range(BLim, BRvel.shape[0] + BLim):
                    temp = BRpad[bb - BLim: bb + BLim + 1, rr - RLim: rr + RLim + 1]  # temp is the matrix of neighbors
                    BRmed[bb][rr] = np.nanmedian(temp)

            # now remove the padding from the array containing the median values
            BRmedtrim = BRmed[BLim:-BLim, RLim:-RLim]

            # calculate velocity minus median of neighbors
            # and put back into single column using the indices saved in BRind
            BRdiff = (
                    BRvel - BRmedtrim
            )  # velocity minus median of neighbors, test these values against current radial_smed_current_difference
            diffcol = self.data["RNGE"] + np.nan  # initialize a single column for the difference results
            for rr in range(BRdiff.shape[1]):
                for bb in range(BRdiff.shape[0]):
                    if not (np.isnan(BRind[bb][rr])):
                        diffcol[BRind[bb][rr]] = BRdiff[bb][rr]
            boolean = diffcol.abs() > smed_current_difference
            self.data[test_str] = self.data[test_str].where(~boolean, other=4)

        except:
            self.data[test_str] = 2  # in case of an error, add column of "not evaluated" flags
            logger.warning(f"qc_qartod_spatial_median did not run due to an error")


        self.metadata['QCTest'][
            test_str] = f"qc_qartod_spatial_median ({test_str}) - Test applies to each row. Thresholds=" \
                        + "[ " + f"range_cell_limit={str(smed_range_cell_limit)} (range cells) " \
                        + f"angular_limit={str(smed_angular_limit)} (degrees) " \
                        + f"current_difference={str(smed_current_difference)} (cm/s) " \
                        + f"]: See results in column {test_str} below"
        self.append_to_tableheader(test_str, "(flag)")

    def qc_qartod_syntax(self):
        """
        Integrated Ocean Observing System (IOOS)
        Quality Assurance of Real-Time Oceanographic Data (QARTOD)
        Syntax (Test 6)

        This test is required to be QARTOD compliant.

        A collection of tests ensuring proper formatting and existence of fields within a radial file.

        The radial file may be tested for proper parsing and content, for file format (hfrweralluv1.0, for example),
        site code, appropriate time stamp, site coordinates, antenna pattern type (measured or ideal, for DF
        systems), and internally consistent row/column specifications.

        ----------------------------------------------------------------------------------------------------------------------
        Fail: One or more fields are corrupt or contain invalid data, If “File Format” ≠ “hfrweralluv1.0”, flag = 4

        Pass: Applies for test pass condition.
        ----------------------------------------------------------------------------------------------------------------------
        Link: https://ioos.noaa.gov/ioos-in-action/manual-real-time-quality-control-high-frequency-radar-surface-current-data/
        """
        test_str = "Q201"

        i = 0

        # check for timestamp in filename
        result = re.search(r"\d{4}_\d{2}_\d{2}_\d{4}", self.file_name)
        if result:
            timestr = dt.datetime.strptime(result.group(), "%Y_%m_%d_%H%M")
            i = i + 1

        # Radial tables must not be empty
        if self.is_valid():
            i = i + 1

        # The following metadata must be defined.
        tmp = self.metadata
        if (
                tmp["FileType"] and tmp["Site"] and tmp["TimeStamp"] and tmp["Origin"] and tmp["PatternType"] and tmp[
            "TimeZone"]
        ):
            filetime = dt.datetime(*map(int, self.metadata["TimeStamp"].split()))
            i = i + 1

        # Filename timestamp must match the timestamp reported within the file.
        if timestr == filetime:
            i = i + 1

        # Radial data table columns stated must match the number of columns reported for each row
        if len(self._tables[1]["TableColumnTypes"].split()) == self.data.shape[1]:
            i = i + 1

        # Make sure site location is within range: -180 <= lon <= 180 & -90 <= lat <= 90
        latlon = re.findall(r"[-+]?\d*\.\d+|\d+", self.metadata["Origin"])
        if (-180 <= float(latlon[1]) <= 180) & (-90 <= float(latlon[0]) <= 90):
            i = i + 1

        if i == 6:
            syntax = 1
        else:
            syntax = 4

        self.data[test_str] = syntax
        self.metadata['QCTest'][
            test_str] = f"qc_qartod_syntax ({test_str}) - Test applies to entire file. Thresholds=[N/A]: See results in column {test_str}"
        self.append_to_tableheader(test_str, "(flag)")

    def qc_qartod_temporal_gradient(self, r0, gradient_temp_fail=54, gradient_temp_warn=36):
        """
        Integrated Ocean Observing System (IOOS)
        Quality Assurance of Real-Time Oceanographic Data (QARTOD)
        Temporal Gradient (Test 11)
        Checks for satisfactory temporal rate of change of radial components

        Test determines whether changes between successive radial velocity measurements at a particular range
        and bearing cell are within an acceptable range. GRADIENT_TEMP = |Rt-1 - Rt|

        Flags Condition Codable Instructions
        Fail = 4 The temporal change between successive radial velocities exceeds the gradient failure threshold.

        If GRADIENT_TEMP ≥ GRADIENT_TEMP_FAIL,
        flag = 4

        Suspect = 3 The temporal change between successive radial velocities is less than the gradient failure threshold but
        exceeds the gradient warn threshold.

        If GRADIENT_TEMP < GRADIENT_TEMP_FAIL & GRADIENT_TEMP ≥ GRADIENT_TEMP_WARN,
        flag = 3

        Pass = 1 The temporal change between successive radial velocities is less than the gradient warn threshold.

        If GRADIENT_TEMP < GRADIENT_TEMP_WARN,
        flag = 1

        Link: https://ioos.noaa.gov/ioos-in-action/manual-real-time-quality-control-high-frequency-radar-surface-current-data/

        Args:
            r0 (str or Path): Full path to the filename of the previous hourly radial.
            gradient_temp_fail (int, optional): Maximum Radial Speed (cm/s). Defaults to 54.
            gradient_temp_warn (int, optional): Warning Radial Speed (cm/s). Defaults to 36.
        """
        test_str = "Q206"
        # self.data[test_str] = data
        self.metadata['QCTest'][
            test_str] = f"qc_qartod_temporal_gradient ({test_str}) - Test applies to each row. Thresholds=" \
                        + "[ " + f"gradient_temp_warn={str(gradient_temp_warn)} (cm/s*hr) " \
                        + f"gradient_temp_fail={str(gradient_temp_fail)} (cm/s*hr) " \
                        + f"]: See results in column {test_str} below"
        self.append_to_tableheader(test_str, "(flag)")

        if os.path.exists(r0):
            r0 = Radial(r0)

            if r0.is_valid():

                merged = self.data.merge(r0.data, on=["LOND", "LATD"], how="left", suffixes=(None, "_x"),
                                         indicator="Exist")
                difference = (merged["VELO"] - merged["VELO_x"]).abs()

                # Add new column to dataframe for test, and set every row as passing, 1, flag
                self.data[test_str] = 1

                # If any point in the recent radial does not exist in the previous radial, set row as a not evaluated, 2, flag
                self.data.loc[merged['Exist'] == 'left_only', test_str] = 2

                # velocity is less than radial_max_speed but greater than radial_high_speed, set row as a warning, 3, flag
                self.data.loc[(difference < gradient_temp_fail) & (difference > gradient_temp_warn), test_str] = 3

                # if velocity is greater than radial_max_speed, set that row as a fail, 4, flag
                self.data.loc[(difference > gradient_temp_fail), test_str] = 4
            else:
                # Add new column to dataframe for test, and set every row as not_evaluated, 2, flag
                self.data[test_str] = 2
                logging.warning(
                    '{} is corrupt or contains no data. Setting column {} to not_evaluated flag'.format(r0, test_str))

        else:
            # Add new column to dataframe for test, and set every row as not_evaluated, 2, flag
            self.data[test_str] = 2
            logging.warning(
                "{} does not exist at specified location. Setting column {} to not_evaluated flag".format(r0, test_str)
            )

    def qc_qartod_stuck_value_version_2(self, resolution=0.01, N=3):
        """
        Integrated Ocean Observing System (IOOS)
        Quality Assurance of Real-Time Oceanographic Data (QARTOD)
        Radial Stuck Value (Test 9)
        Tests for repeating values in radial time series at a location

        If the temporal change between the most recent velocity and each of the previous N-1 radial velocities
        has not exceeded the resolution of the measurement, the most recent velocity is considered a stuck value
        and is assigned a fail flag.  Previous velocities are values obtained from N-1 successive time steps prior to the
        measurement that is under evaluation.

        Flags Condition Codable Instructions

        V is a set of velocities where V at time=0 is the most recent velocity.
        V = [Vt=-(N-1) ...  Vt=-1, Vt=0]

        If all (Vt=0 - [ Vt=-1 ... Vt=-(N-1)]) < R,
        flag = 4

        Pass = 1
        If any (Vt=0 - [ Vt=-1 ... Vt=-(N-1)]) >= R,
        flag = 1

        Link: https://ioos.noaa.gov/ioos-in-action/manual-real-time-quality-control-high-frequency-radar-surface-current-data/

        Args:
            resolution (int, optional): Radial velocity resolution (cm/s). Defaults to 0.01
            N (int, optional): Number of successive time steps to check. Defaults to 3.
        """
        test_str = "Q901"
        # self.data[test_str] = data
        self.metadata['QCTest'][
            test_str] = f"qc_qartod_radial_stuck_value_version_2 ({test_str}) - Test applies to each row. Thresholds=" \
                        + "[ " + f"stuck_value_resolution={str(resolution)} (cm/s) " \
                        + f"stuck_value_number_of_timesteps={str(N)}" \
                        + f"]: See results in column {test_str} below"
        self.append_to_tableheader(test_str, "(flag)")

        # create list of previous files
        i = 1
        f0 = self.full_file
        t0 = get_time(f0)
        t0_str = t0.strftime('%Y_%m_%d_%H')
        filelist = list()

        while i < N:
            prev_time = t0 - dt.timedelta(hours=i)
            prev_time_str = prev_time.strftime('%Y_%m_%d_%H')
            prev_file = f0.replace(t0_str, prev_time_str)
            filelist.append(prev_file)
            i += 1

        # Add new column to dataframe, and set every row as passing, 1, flag
        self.data[test_str] = 1

        # temporary object to add up results while looping through past files
        rtemp = copy.deepcopy(self)
        rtemp.data[test_str] = 0

        for f in filelist:

            if os.path.exists(f):
                r0 = Radial(f)
                if r0.is_valid():
                    merged = self.data.merge(r0.data, on=["LOND", "LATD"], how="left", suffixes=(None, "_x"),
                                             indicator="Exist")
                    difference = (merged["VELO"] - merged["VELO_x"]).abs()

                    # If any point in the recent radial does not exist in the previous radial, set row as 999
                    rtemp.data.loc[merged['Exist'] == 'left_only', test_str] = 999

                    # If velocity difference is less than speed resolution, that's a stuck value to add to the count
                    rtemp.data.loc[(difference < resolution), test_str] = rtemp.data.loc[
                                                                              (difference < resolution), test_str] + 1
                else:
                    # If any of the previous files are invalid, set every row as not_evaluated, 2, flag. (Return, no need to check other files.)
                    self.data[test_str] = 2
                    logging.warning(
                        '{} is corrupt or contains no data. Setting column {} to not_evaluated flag'.format(r0,
                                                                                                            test_str))
                    return
            else:
                # If any of the previous files do not exist, set every row as not_evaluated, 2, flag. (Return, no need to check other files.)
                self.data[test_str] = 2
                logging.warning(
                    "{} does not exist at specified location. Setting column {} to not_evaluated flag".format(r0,
                                                                                                              test_str)
                )
                return

        # If stuck value persisted for N-1 previous files, then set row as a failure, 4, flag
        self.data.loc[rtemp.data[test_str] == N - 1, test_str] = 4

        # If any points in the past radial files did not exist, set row as a not evaluated, 2, flag
        self.data.loc[rtemp.data[test_str] >= 999, test_str] = 2

    def qc_qartod_stuck_value(self, resolution=0.01, N=3):
        """
        Integrated Ocean Observing System (IOOS)
        Quality Assurance of Real-Time Oceanographic Data (QARTOD)
        Radial Stuck Value (Test 9)
        Tests for repeating values in radial time series at a location

        If the difference between successive velocities has not exceeded the resolution of the measurement
        for N time steps, the most recent velocity is considered a stuck value and is assigned a fail flag.

        Flags Condition Codable Instructions

        V is a set of velocities where V at time=0 is the most recent velocity.
        V = [Vt=-(N-1) … Vt=-1, Vt=0]

        IF MAX(ABS(DIFF(V)) < R,
        flag = 4

        Pass = 1
        IF MAX(ABS(DIFF(V)) >= R
        flag = 1

        Link: https://ioos.noaa.gov/ioos-in-action/manual-real-time-quality-control-high-frequency-radar-surface-current-data/

        Args:
            resolution (int, optional): Radial velocity resolution (cm/s). Defaults to 0.01
            N (int, optional): Number of successive time steps to check. Defaults to 3.
        """
        test_str = "Q209"
        # self.data[test_str] = data
        self.metadata['QCTest'][
            test_str] = f"qc_qartod_radial_stuck_value ({test_str}) - Test applies to each row. Thresholds=" \
                        + "[ " + f"stuck_value_resolution={str(resolution)} (cm/s) " \
                        + f"stuck_value_number_of_timesteps={str(N)}" \
                        + f"]: See results in column {test_str} below"
        self.append_to_tableheader(test_str, "(flag)")

        # set all results to passing
        result = np.full(self.data['VELO'].shape, 1)

        # create list of previous N files
        i = 1
        f0 = self.full_file
        t0 = get_time(f0)
        t0_str = t0.strftime('%Y_%m_%d_%H')
        filelist = list()
        filelist.append(f0)
        merged = copy.deepcopy(self.data[['LOND', 'LATD', 'VELO']])
        vlist = list()
        vlist.append('VELO')

        while i < N:
            prev_time = t0 - dt.timedelta(hours=i)
            prev_time_str = prev_time.strftime('%Y_%m_%d_%H')
            prev_file = f0.replace(t0_str, prev_time_str)
            filelist.append(prev_file)
            if os.path.isfile(prev_file):
                previous_r = Radial(prev_file)
                if previous_r.is_valid():
                    suf = "_" + str(i)
                    merged = merged.merge(previous_r.data[['LOND', 'LATD', 'VELO']], on=["LOND", "LATD"], how="left",
                                          suffixes=(None, suf),
                                          indicator=False)
                    # build list of merged velocity columns
                    vlist.append('VELO' + suf)
            i += 1

        if len(vlist) < N:
            # set all results to not evaluated
            result = np.full(self.data['VELO'].shape, 2)
        else:
            # pull out subset of data with only the velocity columns
            ts = merged[vlist]

            # ts_diff = ts.diff(axis=1)
            # ts_absdiff = ts.diff(axis=1).abs()
            # ts_absdiffmax = ts.diff(axis=1).abs().max(axis=1)
            is_stuck = ts.diff(axis=1).abs().max(axis=1) < resolution
            stuck_rows_index = ts.index[is_stuck]
            # If stuck value persisted for N-1 previous files, then set row as a failure, 4, flag
            result[stuck_rows_index] = 4

            # Need to handle lat/lon locations that do not contain data across all time steps. Set those locations to not evaluated, 2 flag
            # It's important to do this after the above process of finding stuck rows since some of those
            # failures may exist only due to presence of nans.
            has_nan = ts.isna().any(axis=1)
            nan_rows_index = ts.index[has_nan]
            result[nan_rows_index] = 2

        # Add new column to dataframe
        self.data[test_str] = result

    def qc_qartod_primary_flag(self, include=None):
        """
        A primary flag is a single flag set to the worst case of all QC flags within the data record.

        Args:
            include (list, optional):
                list of quality control tests which should be included in the primary flag.
                Defaults to None, which includes all tests.
        """
        test_str = "PRIM"

        # Set summary flag column all equal to 1
        self.data[test_str] = 1

        # generate dictionary of executed qc tests found in the header
        executed = dict()
        for b in [x.split("-")[0].strip() for x in self.metadata["QCTest"].values()]:
            i = b.split(" ")
            executed[i[0]] = re.sub(r"[()]", "", i[1])

        if include:
            # only add qartod tests which were set by user to executed dictionary
            included_tests = list({key: value for key, value in executed.items() if key in include}.values())
        else:
            included_tests = list(executed.values())

        equals_3 = self.data[included_tests].eq(3).any(axis=1)
        self.data[test_str] = self.data[test_str].where(~equals_3, other=3)

        equals_4 = self.data[included_tests].eq(4).any(axis=1)
        self.data[test_str] = self.data[test_str].where(~equals_4, other=4)

        included_test_strs = ", ".join(included_tests)
        self.metadata['QCTest'][
            test_str] = f'qc_qartod_primary_flag ({test_str}) - Primary Flag - Highest flag value of {included_test_strs}' + '("not_evaluated" flag results ignored)'
        self.append_to_tableheader(test_str, "(flag)")
        # %QCFlagDefinitions: 1=pass 2=not_evaluated 3=suspect 4=fail 9=missing_data

    def append_to_tableheader(self, test_string, test_unit):
        """
        Append key, value of metadata in the radial header and footer to the _TableHeader

        Args:
            test_string (str): key
            test_unit (str): value
        """
        self._tables[1]["_TableHeader"][0].append(test_string)
        self._tables[1]["_TableHeader"][1].append(test_unit)

    # EUROPEAN HFR NODE (EHN) QC Tests

    def qc_ehn_avg_radial_bearing(self, minBear=0, maxBear=360):
        """

        This QC test labels radial velocity vectors with a ‘good_data” flag if the average radial bearing
        of all the vectors contained in the radial lie within the specified range for normal operations.
        Otherwise, the vectors are labeled with a “bad_data” flag.
        The ARGO QC flagging scale is used.

        This test is applicable only to DF systems.
        Data files from BF systems will have this variable filled with a “good_data” flag.

        This test was defined in the framework of the EuroGOOS HFR Task Team based on the
        Average Radial Bearing test (QC207) from the Integrated Ocean Observing System (IOOS)
        Quality Assurance of Real-Time Oceanographic Data (QARTOD).

        INPUTS:
            minBear: minimum angle in degrees of the specified range for normal operations
            maxBear: maximum angle in degrees of the specified range for normal operations


        """
        # Set the test name
        testName = 'AVRB_QC'

        # Evaluate the average bearing
        if self.is_wera:
            self.data.loc[:, testName] = 1
        else:
            avgBear = self.data['BEAR'].mean()

            if avgBear >= minBear and avgBear <= maxBear:
                self.data.loc[:, testName] = 1
            else:
                self.data.loc[:, testName] = 4

        self.metadata['QCTest'][testName] = 'Average Radial Bearing QC Test - Test applies to entire file. ' \
                                            + 'Thresholds=[' + f'minimum bearing={minBear} (degrees) - ' + f'maximum bearing={maxBear} (degrees)]'

    def qc_ehn_radial_count(self, radMinCount=150):
        """

        This test labels radial velocity vectors with a “good data” flag if the number of velocity vectors contained
        in the radial is bigger than the minimum count specified for normal operations.
        Otherwise the vectors are labeled with a “bad data” flag.
        The ARGO QC flagging scale is used.

        This test was defined in the framework of the EuroGOOS HFR Task Team based on the
        Radial Count test from (QC204) the Integrated Ocean Observing System (IOOS) Quality Assurance of
        Real-Time Oceanographic Data (QARTOD).

        INPUTS:
            radMinCount: minimum number of velocity vectors for normal operations

        """
        # Set the test name
        testName = 'RDCT_QC'

        # Get the number of velocity vectors contained in the radial
        numRad = len(self.data)

        if numRad >= radMinCount:
            self.data.loc[:, testName] = 1
        else:
            self.data.loc[:, testName] = 4

        self.metadata['QCTest'][testName] = 'Radial Count QC Test - Test applies to entire file. ' \
                                            + 'Threshold=[' + f'minimum number of radial vectors={radMinCount}]'

    def qc_ehn_maximum_velocity(self, radMaxSpeed=1.2):
        """
        This test labels radial velocity vectors whose module is smaller than a maximum velocity threshold
        with a “good data” flag. Otherwise the vectors are labeled with a “bad data” flag.
        The ARGO QC flagging scale is used.

        This test was defined in the framework of the EuroGOOS HFR Task Team based on the
        Max Threshold test (QC202) from the Integrated Ocean Observing System (IOOS) Quality Assurance of
        Real-Time Oceanographic Data (QARTOD).

        INPUTS:
            radMaxSpeed: maximum velocity in m/s for normal operations
        """
        # Set the test name
        testName = 'CSPD_QC'

        # Add new column to the DataFrame for QC data by setting every row as passing the test (flag = 1)
        self.data.loc[:, testName] = 1

        # set bad flag for velocities not passing the test
        if self.is_wera:
            self.data.loc[(self.data['VELO'].abs() > radMaxSpeed), testName] = 4  # velocity in m/s (CRAD)
        else:
            self.data.loc[(self.data['VELO'].abs() > radMaxSpeed * 100), testName] = 4  # velocity in cm/s (LLUV)

        self.metadata['QCTest'][testName] = 'Velocity Threshold QC Test - Test applies to each vector. ' \
                                            + 'Threshold=[' + f'maximum velocity={radMaxSpeed} (m/s)]'

    def qc_ehn_median_filter(self, dLim=10, curLim=0.5):
        """
        This test evaluates, for each radial vector, the median of the modules of all vectors lying
        within a distance of <dLim> km.
        Each vector for which the difference between its module and the median is smaller than
        the specified threshold for normal operations (curLim), is labeled with a "good data" flag.
        Otherwise the vector is labeled with a “bad data” flag.
        The ARGO QC flagging scale is used.
        The native CRS of the Radial is used for distance calculations.

        This test was defined in the framework of the EuroGOOS HFR Task Team based on the
        Spatial Median test (QC205) from the Integrated Ocean Observing System (IOOS) Quality
        Assurance of Real-Time Oceanographic Data (QARTOD).

        INPUTS:
            dLim: distance limit in km for selecting vectors for median calculation
            curLim: velocity-median difference threshold in m/s for normal opertions
        """
        # Set the test name
        testName = 'MDFL_QC'

        # Add new column to the DataFrame for QC data by setting every row as passing the test (flag = 1)
        self.data.loc[:, testName] = 1

        # Get the ellipsoid of the radial coordinate reference system
        radEllps = self.metadata['GreatCircle'].split()[0].replace('"', '')

        # Create Geod object with the retrieved ellipsoid
        g = Geod(ellps=radEllps)

        # Evaluate the median of velocities within dLim for each radial vector
        median = self.data.loc[:, ['LOND', 'LATD']].apply(lambda x: velocityMedianInDistLimits(x, self.data, dLim, g),
                                                          axis=1)

        # set bad flag for vectors not passing the test
        if self.is_wera:
            self.data.loc[abs(self.data['VELO'] - median) > curLim, testName] = 4  # velocity in m/s (CRAD)
        else:
            self.data.loc[abs(self.data['VELO'] - median) > curLim * 100, testName] = 4  # velocity in cm/s (LLUV)

        self.metadata['QCTest'][testName] = 'Median Filter QC Test - Test applies to each vector. ' \
                                            + 'Thresholds=[' + f'distance limit={str(dLim)} (km) ' + f'velocity-median difference threshold={str(curLim)} (m/s)]'

    def qc_ehn_over_water(self):
        """
        This test labels radial velocity vectors that lie on water with a good data” flag.
        Otherwise the vectors are labeled with a “bad data” flag.
        The ARGO QC flagging scale is used.

        Radial vector coordinates are checked against a reference file containing information
        about which locations are over land or in an unmeasurable area (for example, behind an
        island or point of land).
        The GeoPandas "naturalearth_lowres" is used as reference.
        CODAR radials are beforehand flagged based on the "VFLAG" variable (+128 means "on land").

        This test was defined in the framework of the EuroGOOS HFR Task Team based on the
        Valid Location (Test 203) from the Integrated Ocean Observing System (IOOS) Quality
        Assurance of Real-Time Oceanographic Data (QARTOD).
        """
        # Set the test name
        testName = 'OWTR_QC'

        # Add new column to the DataFrame for QC data by setting every row as passing the test (flag = 1)
        self.data.loc[:, testName] = 1

        # Set the "vector-on-land" flag name for CODAR data
        vectorOnLandFlag = 'VFLG'

        # Set bad flag where land is flagged by CODAR manufacturer
        if not self.is_wera:
            if vectorOnLandFlag in self.data:
                self.data.loc[(self.data[vectorOnLandFlag] == 128), testName] = 4

                # Set bad flag where land is flagged (mask_over_land method)
        self.data.loc[~self.mask_over_land(subset=False), testName] = 4

        self.metadata['QCTest'][testName] = 'Over Water QC Test - Test applies to each vector. ' \
                                            + 'Thresholds=[' + 'GeoPandas "naturalearth_lowres"]'

    def qc_ehn_maximum_variance(self, radMaxVar=1):
        """
        This test labels radial velocity vectors whose temporal variance is smaller than
        a maximum variance threshold with a “good data” flag.
        Otherwise the vectors are labeled with a “bad data” flag.
        The ARGO QC flagging scale is used.

        This test was defined in the framework of the EuroGOOS HFR Task Team based on the
        U Component Uncertainty and V Component Uncertainty tests (QC306 and QC307) from the
        Integrated Ocean Observing System (IOOS) Quality Assurance of Real-Time Oceanographic
        Data (QARTOD).

        This test is NOT RECOMMENDED for CODAR data because the parameter defining the variance
        is computed at each time step, and therefore considered not statistically solid
        (as documented in the fall 2013 CODAR Currents Newsletter).

        INPUTS:
            radMaxVar: maximum variance in m2/s2 for normal operations
        """
        # Set the test name
        testName = 'VART_QC'

        # Add new column to the DataFrame for QC data by setting every row as passing the test (flag = 1)
        self.data.loc[:, testName] = 1

        # Set bad flag for variances not passing the test
        if self.is_wera:
            self.data.loc[(self.data['HCSS'] > radMaxVar), testName] = 4  # HCSS is the temporal variance for WERA data
        else:
            self.data.loc[((self.data[
                                'ETMP'] / 100) ** 2 > radMaxVar), testName] = 4  # ETMP is the temporal standard deviation in cm/s for CODAR data

        self.metadata['QCTest'][testName] = 'Variance Threshold QC Test - Test applies to each vector. ' \
                                            + 'Threshold=[' + f'maximum variance={radMaxVar} (m2/s2)]'

    def qc_ehn_temporal_derivative(self, r0, tempDerThr=1):
        """
        This test compares the velocity of each radial vector with the velocity of the radial vector
        measured in the previous timestamp at the same location.
        Each vector for which the velocity difference is smaller than the specified threshold for normal
        operations (tempDerThr), is labeled with a "good data" flag.
        Otherwise the vector is labeled with a “bad data” flag.
        The ARGO QC flagging scale is used.

        This test was defined in the framework of the EuroGOOS HFR Task Team based on the
        Temporal Gradient test (QC206) from the Integrated Ocean Observing System (IOOS) Quality
        Assurance of Real-Time Oceanographic Data (QARTOD).

        INPUTS:
            r0: Radial object of the previous timestamp
            tempDerThr: velocity difference threshold in m/s for normal operations
        """
        # Set the test name
        testName = 'VART_QC'

        # Check if the previous timestamp radial file exists
        if not r0 is None:
            # Merge the data DataFrame of the two Radials and evaluate velocity differences at each location
            mergedDF = self.data.merge(r0.data, on=['LOND', 'LATD'], how='left', suffixes=(None, '_x'),
                                       indicator='Exist')
            velDiff = (mergedDF['VELO'] - mergedDF['VELO_x']).abs()

            # Add new column to the DataFrame for QC data by setting every row as passing the test (flag = 1)
            self.data.loc[:, testName] = 1

            # Set rows of the DataFrame for QC data as not evaluated (flag = 0) for locations existing in the current radial but not in the previous one
            self.data.loc[mergedDF['Exist'] == 'left_only', testName] = 0

            # Set bad flag for vectors not passing the test
            if self.is_wera:
                self.data.loc[(velDiff > tempDerThr), testName] = 4  # velocity in m/s (CRAD)
            else:
                self.data.loc[(velDiff > tempDerThr * 100), testName] = 4  # velocity in cm/s (LLUV)

        else:
            # Add new column to the DataFrame for QC data by setting every row as not evaluated (flag = 0)
            self.data.loc[:, testName] = 0

        self.metadata['QCTest'][testName] = 'Temporal Derivative QC Test - Test applies to each vector. ' \
                                            + 'Threshold=[' + f'velocity difference threshold={str(tempDerThr)} (m/s)]'

    def qc_ehn_overall_qc_flag(self):
        """

        This QC test labels radial velocity vectors with a ‘good_data” flag if all QC tests are passed.
        Otherwise, the vectors are labeled with a “bad_data” flag.
        The ARGO QC flagging scale is used.

        INPUTS:


        """
        # Set the test name
        testName = 'QCflag'

        # Add new column to the DataFrame for QC data by setting every row as not passing the test (flag = 4)
        self.data.loc[:, testName] = 4

        # Set good flags for vectors passing all QC tests
        self.data.loc[self.data.loc[:, self.data.columns.str.contains('_QC')].eq(1).all(axis=1), testName] = 1

        self.metadata['QCTest'][
            testName] = 'Overall QC Flag - Test applies to each vector. Test checks if all QC tests are passed.'

    def reset(self):
        """
        Reset data variable of object, r.data, back to original dataset
        """
        logging.info("Resetting instance data variable to original dataset")
        self.data = self._tables[1]["data"]


if __name__ == "__main__":
    from pathlib import Path

    data_path = (Path(__file__).parent.with_name("examples") / "data").resolve()
    f = data_path / "radials" / "ruv" / "SEAB" / "RDLi_SEAB_2019_01_01_0100.ruv"

    r = Radial(f, replace_invalid=True, vflip=True)
    # r.mask_over_land()
    ds = r.to_xarray('gridded')
    print(ds)
