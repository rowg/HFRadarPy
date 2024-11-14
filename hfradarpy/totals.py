import logging
import datetime as dt
from dateutil.relativedelta import relativedelta
import math
import numpy as np
import xarray as xr
import netCDF4
import pandas as pd
from pyproj import Geod
from shapely.geometry import Point
import geopandas as gpd
import re
import io
import os
from pathlib import Path
from common import fileParser, addBoundingBoxMetadata
from collections import OrderedDict
from calc import true2mathAngle, dms2dd, evaluateGDOP, createLonLatGridFromBB, createLonLatGridFromBBwera, \
    createLonLatGridFromTopLeftPointWera
import json
import fnmatch
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


# def concatenate_totals(radial_list):
#     """
#     This function takes a list of radial files. Loads them all separately using the Radial object and then combines
#     them along the time dimension using xarrays built-in concatenation routines.
#     :param radial_list: list of radial files that you want to concatenate
#     :return: radial files concatenated into an xarray dataset by range, bearing, and time
#     """
#
#     totals_dict = {}
#     for each in sorted(radial_list):
#         total = Total(each, multi_dimensional=True)
#         totals_dict[total.file_name] = total.ds
#
#     ds = xr.concat(totals_dict.values(), 'time')
#     return ds

def buildUStotal(ts, pts, USxds, networkData, stationData):
    """
    This function builds the Total object for the input timestamp from an xarray DataSet
    containing the gridded total data of a US network.

    INPUT:
        ts: timestamp as datetime object
        pts: Geoseries containing the lon/lat positions of the data geographical grid
        USxds: xarray DataSet containing gridded total data related to the input timestamp
        networkData: DataFrame containing the information of the network to which the total belongs
        stationData: DataFrame containing the information of the radial sites that produced the total

    OUTPUT:
        Tus: Total object containing the US network total data
    """

    #####
    # Create the Total object for the input timestamp
    #####

    # Create empty total with grid
    Tus = Total(grid=pts)

    # Add timestamp
    Tus.time = ts

    # Fill the Total object with data
    if USxds.u.units.strip() == 'm s-1':
        Tus.data['VELU'] = USxds.u.values[0, :, :].flatten() * 100
    elif USxds.u.units.strip() == 'cm s-1':
        Tus.data['VELU'] = USxds.u.values[0, :, :].flatten()
    if USxds.v.units.strip() == 'm s-1':
        Tus.data['VELV'] = USxds.v.values[0, :, :].flatten() * 100
    elif USxds.v.units.strip() == 'cm s-1':
        Tus.data['VELV'] = USxds.v.values[0, :, :].flatten()
    Tus.data['GDOP'] = USxds.hdop.values[0, :, :].flatten()
    Tus.data['NRAD'] = USxds.number_of_radials.values[0, :, :].flatten()
    Tus.data['VELO'] = np.sqrt(Tus.data['VELU'] ** 2 + Tus.data['VELV'] ** 2)
    Tus.data['HEAD'] = (360 + np.arctan2(Tus.data['VELU'], Tus.data['VELV']) * 180 / np.pi) % 360

    # Set is_wera attribute to False (all velocities are expressed in cm/s)
    Tus.is_wera = False

    # Set is_combined attribute to False
    Tus.is_combined = False

    # Get the indexes of rows without measurements (i.e. containing nan values)
    indexNames = Tus.data.loc[pd.isna(Tus.data['VELU']), :].index
    # Delete these row indexes from DataFrame
    Tus.data.drop(indexNames, inplace=True)
    Tus.data.reset_index(level=None, drop=False, inplace=True)

    # Get the list of contributing radial files
    contrRadFilesDF = pd.Series(data=USxds.radial_metadata.files_loaded.split('\n'))
    # Add contributing radial sites
    contrRadSites = contrRadFilesDF.apply(lambda x: x.split('_')[-5]).to_list()
    # Keep only radial sites registered on the datbase
    contrRadSites = stationData.loc[stationData['station_id'].isin(contrRadSites)]

    # Insert contributing radial sites into Total object
    Tus.site_source = pd.DataFrame(index=contrRadSites.index,
                                   columns=['#', 'Name', 'Lat', 'Lon', 'Coverage(s)', 'RngStep(km)', 'Pattern',
                                            'AntBearing(NCW)'])
    Tus.site_source['#'] = contrRadSites.index
    Tus.site_source['Name'] = contrRadSites.loc[:]['station_id']
    Tus.site_source['Lat'] = contrRadSites.loc[:]['site_lat']
    Tus.site_source['Lon'] = contrRadSites.loc[:]['site_lon']

    # Get spatial resolution from TDS url
    if 'Resolution' in networkData.loc[0]['TDS_root_url'].split('_'):
        sptRes = networkData.loc[0]['TDS_root_url'].split('_')[
            networkData.loc[0]['TDS_root_url'].split('_').index('Resolution') - 1]
        sptRes = re.sub(r"([0-9]+(\.[0-9]+)?)", r" \1 ", sptRes).strip()
    elif 'hourly' in networkData.loc[0]['TDS_root_url'].split('/'):
        sptRes = networkData.loc[0]['TDS_root_url'].split('/')[
            networkData.loc[0]['TDS_root_url'].split('/').index('hourly') - 1]
        sptRes = re.sub(r"([0-9]+(\.[0-9]+)?)", r" \1 ", sptRes).strip()
    else:
        sptRes = None

    # Add metadata
    Tus.metadata['TimeZone'] = '"UTC" +0.000 0 "GMT"'
    Tus.metadata['AveragingRadius'] = str(
        USxds.processing_parameters.grid_search_radius) + ' ' + USxds.processing_parameters.grid_search_radius_units
    Tus.metadata['CurrentVelocityLimit'] = str(
        USxds.processing_parameters.max_rtv_speed) + ' ' + USxds.processing_parameters.max_rtv_speed_units
    Tus.metadata['GridAxisOrientation'] = '0.0 DegNCW'
    if sptRes:
        Tus.metadata['GridSpacing'] = sptRes
    Tus = addBoundingBoxMetadata(Tus, USxds.attrs['geospatial_lon_min'], USxds.attrs['geospatial_lon_max'],
                                 USxds.attrs['geospatial_lat_min'], USxds.attrs['geospatial_lat_max'])

    return Tus


def convertEHNtoINSTACtotalDatamodel(tDS, networkData, stationData, version):
    """
    This function applies the Copernicus Marine Service data model to the input xarray
    dataset containing total (either non temporally aggregated or temporally aggregated)
    data. The input dataset must follow the European standard data model.
    Variable data types and data packing information are collected from
    "Data_Models/CMEMS_IN_SITU_TAC/Totals/Total_Data_Packing.json" file.
    Variable attribute schema is collected from
    "Data_Models/CMEMS_IN_SITU_TAC/Totals/Total_Variables.json" file.
    Global attribute schema is collected from
    "Data_Models/CMEMS_IN_SITU_TAC/Global_Attributes.json" file.
    Global attributes are created starting from the input dataset and from
    DataFrames containing the information about HFR network and radial station
    read from the EU HFR NODE database.
    The function returns an xarray dataset compliant with the Copernicus Marine Service
    In Situ TAC data model.

    INPUT:
        tDS: xarray DataSet containing total (either non temporally aggregated or temporally aggregated)
             data.
        networkData: DataFrame containing the information of the network to which the total belongs
        stationData: DataFrame containing the information of the radial sites that produced the total
        version: version of the data model

    OUTPUT:
        instacDS: xarray dataset compliant with the Copernicus Marine Service In Situ TAC data model
    """
    # Get data packing information per variable
    f = open('Data_Models/CMEMS_IN_SITU_TAC/Totals/Total_Data_Packing.json')
    dataPacking = json.loads(f.read())
    f.close()

    # Get variable attributes
    f = open('Data_Models/CMEMS_IN_SITU_TAC/Totals/Total_Variables.json')
    totVariables = json.loads(f.read())
    f.close()

    # Get global attributes
    f = open('Data_Models/CMEMS_IN_SITU_TAC/Global_Attributes.json')
    globalAttributes = json.loads(f.read())
    f.close()

    # Create the output dataset
    instacDS = tDS
    instacDS.encoding = {}

    # Evaluate time coverage start, end, resolution and duration
    timeCoverageStart = pd.Timestamp(instacDS['TIME'].values.min()).to_pydatetime() - relativedelta(
        minutes=networkData.iloc[0]['temporal_resolution'] / 2)
    timeCoverageEnd = pd.Timestamp(instacDS['TIME'].values.max()).to_pydatetime() + relativedelta(
        minutes=networkData.iloc[0]['temporal_resolution'] / 2)

    timeCoverageDuration = pd.Timedelta(timeCoverageEnd - timeCoverageStart).isoformat()

    # Build the file id
    ID = 'GL_TV_HF_' + tDS.attrs['platform_code'] + '_' + pd.Timestamp(
        instacDS['TIME'].values.max()).to_pydatetime().strftime('%Y%m%d')

    # Get the TIME variable values
    xdsTime = instacDS.TIME.values
    # Convert them to datetime datetimes
    dtTime = np.array([pd.Timestamp(t).to_pydatetime() for t in xdsTime])

    # Evaluate timestamp as number of days since 1950-01-01T00:00:00Z
    timeDelta = dtTime - dt.datetime.strptime('1950-01-01T00:00:00Z', '%Y-%m-%dT%H:%M:%SZ')
    ncTime = np.array([t.days + t.seconds / (60 * 60 * 24) for t in timeDelta])

    # Replace TIME variable values with timestamp as number of days since 1950-01-01T00:00:00Z
    instacDS = instacDS.assign_coords(TIME=ncTime)

    # Rename DEPTH_QC variable to DEPH_QC
    instacDS = instacDS.rename({'DEPTH_QC': 'DEPH_QC'})

    # Add DEPH variable
    instacDS['DEPH'] = xr.DataArray(0,
                                    dims={'DEPTH': 1},
                                    coords={'DEPTH': [0]})

    # Remove DEPTH variable
    instacDS = instacDS.drop_vars('DEPTH')

    # Remove crs variable (it's time-varying because of the temporal aggregation)
    instacDS = instacDS.drop_vars('crs')

    # Add time-independent crs variable
    instacDS['crs'] = xr.DataArray(int(0), )

    # Remove encoding for data variables
    for vv in instacDS:
        if 'char_dim_name' in instacDS[vv].encoding.keys():
            instacDS[vv].encoding = {'char_dim_name': instacDS[vv].encoding['char_dim_name']}
        else:
            instacDS[vv].encoding = {}

    # Remove attributes and encoding from coordinates variables
    for cc in instacDS.coords:
        instacDS[cc].attrs = {}
        instacDS[cc].encoding = {}

    # Add units and calendar to encoding for coordinate TIME
    instacDS['TIME'].encoding['units'] = 'days since 1950-01-01T00:00:00Z'
    instacDS['TIME'].encoding['calendar'] = 'standard'

    # Add data variable attributes to the DataSet
    for vv in instacDS:
        instacDS[vv].attrs = totVariables[vv]

    # Add coordinate variable attributes to the DataSet
    for cc in instacDS.coords:
        instacDS[cc].attrs = totVariables[cc]

    # Update QC variable attribute "comment" for inserting test thresholds
    for qcv in list(tDS.keys()):
        if 'QC' in qcv:
            if not qcv in ['TIME_QC', 'POSITION_QC', 'DEPTH_QC']:
                instacDS[qcv].attrs['comment'] = instacDS[qcv].attrs['comment'] + ' ' + tDS[qcv].attrs['comment']

                # Update QC variable attribute "flag_values" for assigning the right data type
    for qcv in instacDS:
        if 'QC' in qcv:
            instacDS[qcv].attrs['flag_values'] = list(
                np.int_(instacDS[qcv].attrs['flag_values']).astype(dataPacking[qcv]['dtype']))

    # Fill global attributes
    globalAttributes['site_code'] = tDS.attrs['site_code']
    globalAttributes['platform_code'] = tDS.attrs['platform_code']
    globalAttributes['platform_name'] = globalAttributes['platform_code']
    globalAttributes['doa_estimation_method'] = tDS.attrs['doa_estimation_method'].replace(',', ';')
    globalAttributes['calibration_type'] = tDS.attrs['calibration_type'].replace(',', ';')
    globalAttributes['last_calibration_date'] = tDS.attrs['last_calibration_date'].replace(',', ';')
    globalAttributes['calibration_link'] = tDS.attrs['calibration_link'].replace(',', ';')
    globalAttributes['summary'] = tDS.attrs['summary']
    globalAttributes['institution'] = tDS.attrs['institution'].replace(',', ';')
    globalAttributes['institution_edmo_code'] = tDS.attrs['institution_edmo_code'].replace(', ', ' ')
    globalAttributes['institution_references'] = tDS.attrs['institution_references'].replace(', ', ' ')
    # globalAttributes['institution_abreviated'] = tDS.attrs['institution_abreviated'].replace(', ',' ')
    # globalAttributes['institution_country'] = tDS.attrs['institution_country'].replace(',',';')
    globalAttributes['id'] = ID
    globalAttributes['project'] = tDS.attrs['project']
    globalAttributes['comment'] = tDS.attrs['comment']
    globalAttributes['network'] = tDS.attrs['network']
    globalAttributes['geospatial_lat_min'] = tDS.attrs['geospatial_lat_min']
    globalAttributes['geospatial_lat_max'] = tDS.attrs['geospatial_lat_max']
    globalAttributes['geospatial_lat_resolution'] = tDS.attrs['geospatial_lat_resolution']
    globalAttributes['geospatial_lon_min'] = tDS.attrs['geospatial_lon_min']
    globalAttributes['geospatial_lon_max'] = tDS.attrs['geospatial_lon_max']
    globalAttributes['geospatial_lon_resolution'] = tDS.attrs['geospatial_lon_resolution']
    globalAttributes['geospatial_vertical_max'] = tDS.attrs['geospatial_vertical_max']
    globalAttributes['geospatial_vertical_resolution'] = tDS.attrs['geospatial_vertical_resolution']
    globalAttributes['spatial_resolution'] = str(networkData.iloc[0]['grid_resolution'])
    globalAttributes['time_coverage_start'] = timeCoverageStart.strftime('%Y-%m-%dT%H:%M:%SZ')
    globalAttributes['time_coverage_end'] = timeCoverageEnd.strftime('%Y-%m-%dT%H:%M:%SZ')
    globalAttributes['time_coverage_resolution'] = tDS.attrs['time_coverage_resolution']
    globalAttributes['time_coverage_duration'] = timeCoverageDuration
    globalAttributes['area'] = tDS.attrs['area']
    globalAttributes['citation'] += networkData.iloc[0]['citation_statement']
    globalAttributes['processing_level'] = tDS.attrs['processing_level']
    globalAttributes['manufacturer'] = tDS.attrs['manufacturer']
    globalAttributes['sensor_model'] = tDS.attrs['manufacturer']

    creationDate = dt.datetime.utcnow()
    globalAttributes['date_created'] = creationDate.strftime('%Y-%m-%dT%H:%M:%SZ')
    globalAttributes['date_modified'] = creationDate.strftime('%Y-%m-%dT%H:%M:%SZ')
    globalAttributes['history'] = 'Data measured from ' + timeCoverageStart.strftime('%Y-%m-%dT%H:%M:%SZ') + ' to ' \
                                  + timeCoverageEnd.strftime('%Y-%m-%dT%H:%M:%SZ') + '. netCDF file created at ' \
                                  + creationDate.strftime('%Y-%m-%dT%H:%M:%SZ') + ' by the European HFR Node.'

    # Add global attributes to the DataSet
    instacDS.attrs = globalAttributes

    # Encode data types, data packing and _FillValue for the data variables of the DataSet
    for vv in instacDS:
        if vv in dataPacking:
            if 'dtype' in dataPacking[vv]:
                instacDS[vv].encoding['dtype'] = dataPacking[vv]['dtype']
            if 'scale_factor' in dataPacking[vv]:
                instacDS[vv].encoding['scale_factor'] = dataPacking[vv]['scale_factor']
            if 'add_offset' in dataPacking[vv]:
                instacDS[vv].encoding['add_offset'] = dataPacking[vv]['add_offset']
            if 'fill_value' in dataPacking[vv]:
                if not vv in ['SCDR', 'SCDT']:
                    instacDS[vv].encoding['_FillValue'] = netCDF4.default_fillvals[
                        np.dtype(dataPacking[vv]['dtype']).kind + str(np.dtype(dataPacking[vv]['dtype']).itemsize)]
                else:
                    instacDS[vv].encoding['_FillValue'] = b' '

            else:
                instacDS[vv].encoding['_FillValue'] = None

    # Update valid_min and valid_max variable attributes according to data packing for data variables
    for vv in instacDS:
        if 'valid_min' in totVariables[vv]:
            if ('scale_factor' in dataPacking[vv]) and ('add_offset' in dataPacking[vv]):
                instacDS[vv].attrs['valid_min'] = np.float_(((totVariables[vv]['valid_min'] - dataPacking[vv][
                    'add_offset']) / dataPacking[vv]['scale_factor'])).astype(dataPacking[vv]['dtype'])
            else:
                instacDS[vv].attrs['valid_min'] = np.float_(totVariables[vv]['valid_min']).astype(
                    dataPacking[vv]['dtype'])
        if 'valid_max' in totVariables[vv]:
            if ('scale_factor' in dataPacking[vv]) and ('add_offset' in dataPacking[vv]):
                instacDS[vv].attrs['valid_max'] = np.float_(((totVariables[vv]['valid_max'] - dataPacking[vv][
                    'add_offset']) / dataPacking[vv]['scale_factor'])).astype(dataPacking[vv]['dtype'])
            else:
                instacDS[vv].attrs['valid_max'] = np.float_(totVariables[vv]['valid_max']).astype(
                    dataPacking[vv]['dtype'])

    # Encode data types, data packing and _FillValue for the coordinate variables of the DataSet
    for cc in instacDS.coords:
        if cc in dataPacking:
            if 'dtype' in dataPacking[cc]:
                instacDS[cc].encoding['dtype'] = dataPacking[cc]['dtype']
            if 'scale_factor' in dataPacking[cc]:
                instacDS[cc].encoding['scale_factor'] = dataPacking[cc]['scale_factor']
            if 'add_offset' in dataPacking[cc]:
                instacDS[cc].encoding['add_offset'] = dataPacking[cc]['add_offset']
            if 'fill_value' in dataPacking[cc]:
                instacDS[cc].encoding['_FillValue'] = netCDF4.default_fillvals[
                    np.dtype(dataPacking[cc]['dtype']).kind + str(np.dtype(dataPacking[cc]['dtype']).itemsize)]
            else:
                instacDS[cc].encoding['_FillValue'] = None

    # Update valid_min and valid_max variable attributes according to data packing for coordinate variables
    for cc in instacDS.coords:
        if 'valid_min' in totVariables[cc]:
            if ('scale_factor' in dataPacking[cc]) and ('add_offset' in dataPacking[cc]):
                instacDS[cc].attrs['valid_min'] = np.float_(((totVariables[cc]['valid_min'] - dataPacking[cc][
                    'add_offset']) / dataPacking[cc]['scale_factor'])).astype(dataPacking[cc]['dtype'])
            else:
                instacDS[cc].attrs['valid_min'] = np.float_(totVariables[cc]['valid_min']).astype(
                    dataPacking[cc]['dtype'])
        if 'valid_max' in totVariables[cc]:
            if ('scale_factor' in dataPacking[cc]) and ('add_offset' in dataPacking[cc]):
                instacDS[cc].attrs['valid_max'] = np.float_(((totVariables[cc]['valid_max'] - dataPacking[cc][
                    'add_offset']) / dataPacking[cc]['scale_factor'])).astype(dataPacking[cc]['dtype'])
            else:
                instacDS[cc].attrs['valid_max'] = np.float_(totVariables[cc]['valid_max']).astype(
                    dataPacking[cc]['dtype'])

    return instacDS


def buildINSTACtotalFilename(networkID, ts, ext):
    """
    This function builds the filename for total files according to
    the structure used by Copernicus Marine Service In Situ TAC,
    i.e. networkID-stationID_YYYY_MM_DD.

    INPUT:
        networkID: ID of the HFR network
        ts: timestamp as datetime object
        ext: file extension

    OUTPUT:
        totFilename: filename for total file.
    """
    # Get the time related part of the filename
    timeStr = ts.strftime("_%Y%m%d")

    # Build the filename
    totFilename = 'GL_TV_HF_' + networkID + '-Total' + timeStr + ext

    return totFilename


def buildINSTACtotalFolder(basePath, networkID, vers):
    """
    This function builds the folder structure for storing total files according to
    the structure used by Copernicus Marine Service In Situ TAC,
    i.e. networkID-stationID_YYYY_MM_DD.

    INPUT:
        basePath: base path
        networkID: ID of the HFR network
        vers: version of the data model

    OUTPUT:
        totFolder: folder path for storing total files.
    """
    # Strip trailing slash character
    basePath = basePath.rstrip('/')

    # Build the folder path
    totFolder = basePath + '/' + networkID + '/Totals/' + vers + '/'

    return totFolder


def buildEHNtotalFilename(networkID, ts, ext):
    """
    This function builds the filename for total files according to
    the structure used by the European HFR Node, i.e. networkID-Total_YYYY_MM_DD_hhmm.

    INPUT:
        networkID: ID of the HFR network
        ts: timestamp as datetime object
        ext: file extension

    OUTPUT:
        totFilename: filename for total file.
    """
    # Get the time related part of the filename
    timeStr = ts.strftime("_%Y_%m_%d_%H%M")

    # Build the filename
    totFilename = networkID + '-Total' + timeStr + ext

    return totFilename


def buildEHNtotalFolder(basePath, ts, vers):
    """
    This function builds the folder structure for storing total files according to
    the structure used by the European HFR Node, i.e. YYYY/YYYY_MM/YYYY_MM_DD/.

    INPUT:
        basePath: base path
        ts: timestamp as datetime object
        vers: version of the data model

    OUTPUT:
        totFolder: folder path for storing total files.
    """
    # Strip trailing slash character
    basePath = basePath.rstrip('/')

    # Get the time related part of the path
    timeStr = ts.strftime("/%Y/%Y_%m/%Y_%m_%d/")

    # Build the folder path
    totFolder = basePath + '/' + vers + timeStr

    return totFolder


def radBinsInSearchRadius(cell, radial, sR, g):
    """
    This function finds out which radial bins are within the spatthresh of the
    origin grid cell.
    The WGS84 CRS is used for distance calculations.

    INPUT:
        cell: Series containing longitudes and latitudes of the origin grid cells
        radial: Radial object
        sR: search radius in meters
        g: Geod object according to the Total CRS

    OUTPUT:
        radInSr: list of the radial bins falling within the search radius of the
                 origin grid cell.
    """
    # Convert grid cell Series and radial bins DataFrame to numpy arrays
    cell = cell.to_numpy()
    radLon = radial.data['LOND'].to_numpy()
    radLat = radial.data['LATD'].to_numpy()
    # Evaluate distances between origin grid cells and radial bins
    az12, az21, cellToRadDist = g.inv(len(radLon) * [cell[0]], len(radLat) * [cell[1]], radLon, radLat)
    # Figure out which radial bins are within the spatthresh of the origin grid cell
    radInSR = np.where(cellToRadDist < sR)[0].tolist()

    return radInSR


def totalLeastSquare(VelHeadStd):
    """
    This function calculates the u/v components of a total vector from 2 to n
    radial vector components using weighted Least Square method.

    INPUT:
        VelHeadStd: DataFrame containing contributing radial velocities, bearings
                    and standard deviations

    OUTPUT:
        u: U component of the total vector
        v: V component of the total vector
        C: covariance matrix
        Cgdop: covariance matrix assuming uniform unit errors for all radials (i.e. all radial std=1)
    """
    # Convert angles from true convention to math convention
    VelHeadStd['HEAD'] = true2mathAngle(VelHeadStd['HEAD'].to_numpy())

    # Form the design matrix (i.e. the angle matrix)
    A = np.stack((np.array([np.cos(np.deg2rad(VelHeadStd['HEAD'])) / VelHeadStd['STD']]),
                  np.array([np.sin(np.deg2rad(VelHeadStd['HEAD'])) / VelHeadStd['STD']])), axis=-1)[0, :, :]

    # Form the velocity vector
    b = (VelHeadStd['VELO'].to_numpy()) / VelHeadStd['STD']

    # Evaluate the covariance matrix C (variance(U) = C(1,1) and variance(V) = C(2,2))
    A2 = np.matmul(A.T, A)
    if np.linalg.det(A2) > 0:
        C = np.linalg.inv(A2)

        # Calculate the u and v for the total vector
        a = np.matmul(C, np.matmul(A.T, b))
        u = a[0]
        v = a[1]

        # Form the design matrix for GDOP evaluation (i.e. setting all radial std to 1)
        Agdop = np.stack(
            (np.array([np.cos(np.deg2rad(VelHeadStd['HEAD']))]), np.array([np.sin(np.deg2rad(VelHeadStd['HEAD']))])),
            axis=-1)[0, :, :]

        # Evaluate the covariance matrix Cgdop for GDOP evaluation (i.e. setting all radial std to 1)
        Agdop2 = np.matmul(Agdop.T, Agdop)
        if np.linalg.det(Agdop2):
            Cgdop = np.linalg.inv(Agdop2)

            return u, v, C, Cgdop

        else:
            u = np.nan
            v = np.nan
            C = np.nan
            Cgdop = np.nan
            return u, v, C, Cgdop

    else:
        u = np.nan
        v = np.nan
        C = np.nan
        Cgdop = np.nan
        return u, v, C, Cgdop


def makeTotalVector(rBins, rDF):
    """
    This function combines radial contributions to get the total vector for each
    grid cell.
    The weighted Least Square method is used for combination.

    INPUT:
        rBins: Series containing contributing radial indices.
        rDF: DataFrame containing input Radials.

    OUTPUT:
        totalData: Series containing u/v components and related errors of
                   total vector for each grid cell.
    """
    # set minimum number of contributing radial sites
    minContrSites = 2
    # set minimum number of contributing radial vectors
    minContrRads = 3

    # create output total Series
    totalData = pd.Series(np.nan, index=range(9))
    # only consider contributing radial sites
    contrRad = rBins[rBins.str.len() != 0]
    # check if there are at least two contributing radial sites
    if contrRad.size >= minContrSites:
        # loop over contributing radial indices for collecting velocities and angles
        contributions = pd.DataFrame()
        for idx in contrRad.index:
            contrVel = rDF.loc[idx]['Radial'].data.VELO[contrRad[idx]]  # pandas Series
            contrHead = rDF.loc[idx]['Radial'].data.HEAD[contrRad[idx]]  # pandas Series
            if 'ETMP' in rDF.loc[idx]['Radial'].data.columns:
                contrStd = rDF.loc[idx]['Radial'].data.ETMP[contrRad[idx]]  # pandas Series
            elif 'HCSS' in rDF.loc[idx]['Radial'].data.columns:
                contrStd = rDF.loc[idx]['Radial'].data.HCSS[contrRad[idx]].apply(
                    lambda x: math.sqrt(x))  # pandas Series
            else:
                contrStd = np.nan
            contrStd = contrStd.rename("STD")  # pandas Series
            contributions = pd.concat(
                [contributions, pd.concat([contrVel, contrHead, contrStd], axis=1)])  # pandas DataFrame

        # Only keep contributing radials with valid standard deviation values (i.e. different from NaN and 0)
        contributions = contributions[contributions.STD.notnull()]
        contributions = contributions[contributions.STD != 0]

        # check if there are at least three contributing radial vectors
        if len(contributions.index) >= minContrRads:
            # combine radial contributions to get total vector for the current grid cell
            u, v, C, Cgdop = totalLeastSquare(contributions)

            if not math.isnan(u):
                # populate Total Series
                totalData.loc[0] = u  # VELU
                totalData.loc[1] = v  # VELV
                totalData.loc[2] = np.sqrt(u ** 2 + v ** 2)  # VELO
                totalData.loc[3] = (360 + np.arctan2(u, v) * 180 / np.pi) % 360  # HEAD
                totalData.loc[4] = math.sqrt(C[0, 0])  # UQAL
                totalData.loc[5] = math.sqrt(C[1, 1])  # VQAL
                totalData.loc[6] = C[0, 1]  # CQAL
                totalData.loc[7] = math.sqrt(np.abs(Cgdop.trace()))  # GDOP
                totalData.loc[8] = len(contributions.index)  # NRAD

    return totalData


def combineRadials(rDF, gridGS, sRad, gRes, tStp, minContrSites=2):
    """
    This function generataes total vectors from radial measurements using the
    weighted Least Square method for combination.

    INPUT:
        rDF: DataFrame containing input Radials; indices must be the site codes.
        gridGS: GeoPandas GeoSeries containing the longitude/latitude pairs of all
            the points in the grid
        sRad: search radius for combination in meters.
        gRes: grid resoultion in meters
        tStp: timestamp in datetime format (YYYY-MM-DD hh:mm:ss)
        minContrSites: minimum number of contributing radial sites (default to 2)

    OUTPUT:
        Tcomb: Total object generated by the combination
        warn: string containing warnings related to the success of the combination
    """
    # Initialize empty warning string
    warn = ''

    # Create empty total with grid
    Tcomb = Total(grid=gridGS)

    # Check if there are enough contributing radial sites
    if rDF.size >= minContrSites:
        # Fill site_source DataFrame with contributing radials information
        siteNum = 0  # initialization of site number
        for Rindex, Rrow in rDF.iterrows():
            siteNum = siteNum + 1
            rad = Rrow['Radial']
            thisRadial = pd.DataFrame(index=[Rindex],
                                      columns=['#', 'Name', 'Lat', 'Lon', 'Coverage(s)', 'RngStep(km)', 'Pattern',
                                               'AntBearing(NCW)'])
            thisRadial['#'] = siteNum
            thisRadial['Name'] = Rindex
            if rad.is_wera:
                if 'Longitude(dd)OfTheCenterOfTheReceiveArray' in rad.metadata.keys():
                    thisRadial['Lon'] = float(rad.metadata['Longitude(dd)OfTheCenterOfTheReceiveArray'][:-1])
                    if rad.metadata['Longitude(dd)OfTheCenterOfTheReceiveArray'][-1] == 'W':
                        thisRadial['Lon'] = -thisRadial['Lon']
                elif 'Longitude(deg-min-sec)OfTheCenterOfTheReceiveArray' in rad.metadata.keys():
                    thisRadial['Lon'] = dms2dd(list(
                        map(int, rad.metadata['Longitude(deg-min-sec)OfTheCenterOfTheReceiveArray'][:-2].split('-'))))
                    if rad.metadata['Longitude(deg-min-sec)OfTheCenterOfTheReceiveArray'][-1] == 'W':
                        thisRadial['Lon'] = -thisRadial['Lon']
                if 'Latitude(dd)OfTheCenterOfTheReceiveArray' in rad.metadata.keys():
                    thisRadial['Lat'] = float(rad.metadata['Latitude(dd)OfTheCenterOfTheReceiveArray'][:-1])
                    if rad.metadata['Latitude(dd)OfTheCenterOfTheReceiveArray'][-1] == 'S':
                        thisRadial['Lat'] = -thisRadial['Lat']
                elif 'Latitude(deg-min-sec)OfTheCenterOfTheReceiveArray' in rad.metadata.keys():
                    thisRadial['Lat'] = dms2dd(list(
                        map(int, rad.metadata['Latitude(deg-min-sec)OfTheCenterOfTheReceiveArray'][:-2].split('-'))))
                    if rad.metadata['Latitude(deg-min-sec)OfTheCenterOfTheReceiveArray'][-1] == 'S':
                        thisRadial['Lat'] = -thisRadial['Lat']
                thisRadial['Coverage(s)'] = float(rad.metadata['ChirpRate'].replace('S', '')) * int(
                    rad.metadata['Samples'])
                thisRadial['RngStep(km)'] = float(rad.metadata['Range'].split()[0])
                thisRadial['Pattern'] = 'Internal'
                thisRadial['AntBearing(NCW)'] = float(rad.metadata['TrueNorth'].split()[0])
            else:
                thisRadial['Lat'] = float(rad.metadata['Origin'].split()[0])
                thisRadial['Lon'] = float(rad.metadata['Origin'].split()[1])
                thisRadial['Coverage(s)'] = float(rad.metadata['TimeCoverage'].split()[0])
                if 'RangeResolutionKMeters' in rad.metadata:
                    thisRadial['RngStep(km)'] = float(rad.metadata['RangeResolutionKMeters'].split()[0])
                elif 'RangeResolutionMeters' in rad.metadata:
                    thisRadial['RngStep(km)'] = float(rad.metadata['RangeResolutionKMeters'].split()[0]) * 0.001
                thisRadial['Pattern'] = rad.metadata['PatternType'].split()[0]
                thisRadial['AntBearing(NCW)'] = float(rad.metadata['AntennaBearing'].split()[0])
            Tcomb.site_source = pd.concat([Tcomb.site_source, thisRadial])

        # Insert timestamp
        Tcomb.time = tStp

        # Fill Total with some metadata
        Tcomb.metadata['TimeZone'] = rad.metadata[
            'TimeZone']  # trust all radials have the same, pick from the last radial
        Tcomb.metadata['AveragingRadius'] = str(sRad / 1000) + ' km'
        Tcomb.metadata['GridAxisOrientation'] = '0.0 DegNCW'
        Tcomb.metadata['GridSpacing'] = str(gRes / 1000) + ' km'

        # Create Geod object according to the Total CRS
        g = Geod(ellps=Tcomb.metadata['GreatCircle'].split()[0])

        # Create DataFrame for storing indices of radial bins falling within the search radius of each grid cell
        combineRadBins = pd.DataFrame(columns=range(len(Tcomb.data.index)))

        # Figure out which radial bins are within the spatthresh of each grid cell
        for Rindex, Rrow in rDF.iterrows():
            rad = Rrow['Radial']
            thisRadBins = Tcomb.data.loc[:, ['LOND', 'LATD']].apply(lambda x: radBinsInSearchRadius(x, rad, sRad, g),
                                                                    axis=1)
            combineRadBins.loc[Rindex] = thisRadBins

        # Loop over grid points and pull out contributing radial vectors
        combineRadBins = combineRadBins.T
        totData = combineRadBins.apply(lambda x: makeTotalVector(x, rDF), axis=1)

        # Assign column names to the combination DataFrame
        totData.columns = ['VELU', 'VELV', 'VELO', 'HEAD', 'UQAL', 'VQAL', 'CQAL', 'GDOP', 'NRAD']

        # Fill Total with combination results
        Tcomb.data[['VELU', 'VELV', 'VELO', 'HEAD', 'UQAL', 'VQAL', 'CQAL', 'GDOP', 'NRAD']] = totData

        # Mask out vectors on land
        Tcomb.mask_over_land(subset=True)

        # Get the indexes of grid cells without total vectors
        indexNoVec = Tcomb.data[Tcomb.data['VELU'].isna()].index
        # Delete these row indexes from DataFrame
        Tcomb.data.drop(indexNoVec, inplace=True)
        Tcomb.data.reset_index(level=None, drop=False,
                               inplace=True)  # Set drop=True if the former indices are not necessary

        if Tcomb.data.empty:
            warn = 'No combination performed: no overlap in radial coverages'

    else:
        warn = 'No combination performed: not enough contributing radial sites'

    return Tcomb, warn


class Total(fileParser):
    """
    Totals Subclass.

    This class should be used when loading CODAR (.tuv) and WERA (.cur_asc) total files.
    This class utilizes the generic LLUV and CUR classes.
    """

    def __init__(self, fname='', replace_invalid=True, grid=gpd.GeoSeries(), empty_total=False):

        if not fname:
            empty_total = True
            replace_invalid = False

        super().__init__(fname)
        for key in self._tables.keys():
            table = self._tables[key]
            if 'LLUV' in table['TableType']:
                self.data = table['data']
            elif 'src' in table['TableType']:
                self.diagnostics_source = table['data']
            elif 'CUR' in table['TableType']:
                self.cur_data = table['data']

        if 'SiteSource' in self.metadata.keys():
            if not self.is_wera:
                table_data = u''
                for ss in self.metadata['SiteSource']:
                    if '%%' in ss:
                        rep = {' comp': '_comp', ' Distance': '_Distance', ' Ratio': '_Ratio', ' (dB)': '_(dB)',
                               ' Width': '_Width', ' Resp': '_Resp', 'Value ': 'Value_', 'FOL ': 'FOL_'}
                        rep = dict((re.escape(k), v) for k, v in rep.items())
                        pattern = re.compile('|'.join(rep.keys()))
                        ss_header = pattern.sub(lambda m: rep[re.escape(m.group(0))], ss).strip('%% SiteSource \n')
                    else:
                        ss = ss.replace('%SiteSource:', '').strip()
                        ss = ss.replace('Radial', '').strip()
                        table_data += '{}\n'.format(ss)
                # use pandas read_csv because it interprets the datatype for each column of the csv
                tdf = pd.read_csv(
                    io.StringIO(table_data),
                    sep=' ',
                    header=None,
                    names=ss_header.split(),
                    skipinitialspace=True
                )
                self.site_source = tdf

        # Evaluate GDOP for total files
        if hasattr(self, 'site_source'):
            # Get the site coordinates
            siteLon = self.site_source['Lon'].values.tolist()
            siteLat = self.site_source['Lat'].values.tolist()
            # Create Geod object according to the Total CRS, if defined. Otherwise use WGS84 ellipsoid
            if self.metadata['GreatCircle']:
                g = Geod(ellps=self.metadata['GreatCircle'].split()[0].replace('"', ''))
            else:
                g = Geod(ellps='WGS84')
                self.metadata['GreatCircle'] = '"WGS84"' + ' ' + str(g.a) + '  ' + str(1 / g.f)
            self.data['GDOP'] = self.data.loc[:, ['LOND', 'LATD']].apply(lambda x: evaluateGDOP(x, siteLon, siteLat, g),
                                                                         axis=1)
        elif hasattr(self, 'data'):
            self.data['GDOP'] = np.nan

        # Evaluate the number of contributing radials (NRAD) for CODAR total files
        if hasattr(self, 'is_wera'):
            if not self.is_wera:
                if hasattr(self, 'data'):
                    self.data['NRAD'] = self.data.loc[:, self.data.columns.str.contains('S.*CN')].sum(axis=1)

        if replace_invalid:
            self.replace_invalid_values()

        if empty_total:
            self.empty_total()

        # if mask_over_land:
        #     self.mask_over_land()

        if not grid.empty:
            self.initialize_grid(grid)

    def empty_total(self):
        """
        Create an empty Total object. The empty Total object can be created by setting
        the geographical grid.
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
            elif 'rads' in table['TableType']:
                self.diagnostics_radial.drop(self.diagnostics_radial.index[:], inplace=True)
                self._tables[key]['data'] = self.diagnostics_radial
            elif 'rcvr' in table['TableType']:
                self.diagnostics_hardware.drop(self.diagnostics_hardware.index[:], inplace=True)
                self._tables[key]['data'] = self.diagnostics_hardware
            elif 'RINF' in table['TableType']:
                self.range_information.drop(self.range_information.index[:], inplace=True)
                self._tables[key]['data'] = self.range_information
            elif 'CUR' in table['TableType']:
                self.cur_data.drop(self.cur_data.index[:], inplace=True)
                self._tables[key]['data'] = self.cur_data

        if not hasattr(self, 'data'):
            self.data = pd.DataFrame()

        if hasattr(self, 'site_source'):
            self.site_source.drop(self.site_source.index[:], inplace=True)
        else:
            self.site_source = pd.DataFrame()

    def initialize_grid(self, gridGS):
        """
        Initialize the geogprahic grid for filling the LOND and LATD columns of the
        Total object data DataFrame.

        INPUT:
            gridGS: GeoPandas GeoSeries containing the longitude/latitude pairs of all
                the points in the grid

        OUTPUT:
            DataFrame with filled LOND and LATD columns.
        """

        # initialize data DataFrame with column names
        self.data = pd.DataFrame(
            columns=['LOND', 'LATD', 'VELU', 'VELV', 'VELO', 'HEAD', 'UQAL', 'VQAL', 'CQAL', 'GDOP', 'NRAD'])

        # extract longitudes and latitude from grid GeoSeries and insert them into data DataFrame
        self.data['LOND'] = gridGS.x
        self.data['LATD'] = gridGS.y

        # add metadata about datum and CRS
        self.metadata = OrderedDict()
        self.metadata['GreatCircle'] = ''.join(gridGS.crs.ellipsoid.name.split()) + ' ' + str(
            gridGS.crs.ellipsoid.semi_major_metre) + '  ' + str(gridGS.crs.ellipsoid.inverse_flattening)

    def mask_over_land(self, subset=False, res='high'):
        """
        This function masks the total vectors lying on land.
        Total vector coordinates are checked against a reference file containing information
        about which locations are over land or in an unmeasurable area (for example, behind an
        island or point of land).
        The Natural Earth public domain maps are used as reference.
        If "res" option is set to "high", the map with 10 m resolution is used, otherwise the map with 110 m resolution is used.
        The EPSG:4326 CRS is used for distance calculations.
        If "subset" option is set to True, the total vectors lying on land are removed.

        INPUT:
            subset: option enabling the removal of total vectors on land (if set to True)
            res: resolution of the www.naturalearthdata.com dataset used to perform the masking; None or 'low' or 'high'. Defaults to 'high'.

        OUTPUT:
            waterIndex: list containing the indices of total vectors lying on water.
        """
        # Load the reference file (GeoPandas "naturalearth_lowres")
        mask_dir = '.hfradarpy'
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
        This function plots the current total velocity field (i.e. VELU and VELV components) on a
        Cartesian grid. The grid is defined either from the input values or from the Total object
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
        siteLon = self.site_source['Lon'].to_numpy()
        siteLat = self.site_source['Lat'].to_numpy()
        siteCode = self.site_source['Name'].tolist()

        # Compute the native map projection coordinates for the stations
        xS, yS = m(siteLon, siteLat)

        # Plot radial stations
        m.plot(xS, yS, 'rD')
        for label, xs, ys in zip(siteCode, xS, yS):
            plt.text(xs, ys, label, fontdict={'fontsize': 22, 'fontweight': 'bold'})

        # Plot velocity field
        if shade:
            self.to_xarray_multidimensional()

            # Create grid from longitude and latitude
            [longitudes, latitudes] = np.meshgrid(self.xdr['LONGITUDE'].data, self.xdr['LATITUDE'].data)

            # Compute the native map projection coordinates for the pseudo-color cells
            X, Y = m(longitudes, latitudes)

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
        plt.title(self.file_name + ' total velocity field', fontdict={'fontsize': 30, 'fontweight': 'bold'})

        if show:
            plt.show()

        return fig

    def plot(self, lon_min=None, lon_max=None, lat_min=None, lat_max=None, shade=False, show=True):
        """
        This function plots the current total velocity field (i.e. VELU and VELV components) on a
        Cartesian grid. The grid is defined either from the input values or from the Total object
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

                # Initialize the figure
        cmap = 'jet'  # set the colorbar
        norm = colors.Normalize(vmin=0, vmax=1)  # set colorbar limits
        extent = [lon_min, lon_max, lat_min, lat_max]  # set the map extent
        fig = plt.figure(num=None, figsize=(24, 24), dpi=100, facecolor='w', edgecolor='k')
        ax = plt.axes(projection=ccrs.Mercator())  # set the map projection
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5,
                          linestyle='--')  # add grid lines
        ax.add_feature(cfeature.LAND)  # add land
        ax.add_feature(cfeature.OCEAN)  # add ocean
        ax.add_feature(cfeature.COASTLINE)  # add coastline
        ax.set_extent(extent)  # apply the map extent

        # Plot radial stations
        for Rindex, site in self.site_source.iterrows():
            plt.plot(site['Lon'], site['Lat'], color='r', markeredgecolor='k', marker='o', markersize=10,
                     transform=ccrs.Geodetic())
            ax.text(site['Lon'], site['Lat'], site['Name'], transform=ccrs.Geodetic(),
                    fontdict={'fontsize': 22, 'fontweight': 'bold'})

        # Scale the velocity component variables
        if self.is_wera:
            u = self.data.VELU
            v = self.data.VELV
        else:
            u = self.data.VELU / 100  # CODAR velocities are in cm/s
            v = self.data.VELV / 100  # CODAR velocities are in cm/s

        # Plot velocity field
        if shade:
            self.to_xarray_multidimensional()

            # Create grid from longitude and latitude
            [X, Y] = np.meshgrid(self.xdr['LONGITUDE'].data, self.xdr['LATITUDE'].data)

            # Create velocity variable in the shape of the grid
            M = abs(self.xdr['VELO'][0, 0, :, :].to_numpy())
            # V = V[:-1,:-1]

            # Make the pseudo-color plot
            warnings.simplefilter("ignore", category=UserWarning)
            ax.pcolormesh(X, Y, M, transform=ccrs.PlateCarree(), shading='nearest', cmap=cmap, vmin=0, vmax=1)

            # Get the longitude, latitude and velocity components
            x, y, U, V = self.data['LOND'].values, self.data['LATD'].values, u.values, v.values

            # Evaluate the velocity magnitude
            m = (U ** 2 + V ** 2) ** 0.5

            # Normalize velocity components
            Un, Vn = U / m, V / m

            # Make the quiver plot
            Q = ax.quiver(x, y, Un, Vn, transform=cartopy.crs.PlateCarree(), cmap=cmap, norm=norm, scale=20)

            warnings.simplefilter("default", category=UserWarning)

        else:
            # Get the longitude, latitude and velocity components
            x, y, U, V = self.data['LOND'].values, self.data['LATD'].values, u.values, v.values

            # Evaluate the velocity magnitude
            m = (U ** 2 + V ** 2) ** 0.5

            # Make the quiver plot
            Q = ax.quiver(x, y, U, V, m, transform=cartopy.crs.PlateCarree(), cmap=cmap, norm=norm, scale=5)
            # Add the reference arrow
            ax.quiverkey(Q, 0.1, 0.9, 0.5, r'$0.5 m/s$', fontproperties={'size': 12, 'weight': 'bold'})

        # Add colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        plt.colorbar(sm, ax=ax, orientation='vertical', pad=0.1).ax.set_xlabel('velocity (m/s)', labelpad=10, )

        # Add title
        plt.title(self.file_name + ' total velocity field', fontdict={'fontsize': 30, 'fontweight': 'bold'}, pad=25)

        if show:
            plt.show()

        return fig

    def to_xarray_multidimensional(self, lon_min=None, lon_max=None, lat_min=None, lat_max=None, grid_res=None):
        """
        This function creates a dictionary of xarray DataArrays containing the variables
        of the Total object multidimensionally expanded along the coordinate axes (T,Z,Y,X).
        The coordinate axes are set as (TIME, DEPTH, LATITUDE, LONGITUDE).
        The coordinate limits and steps are taken from total metadata when possible.
        Only longitude and latitude limits and grid resolution can be specified by the user.
        Some refinements are performed on Codar data in order to comply CF convention
        for positive velocities and to have velocities in m/s.
        The generated dictionary is attached to the Total object, named as xdr.

        INPUT:
            lon_min: minimum longitude value in decimal degrees (if None it is taken from Total metadata)
            lon_max: maximum longitude value in decimal degrees (if None it is taken from Total metadata)
            lat_min: minimum latitude value in decimal degrees (if None it is taken from Total metadata)
            lat_max: maximum latitude value in decimal degrees (if None it is taken from Total metadata)
            grid_res: grid resolution in meters (if None it is taken from Total metadata)

        OUTPUT:

        """
        # Initialize empty dictionary
        xdr = OrderedDict()

        # process Codar totals
        if not self.is_wera:
            # Get longitude limits
            if lon_min is None:
                if 'BBminLongitude' in self.metadata:
                    lon_min = float(self.metadata['BBminLongitude'].split()[0])
                else:
                    lon_min = self.data.LOND.min()
            if lon_max is None:
                if 'BBmaxLongitude' in self.metadata:
                    lon_max = float(self.metadata['BBmaxLongitude'].split()[0])
                else:
                    lon_max = self.data.LOND.max()

            # Get latitude limits
            if lat_min is None:
                if 'BBminLatitude' in self.metadata:
                    lat_min = float(self.metadata['BBminLatitude'].split()[0])
                else:
                    lat_min = self.data.LATD.min()
            if lat_max is None:
                if 'BBmaxLatitude' in self.metadata:
                    lat_max = float(self.metadata['BBmaxLatitude'].split()[0])
                else:
                    lat_max = self.data.LATD.max()

                    # Get grid resolution
            if grid_res is None:
                if 'GridSpacing' in self.metadata:
                    grid_res = float(self.metadata['GridSpacing'].split()[0]) * 1000
                else:
                    grid_res = float(1)

                    # Generate grid coordinates
            gridGS = createLonLatGridFromBB(lon_min, lon_max, lat_min, lat_max, grid_res)

        # process WERA totals
        else:
            if not self.is_combined:
                # Get longitude limits and step
                if 'TopLeftLongitude' in self.metadata:
                    top_left_lon = float(self.metadata['TopLeftLongitude'].split()[0])
                else:
                    top_left_lon = self.data.LOND.min()
                if 'NX' in self.metadata:
                    cells_lon = int(self.metadata['NX'].split()[0])
                else:
                    cells_lon = 100

                # Get latitude limits and step
                if 'TopLeftLatitude' in self.metadata:
                    top_left_lat = float(self.metadata['TopLeftLatitude'].split()[0])
                else:
                    top_left_lat = self.data.LATD.max()
                if 'NY' in self.metadata:
                    cells_lat = int(self.metadata['NY'].split()[0])
                else:
                    cells_lat = 100

                # Get cell size in km
                if 'DGT' in self.metadata:
                    cell_size = float(self.metadata['DGT'].split()[0])
                else:
                    cell_size = float(1)

                # Generate grid coordinates
                gridGS = createLonLatGridFromTopLeftPointWera(top_left_lon, top_left_lat, cell_size, cells_lon,
                                                              cells_lat)

            else:
                # Get longitude limits
                if lon_min is None:
                    if 'BBminLongitude' in self.metadata:
                        lon_min = float(self.metadata['BBminLongitude'].split()[0])
                    else:
                        lon_min = self.data.LOND.min()
                if lon_max is None:
                    if 'BBmaxLongitude' in self.metadata:
                        lon_max = float(self.metadata['BBmaxLongitude'].split()[0])
                    else:
                        lon_max = self.data.LOND.max()

                # Get latitude limits
                if lat_min is None:
                    if 'BBminLatitude' in self.metadata:
                        lat_min = float(self.metadata['BBminLatitude'].split()[0])
                    else:
                        lat_min = self.data.LATD.min()
                if lat_max is None:
                    if 'BBmaxLatitude' in self.metadata:
                        lat_max = float(self.metadata['BBmaxLatitude'].split()[0])
                    else:
                        lat_max = self.data.LATD.max()

                        # Get grid resolution
                if grid_res is None:
                    if 'GridSpacing' in self.metadata:
                        grid_res = float(self.metadata['GridSpacing'].split()[0]) * 1000
                    else:
                        grid_res = float(1)

                        # Generate grid coordinates
                gridGS = createLonLatGridFromBBwera(lon_min, lon_max, lat_min, lat_max, grid_res)

        # extract longitudes and latitude from grid GeoSeries and insert them into numpy arrays
        lon_dim = np.unique(gridGS.x.to_numpy())
        lat_dim = np.unique(gridGS.y.to_numpy())
        # manage antimeridian crossing
        lon_dim = np.concatenate((lon_dim[lon_dim >= 0], lon_dim[lon_dim < 0]))

        # # Get the longitude and latitude values of the total measurements
        # unqLon = np.sort(np.unique(self.data['LOND']))
        # unqLat = np.sort(np.unique(self.data['LATD']))

        # # Insert unqLon and unqLat values to replace the closest in lon_dim and lat_dim
        # replaceIndLon = abs(unqLon[None, :] - lon_dim[:, None]).argmin(axis=0).tolist()
        # replaceIndLat = abs(unqLat[None, :] - lat_dim[:, None]).argmin(axis=0).tolist()
        # lon_dim[replaceIndLon] = unqLon
        # lat_dim[replaceIndLat] = unqLat

        # Create total grid from longitude and latitude
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

        # create dictionary containing variables from dataframe in the shape of total grid
        d = {key: np.tile(np.nan, longitudes.shape) for key in self.data.keys()}

        # Remap all variables
        for k, v in d.items():
            v[Y_map_idx.astype(int), X_map_idx.astype(int)] = self.data[k]
            d[k] = v

        # Add extra dimensions for time (T) and depth (Z) - CF Standard: T, Z, Y, X -> T=axis0, Z=axis1
        d = {k: np.expand_dims(np.float32(v), axis=(0, 1)) for (k, v) in d.items()}

        # Drop LOND and LATD variables (they are set as coordinates of the DataSet)
        d.pop('LOND')
        d.pop('LATD')

        # Refine Codar and combined data
        if not self.is_wera:
            # Scale velocities to be in m/s (only for Codar or combined totals with is_wera attribute set to False)
            toMs = ['VELU', 'VELV', 'VELO', 'UQAL', 'VQAL', 'CQAL']
            for t in toMs:
                if t in d:
                    d[t] = d[t] * 0.01

        # Evaluate timestamp as number of days since 1950-01-01T00:00:00Z
        timeDelta = self.time - dt.datetime.strptime('1950-01-01T00:00:00Z', '%Y-%m-%dT%H:%M:%SZ')
        ncTime = timeDelta.days + timeDelta.seconds / (60 * 60 * 24)

        # Add all variables as xarray
        for k, v in d.items():
            xdr[k] = xr.DataArray(v,
                                  dims={'TIME': v.shape[0], 'DEPTH': v.shape[1], 'LATITUDE': v.shape[2],
                                        'LONGITUDE': v.shape[3]},
                                  coords={'TIME': [ncTime],
                                          'DEPTH': [0],
                                          'LATITUDE': lat_dim,
                                          'LONGITUDE': lon_dim})

        # Add DataArray for coordinate variables
        xdr['TIME'] = xr.DataArray(ncTime,
                                   dims={'TIME': len(pd.date_range(self.time, periods=1))},
                                   coords={'TIME': [ncTime]})
        xdr['DEPTH'] = xr.DataArray(0,
                                    dims={'DEPTH': 1},
                                    coords={'DEPTH': [0]})
        xdr['LATITUDE'] = xr.DataArray(lat_dim,
                                       dims={'LATITUDE': lat_dim},
                                       coords={'LATITUDE': lat_dim})
        xdr['LONGITUDE'] = xr.DataArray(lon_dim,
                                        dims={'LONGITUDE': lon_dim},
                                        coords={'LONGITUDE': lon_dim})

        # Attach the dictionary to the Total object
        self.xdr = xdr

        return

    def check_ehn_mandatory_variables(self):
        """
        This function checks if the Total object contains all the mandatory data variables
        (i.e. not coordinate variables) required by the European standard data model developed in the framework of the
        EuroGOOS HFR Task Team.
        Missing variables are appended to the DataFrame containing data, filled with NaNs.

        INPUT:

        OUTPUT:
        """
        # Set mandatory variables based on the HFR manufacturer
        if self.is_wera:
            chkVars = ['VELU', 'VELV', 'VELO', 'UACC', 'VACC', 'GDOP']
        else:
            chkVars = ['VELU', 'VELV', 'UQAL', 'VQAL', 'CQAL', 'GDOP']

        # Check variables and add missing ones
        for vv in chkVars:
            if vv not in self.data.columns:
                self.data[vv] = np.nan

        return

    def apply_ehn_datamodel(self, network_data, station_data, version):
        """
        This function applies the European standard data model developed in the
        framework of the EuroGOOS HFR Task Team to the Total object.
        The Total object content is stored into an xarray Dataset built from the
        xarray DataArrays created by the Total method to_xarray_multidimensional.
        Variable data types and data packing information are collected from
        "Data_Models/EHN/Totals/Total_Data_Packing.json" file.
        Variable attribute schema is collected from
        "Data_Models/EHN/Totals/Total_Variables.json" file.
        Global attribute schema is collected from
        "Data_Models/EHN/Global_Attributes.json" file.
        Global attributes are created starting from Total object metadata and from
        DataFrames containing the information about HFR network and radial stations
        read from the EU HFR NODE database.
        The generated xarray Dataset is attached to the Total object, named as xds.

        INPUT:
            network_data: DataFrame containing the information of the network to which the radial site belongs
            station_data: DataFrame containing the information of the radial site that produced the radial
            version: version of the data model


        OUTPUT:
        """
        # Set the netCDF format
        ncFormat = 'NETCDF4_CLASSIC'

        # Expand Total object variables along the coordinate axes
        if self.is_combined:
            self.to_xarray_multidimensional()
        else:
            # Get bounding box limits and grid resolution from database
            lonMin = network_data.iloc[0]['geospatial_lon_min']
            lonMax = network_data.iloc[0]['geospatial_lon_max']
            latMin = network_data.iloc[0]['geospatial_lat_min']
            latMax = network_data.iloc[0]['geospatial_lat_max']
            gridRes = network_data.iloc[0]['grid_resolution'] * 1000
            self.to_xarray_multidimensional(lonMin, lonMax, latMin, latMax, gridRes)

        # Set auxiliary coordinate sizes
        maxsiteSize = 150
        refmaxSize = 50
        maxinstSize = 50

        # Get data packing information per variable
        f = open('Data_Models/EHN/Totals/Total_Data_Packing.json')
        dataPacking = json.loads(f.read())
        f.close()

        # Get variable attributes
        f = open('Data_Models/EHN/Totals/Total_Variables.json')
        totVariables = json.loads(f.read())
        f.close()

        # Get global attributes
        f = open('Data_Models/EHN/Global_Attributes.json')
        globalAttributes = json.loads(f.read())
        f.close()

        # Rename velocity related and quality related variables
        self.xdr['EWCT'] = self.xdr.pop('VELU')
        self.xdr['NSCT'] = self.xdr.pop('VELV')
        if 'UQAL' in self.xdr:
            self.xdr['EWCS'] = self.xdr.pop('UQAL')
        if 'VQAL' in self.xdr:
            self.xdr['NSCS'] = self.xdr.pop('VQAL')
        if 'CQAL' in self.xdr:
            self.xdr['CCOV'] = self.xdr.pop('CQAL')

        # Drop unnecessary DataArrays from the DataSet
        toDrop = ['VFLG', 'XDST', 'YDST', 'RNGE', 'BEAR', 'NRAD', 'VELO', 'HEAD', 'index']
        for t in toDrop:
            if t in self.xdr:
                self.xdr.pop(t)
        toDrop = list(self.xdr.keys())
        for t in toDrop:
            if fnmatch.fnmatch(t, 'S*CN'):
                self.xdr.pop(t)
        toDrop = []
        for vv in self.xdr:
            if vv not in totVariables.keys():
                toDrop.append(vv)
        for rv in toDrop:
            self.xdr.pop(rv)

            # Add coordinate reference system to the dictionary
        self.xdr['crs'] = xr.DataArray(int(0), )

        # Add antenna related variables to the dictionary
        # Number of antennas
        contributingSiteNrx = station_data.loc[station_data['station_id'].isin(self.site_source.Name.tolist())][
            'number_of_receive_antennas'].to_numpy()
        nRX = np.asfarray(contributingSiteNrx)
        nRX = np.pad(nRX, (0, maxsiteSize - len(nRX)), 'constant', constant_values=(np.nan, np.nan))
        contributingSiteNtx = station_data.loc[station_data['station_id'].isin(self.site_source.Name.tolist())][
            'number_of_transmit_antennas'].to_numpy()
        nTX = np.asfarray(contributingSiteNtx)
        nTX = np.pad(nTX, (0, maxsiteSize - len(nTX)), 'constant', constant_values=(np.nan, np.nan))
        self.xdr['NARX'] = xr.DataArray([nRX],
                                        dims={'TIME': len(pd.date_range(self.time, periods=1)), 'MAXSITE': maxsiteSize})
        self.xdr['NATX'] = xr.DataArray([nTX],
                                        dims={'TIME': len(pd.date_range(self.time, periods=1)), 'MAXSITE': maxsiteSize})

        # Longitude and latitude of antennas
        contributingSiteLat = station_data.loc[station_data['station_id'].isin(self.site_source.Name.tolist())][
            'site_lat'].to_numpy()
        siteLat = np.pad(contributingSiteLat, (0, maxsiteSize - len(contributingSiteLat)), 'constant',
                         constant_values=(np.nan, np.nan))
        contributingSiteLon = station_data.loc[station_data['station_id'].isin(self.site_source.Name.tolist())][
            'site_lon'].to_numpy()
        siteLon = np.pad(contributingSiteLon, (0, maxsiteSize - len(contributingSiteLon)), 'constant',
                         constant_values=(np.nan, np.nan))
        self.xdr['SLTR'] = xr.DataArray([siteLat],
                                        dims={'TIME': len(pd.date_range(self.time, periods=1)), 'MAXSITE': maxsiteSize})
        self.xdr['SLNR'] = xr.DataArray([siteLon],
                                        dims={'TIME': len(pd.date_range(self.time, periods=1)), 'MAXSITE': maxsiteSize})
        self.xdr['SLTT'] = xr.DataArray([siteLat],
                                        dims={'TIME': len(pd.date_range(self.time, periods=1)), 'MAXSITE': maxsiteSize})
        self.xdr['SLNT'] = xr.DataArray([siteLon],
                                        dims={'TIME': len(pd.date_range(self.time, periods=1)), 'MAXSITE': maxsiteSize})

        # Codes of antennas
        contributingSiteCodeList = station_data.loc[station_data['station_id'].isin(self.site_source.Name.tolist())][
            'station_id'].tolist()
        antCode = np.array([site.encode() for site in contributingSiteCodeList])
        antCode = np.pad(antCode, (0, maxsiteSize - len(contributingSiteCodeList)), 'constant',
                         constant_values=('', ''))
        self.xdr['SCDR'] = xr.DataArray(np.array([antCode]),
                                        dims={'TIME': len(pd.date_range(self.time, periods=1)), 'MAXSITE': maxsiteSize})
        self.xdr['SCDR'].encoding['char_dim_name'] = 'STRING' + str(len(station_data['station_id'].to_numpy()[0]))
        self.xdr['SCDT'] = xr.DataArray(np.array([antCode]),
                                        dims={'TIME': len(pd.date_range(self.time, periods=1)), 'MAXSITE': maxsiteSize})
        self.xdr['SCDT'].encoding['char_dim_name'] = 'STRING' + str(len(station_data['station_id'].to_numpy()[0]))

        # Add SDN namespace variables to the dictionary
        siteCode = ('%s' % network_data.iloc[0]['network_id']).encode()
        self.xdr['SDN_CRUISE'] = xr.DataArray([siteCode], dims={'TIME': len(pd.date_range(self.time, periods=1))})
        self.xdr['SDN_CRUISE'].encoding['char_dim_name'] = 'STRING' + str(len(siteCode))
        platformCode = ('%s' % network_data.iloc[0]['network_id'] + '-Total').encode()
        self.xdr['SDN_STATION'] = xr.DataArray([platformCode], dims={'TIME': len(pd.date_range(self.time, periods=1))})
        self.xdr['SDN_STATION'].encoding['char_dim_name'] = 'STRING' + str(len(platformCode))
        ID = ('%s' % platformCode.decode() + '_' + self.time.strftime('%Y-%m-%dT%H:%M:%SZ')).encode()
        self.xdr['SDN_LOCAL_CDI_ID'] = xr.DataArray([ID], dims={'TIME': len(pd.date_range(self.time, periods=1))})
        self.xdr['SDN_LOCAL_CDI_ID'].encoding['char_dim_name'] = 'STRING' + str(len(ID))
        sdnEDMO = np.asfarray(pd.concat([network_data['EDMO_code'], station_data['EDMO_code']]).unique())
        sdnEDMO = np.pad(sdnEDMO, (0, maxinstSize - len(sdnEDMO)), 'constant', constant_values=(np.nan, np.nan))
        self.xdr['SDN_EDMO_CODE'] = xr.DataArray([sdnEDMO], dims={'TIME': len(pd.date_range(self.time, periods=1)),
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
            self.xds[vv].attrs = totVariables[vv]

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
            self.xds[cc].attrs = totVariables[cc]

        # Evaluate measurement maximum depth
        vertMax = 3e8 / (8 * np.pi * station_data['transmit_central_frequency'].to_numpy().min() * 1e6)

        # Evaluate time coverage start, end, resolution and duration
        timeCoverageStart = self.time - relativedelta(minutes=network_data.iloc[0]['temporal_resolution'] / 2)
        timeCoverageEnd = self.time + relativedelta(minutes=network_data.iloc[0]['temporal_resolution'] / 2)
        timeResRD = relativedelta(minutes=network_data.iloc[0]['temporal_resolution'])
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
        globalAttributes.pop('oceanops_ref')
        globalAttributes.pop('wmo_platform_code')
        globalAttributes.pop('wigos_id')
        globalAttributes['doa_estimation_method'] = ', '.join(
            station_data[["station_id", "DoA_estimation_method"]].apply(": ".join, axis=1))
        globalAttributes['calibration_type'] = ', '.join(
            station_data[["station_id", "calibration_type"]].apply(": ".join, axis=1))
        if 'HFR-US' in network_data.iloc[0]['network_id']:
            station_data['last_calibration_date'] = 'N/A'
            globalAttributes['last_calibration_date'] = ', '.join(
                pd.concat([station_data['station_id'], station_data['last_calibration_date']], axis=1)[
                    ["station_id", "last_calibration_date"]].apply(": ".join, axis=1))
        else:
            globalAttributes['last_calibration_date'] = ', '.join(pd.concat([station_data['station_id'], station_data[
                'last_calibration_date'].apply(lambda x: x.strftime('%Y-%m-%dT%H:%M:%SZ'))], axis=1)[
                                                                      ["station_id", "last_calibration_date"]].apply(
                ": ".join, axis=1))
            globalAttributes['last_calibration_date'] = globalAttributes['last_calibration_date'].replace(
                '1-01-01T00:00:00Z', 'N/A')
        globalAttributes['calibration_link'] = ', '.join(
            station_data[["station_id", "calibration_link"]].apply(": ".join, axis=1))
        # globalAttributes['title'] = network_data.iloc[0]['title']
        globalAttributes['title'] = 'Near Real Time Surface Ocean Total Velocity by ' + globalAttributes[
            'platform_code']
        globalAttributes['summary'] = network_data.iloc[0]['summary']
        globalAttributes['institution'] = ', '.join(
            pd.concat([network_data['institution_name'], station_data['institution_name']]).unique().tolist())
        globalAttributes['institution_edmo_code'] = ', '.join(
            [str(x) for x in pd.concat([network_data['EDMO_code'], station_data['EDMO_code']]).unique().tolist()])
        globalAttributes['institution_references'] = ', '.join(
            pd.concat([network_data['institution_website'], station_data['institution_website']]).unique().tolist())
        globalAttributes['id'] = ID.decode()
        globalAttributes['project'] = network_data.iloc[0]['project']
        globalAttributes['comment'] = network_data.iloc[0]['comment']
        globalAttributes['network'] = network_data.iloc[0]['network_name']
        globalAttributes['data_type'] = globalAttributes['data_type'].replace('current data', 'total current data')
        globalAttributes['geospatial_lat_min'] = str(network_data.iloc[0]['geospatial_lat_min'])
        globalAttributes['geospatial_lat_max'] = str(network_data.iloc[0]['geospatial_lat_max'])
        globalAttributes['geospatial_lat_resolution'] = str(network_data.iloc[0]['grid_resolution'])
        globalAttributes['geospatial_lon_min'] = str(network_data.iloc[0]['geospatial_lon_min'])
        globalAttributes['geospatial_lon_max'] = str(network_data.iloc[0]['geospatial_lon_max'])
        globalAttributes['geospatial_lon_resolution'] = str(network_data.iloc[0]['grid_resolution'])
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
        globalAttributes['processing_level'] = '3B'
        globalAttributes['contributor_name'] = network_data.iloc[0]['contributor_name']
        globalAttributes['contributor_role'] = network_data.iloc[0]['contributor_role']
        globalAttributes['contributor_email'] = network_data.iloc[0]['contributor_email']
        globalAttributes['manufacturer'] = ', '.join(
            station_data[["station_id", "manufacturer"]].apply(": ".join, axis=1))
        globalAttributes['sensor_model'] = ', '.join(
            station_data[["station_id", "manufacturer"]].apply(": ".join, axis=1))
        globalAttributes['software_version'] = version

        creationDate = dt.datetime.utcnow()
        globalAttributes['metadata_date_stamp'] = creationDate.strftime('%Y-%m-%dT%H:%M:%SZ')
        globalAttributes['date_created'] = creationDate.strftime('%Y-%m-%dT%H:%M:%SZ')
        globalAttributes['date_modified'] = creationDate.strftime('%Y-%m-%dT%H:%M:%SZ')
        globalAttributes['history'] = 'Data collected at ' + self.time.strftime(
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
            if 'valid_min' in totVariables[vv]:
                if ('scale_factor' in dataPacking[vv]) and ('add_offset' in dataPacking[vv]):
                    self.xds[vv].attrs['valid_min'] = np.float_(((totVariables[vv]['valid_min'] - dataPacking[vv][
                        'add_offset']) / dataPacking[vv]['scale_factor'])).astype(dataPacking[vv]['dtype'])
                else:
                    self.xds[vv].attrs['valid_min'] = np.float_(totVariables[vv]['valid_min']).astype(
                        dataPacking[vv]['dtype'])
            if 'valid_max' in totVariables[vv]:
                if ('scale_factor' in dataPacking[vv]) and ('add_offset' in dataPacking[vv]):
                    self.xds[vv].attrs['valid_max'] = np.float_(((totVariables[vv]['valid_max'] - dataPacking[vv][
                        'add_offset']) / dataPacking[vv]['scale_factor'])).astype(dataPacking[vv]['dtype'])
                else:
                    self.xds[vv].attrs['valid_max'] = np.float_(totVariables[vv]['valid_max']).astype(
                        dataPacking[vv]['dtype'])

        # Encode data types and avoid data packing, valid_min, valid_max and _FillValue for the coordinate variables of the DataSet
        for cc in self.xds.coords:
            if cc in dataPacking:
                if 'dtype' in dataPacking[cc]:
                    self.xds[cc].encoding['dtype'] = dataPacking[cc]['dtype']
                if 'valid_min' in totVariables[cc]:
                    del self.xds[cc].attrs['valid_min']
                if 'valid_max' in totVariables[cc]:
                    del self.xds[cc].attrs['valid_max']
                self.xds[cc].encoding['_FillValue'] = None

        return

    def initialize_qc(self):
        """
        Initialize dictionary entry for QC metadata.
        """
        # Initialize dictionary entry for QC metadta
        self.metadata['QCTest'] = {}

    def qc_ehn_maximum_velocity(self, totMaxSpeed=1.2):
        """
        This test labels total velocity vectors whose module is smaller than a maximum velocity threshold
        with a “good data” flag. Otherwise the vectors are labeled with a “bad data” flag.
        The ARGO QC flagging scale is used.

        This test was defined in the framework of the EuroGOOS HFR Task Team based on the
        Max Speed Threshold test (QC303) from the Integrated Ocean Observing System (IOOS) Quality Assurance of
        Real-Time Oceanographic Data (QARTOD).

        INPUTS:
            totMaxSpeed: maximum velocity in m/s for normal operations
        """
        # Set the test name
        testName = 'CSPD_QC'

        # Add new column to the DataFrame for QC data by setting every row as passing the test (flag = 1)
        self.data.loc[:, testName] = 1

        # set bad flag for velocities not passing the test
        if self.is_wera:
            self.data.loc[(self.data['VELO'].abs() > totMaxSpeed), testName] = 4  # velocity in m/s (CRAD)
        else:
            self.data.loc[(self.data['VELO'].abs() > totMaxSpeed * 100), testName] = 4  # velocity in cm/s (LLUV)

        self.metadata['QCTest'][testName] = 'Velocity Threshold QC Test - Test applies to each vector. ' \
                                            + 'Threshold=[' + f'maximum velocity={totMaxSpeed} (m/s)]'

    def qc_ehn_maximum_variance(self, totMaxVar=1):
        """
        This test labels total velocity vectors whose temporal variances for both U and V
        components are smaller than a maximum variance threshold with a “good data” flag.
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
            totMaxVar: maximum variance in m2/s2 for normal operations
        """
        # Set the test name
        testName = 'VART_QC'

        # Add new column to the DataFrame for QC data by setting every row as passing the test (flag = 1)
        self.data.loc[:, testName] = 1

        # Set bad flag for variances not passing the test
        if self.is_wera and not self.is_combined:
            self.data.loc[(self.data[
                               'UACC'] ** 2 > totMaxVar), testName] = 4  # UACC is the temporal standard deviation of U component in m/s for WERA non-combined total data
            self.data.loc[(self.data[
                               'VACC'] ** 2 > totMaxVar), testName] = 4  # VACC is the temporal standard deviation of V component in m/s for WERA non-combined total data
        else:
            self.data.loc[((self.data[
                                'UQAL'] / 100) ** 2 > totMaxVar), testName] = 4  # UQAL is the temporal standard deviation of U component in cm/s for CODAR and combined total data
            self.data.loc[((self.data[
                                'VQAL'] / 100) ** 2 > totMaxVar), testName] = 4  # VQAL is the temporal standard deviation of V component in cm/s for CODAR and combined total data

        self.metadata['QCTest'][testName] = 'Variance Threshold QC Test - Test applies to each vector. ' \
                                            + 'Threshold=[' + f'maximum variance={totMaxVar} (m2/s2)]'

    def qc_ehn_gdop_threshold(self, maxGDOP=2):
        """
        This test labels total velocity vectors whose GDOP is smaller than a maximum GDOP threshold
        with a “good data” flag. Otherwise the vectors are labeled with a “bad data” flag.
        The ARGO QC flagging scale is used.

        This test was defined in the framework of the EuroGOOS HFR Task Team based on the
        GDOP Threshold test (QC302) from the Integrated Ocean Observing System (IOOS) Quality Assurance of
        Real-Time Oceanographic Data (QARTOD).

        INPUTS:
            maxGDOP: maximum allowed GDOP for normal operations
        """
        # Set the test name
        testName = 'GDOP_QC'

        # Add new column to the DataFrame for QC data by setting every row as passing the test (flag = 1)
        self.data.loc[:, testName] = 1

        # set bad flag for velocities not passing the test
        self.data.loc[(self.data['GDOP'] > maxGDOP), testName] = 4

        self.metadata['QCTest'][testName] = 'GDOP Threshold QC Test - Test applies to each vector. ' \
                                            + 'Threshold=[' + f'GDOP threshold={maxGDOP}]'

    def qc_ehn_data_density_threshold(self, minContrRad=2):
        """
        This test labels total velocity vectors with a number of contributing radial velocities bigger
        than the minimum number defined for normal operations with a “good data” flag.
        Otherwise the vectors are labeled with a “bad data” flag.
        The ARGO QC flagging scale is used.

        This test was defined in the framework of the EuroGOOS HFR Task Team based on the
        Data Density Threshold test (QC301) from the Integrated Ocean Observing System (IOOS) Quality Assurance of
        Real-Time Oceanographic Data (QARTOD).

        INPUTS:
            minContrRad: minimum number of contributing radial velocities for normal operations
        """
        # Set the test name
        testName = 'DDNS_QC'

        # Add new column to the DataFrame for QC data by setting every row as passing the test (flag = 1)
        self.data.loc[:, testName] = 1

        # set bad flag for velocities not passing the test
        if not self.is_wera:
            if 'NRAD' in self.data.columns:
                self.data.loc[(self.data['NRAD'] < minContrRad), testName] = 4
            else:
                self.data.loc[:, testName] = 0

        self.metadata['QCTest'][testName] = 'Data Density Threshold QC Test - Test applies to each vector. ' \
                                            + 'Threshold=[' + f'minimum number of contributing radial velocities={minContrRad}]'

    def qc_ehn_temporal_derivative(self, t0, tempDerThr=1):
        """
        This test compares the velocity of each total vector with the velocity of the total vector
        measured in the previous timestamp at the same location.
        Each vector for which the velocity difference is smaller than the specified threshold for normal
        operations (tempDerThr), is labeled with a "good data" flag.
        Otherwise the vector is labeled with a “bad data” flag.
        The ARGO QC flagging scale is used.

        This test was defined in the framework of the EuroGOOS HFR Task Team based on the
        Temporal Gradient test (QC206) from the Integrated Ocean Observing System (IOOS) Quality
        Assurance of Real-Time Oceanographic Data (QARTOD).

        INPUTS:
            t0: Total object of the previous timestamp
            tempDerThr: velocity difference threshold in m/s for normal operations
        """
        # Set the test name
        testName = 'VART_QC'

        # Check if the previous timestamp total file exists
        if not t0 is None:
            # Merge the data DataFrame of the two Totals and evaluate velocity differences at each location
            mergedDF = self.data.merge(t0.data, on=['LOND', 'LATD'], how='left', suffixes=(None, '_x'),
                                       indicator='Exist')
            velDiff = (mergedDF['VELO'] - mergedDF['VELO_x']).abs()

            # Add new column to the DataFrame for QC data by setting every row as passing the test (flag = 1)
            self.data.loc[:, testName] = 1

            # Set rows of the DataFrame for QC data as not evaluated (flag = 0) for locations existing in the current total but not in the previous one
            self.data.loc[mergedDF['Exist'] == 'left_only', testName] = 0

            # Set bad flag for vectors not passing the test
            if self.is_wera:
                self.data.loc[(velDiff > tempDerThr), testName] = 4  # velocity in m/s (CUR)
            else:
                self.data.loc[(velDiff > tempDerThr * 100), testName] = 4  # velocity in cm/s (LLUV)

        else:
            # Add new column to the DataFrame for QC data by setting every row as not evaluated (flag = 0)
            self.data.loc[:, testName] = 0

        self.metadata['QCTest'][testName] = 'Temporal Derivative QC Test - Test applies to each vector. ' \
                                            + 'Threshold=[' + f'velocity difference threshold={str(tempDerThr)} (m/s)]'

    def qc_ehn_overall_qc_flag(self):
        """

        This QC test labels total velocity vectors with a ‘good_data” flag if all QC tests are passed.
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

    def file_type(self):
        """
        Return a string representing the type of file this is.
        """
        return 'totals'