import numpy as np
import math
from pyproj import Geod
from shapely.geometry import Point
from geopandas import GeoSeries
import pandas as pd
import datetime as dt


def roundToNearest20(dT):
    """
    This function rounds teh input datetime to the closest 20-minute time. 
    
    INPUT:
        dT: datetime as datetime object or string
        
        
    OUTPUT:
        dT: datetime rounded to the closest 20-minutes time as datetime object
    """
    # Get the minutes from the timestamp
    if not hasattr(dT,'minute'):
        # Convert dT to datetime if needeed
        dT = dt.datetime.strptime(dT,'%Y-%m-%d %H:%M:%S')
    minute = dT.minute
        
    # Round the datetime
    if minute <= 10:
        return dT.replace(minute=0, second=0, microsecond=0)
    elif minute <= 30:
        return dT.replace(minute=20, second=0, microsecond=0)
    elif minute <= 50:
        return dT.replace(minute=40, second=0, microsecond=0)
    else:
        return (dT + pd.Timedelta(minutes=20)).replace(minute=0, second=0, microsecond=0)


def createLonLatGridFromTopLeftPointWera(topLeftLon, topLeftLat, cellSize, nx, ny):
    """
    This function creates a regular logitude/latitude grid given the longitude 
    and latitude of the top-left point, the grid cell size in km and the number of
    grid cells along x and y axes.
    The grid is created following the method developed by Helzel for WERA systems.
    The WGS84 CRS is used.
    
    INPUT:
        topLeftLon: longitude of the top-left point of the grid
        topLeftLat: latitude of the top-left point of the grid
        cellSize: spatial resolution of the grid, expressed in kilometers
        nx: number of grid cells along x axis
        ny: number of grid cells along y axis
        
        
    OUTPUT:
        pts: GeoPandas GeoSeries containing the longitude/latitude pairs of all
            the points in the grid
    """
    # evaluate kilometers per degree (derived from Nautical Mile)
    kmDegree = 1.852 * 60

    # evaluate distance of one lat-grid cell in Degrees
    dLat = cellSize / kmDegree

    # generate latitude list
    latVec = np.arange(topLeftLat,topLeftLat-(ny*dLat),-dLat)
    if len(latVec) > ny:
        latVec = latVec[:ny]

    # evaluate shperical correction
    sphCorr = math.cos(math.radians(topLeftLat - (0.5 *(ny-1)*dLat)))

    # evaluate distance of one lon-grid cell in Degrees
    dLon = cellSize / (kmDegree * sphCorr)

    # generate longitude list
    lonVec = np.arange(topLeftLon,topLeftLon+(nx*dLon),dLon)
    if len(lonVec) > nx:
        lonVec = lonVec[:nx]
        
    # manage antimeridian crossing
    lonVec[lonVec>180] = lonVec[lonVec>180]-360

    # Create grid
    Lon, Lat = np.meshgrid(lonVec, latVec);
    Lonc = Lon.flatten()
    Latc = Lat.flatten()

    # Now convert these points to geo-data
    pts = GeoSeries([Point(x, y) for x, y in zip(Lonc, Latc)])
    pts = pts.set_crs('epsg:4326')
    
    return pts


def createLonLatGridFromBBwera(lonMin, lonMax, latMin, latMax, gridResolution):
    """
    This function creates a regular logitude/latitude grid given the longitude 
    and latitude limits for the bounding box and the spatial resolution in km.
    The grid is created following the method developed by Helzel for WERA systems.
    The WGS84 CRS is used.
    
    INPUT:
        lonMin: westernmost longitude of the bounding box
        lonMax: easternmost longitude of the bounding box
        latMin: southernmost longitude of the bounding box
        latMax: northernmost longitude of the bounding box
        gridResolution: spatial resolution of the grid, expressed in meters
        
    OUTPUT:
        pts: GeoPandas GeoSeries containing the longitude/latitude pairs of all
            the points in the grid
    """
    # manage antimeridian crossing
    antiMeridianCrossing = False
    if lonMin>0 and lonMax<0:
        antiMeridianCrossing = True
        lonMax = lonMax + 360
    
    # evaluate kilometers per degree (derived from Nautical Mile)
    kmDegree = 1.852 * 60

    # evaluate distance of one lat-grid cell in Degrees
    dLat = (gridResolution/1000) / kmDegree

    # generate latitude list
    latVec = np.arange(latMax,latMin,-dLat)

    # evaluate number of latitude cells
    latCells = len(latVec)

    # evaluate shperical correction
    sphCorr = math.cos(math.radians(latMax - (0.5 *(latCells-1)*dLat)))

    # evaluate distance of one lon-grid cell in Degrees
    dLon = (gridResolution/1000) / (kmDegree * sphCorr)

    # generate longitude list
    lonVec = np.arange(lonMin,lonMax,dLon)
    
    # manage antimeridian crossing
    if antiMeridianCrossing:
        lonVec[lonVec>180] = lonVec[lonVec>180]-360

    # Create grid
    Lon, Lat = np.meshgrid(lonVec, latVec);
    Lonc = Lon.flatten()
    Latc = Lat.flatten()

    # Now convert these points to geo-data
    pts = GeoSeries([Point(x, y) for x, y in zip(Lonc, Latc)])
    pts = pts.set_crs('epsg:4326')
    
    return pts


def createLonLatGridFromBB(lonMin, lonMax, latMin, latMax, gridResolution):
    """
    This function creates a regular logitude/latitude grid given the longitude 
    and latitude limits for the bounding box and the spatial resolution in km.
    The WGS84 CRS is used.
    
    INPUT:
        lonMin: westernmost longitude of the bounding box
        lonMax: easternmost longitude of the bounding box
        latMin: southernmost longitude of the bounding box
        latMax: northernmost longitude of the bounding box
        gridResolution: spatial resolution of the grid, expressed in meters
        
    OUTPUT:
        pts: GeoPandas GeoSeries containing the longitude/latitude pairs of all
            the points in the grid
    """
    # manage antimeridian crossing
    antiMeridianCrossing = False
    if lonMin>0 and lonMax<0:
        antiMeridianCrossing = True
        lonMax = lonMax + 360
        
    # Use WGS84 ellipsoid
    g = Geod(ellps='WGS84')

    # Compute forward and back azimuths, plus distance between lower left and upper left corners
    az12,az21,dist_ULtoLL = g.inv(lonMin,latMin,lonMin,latMax)

    # Retrieve the array of distances from lower left points along the latitude axis based on the grid resolution
    dd = np.arange(0,dist_ULtoLL,gridResolution)

    # Compute latitude, longitude and back azimuth of all the points along the latitude axis
    fooLon, Lat, backaz = g.fwd(len(dd)*[lonMin], len(dd)*[latMin], len(dd)*[0], dd)
    Lat = np.array(Lat)

    # Retrieve coordinates of the center grid point
    lonCenter = (lonMin + lonMax) / 2
    latCenter = (latMin + latMax) / 2

    # Evaluate meters distance of small displacement scaled up to 1 degree longitude displacement
    az12,az21,dist_smallDisp1degLon = g.inv(lonCenter,latCenter,lonCenter+1e-4,latCenter)
    dist_smallDisp1degLon = 1e4 * dist_smallDisp1degLon

    # Evaluate the appropriate displacement in longitude
    dd = gridResolution / dist_smallDisp1degLon

    # Compute longitudes of all the points along the longitude axis
    Lon = np.arange(lonMin,lonMax,dd)
    
    # Adjust Lon and Lat so that they are centered in box
    Lon = Lon + (lonMax - Lon[-1]) / 2;
    Lat = Lat + (latMax - Lat[-1]) / 2;
    
    # manage antimeridian crossing
    if antiMeridianCrossing:
        Lon[Lon>180] = Lon[Lon>180]-360

    # Create grid
    Lon, Lat = np.meshgrid(Lon, Lat);
    Lonc = Lon.flatten()
    Latc = Lat.flatten()

    # Now convert these points to geo-data
    pts = GeoSeries([Point(x, y) for x, y in zip(Lonc, Latc)])
    pts = pts.set_crs('epsg:4326')
    
    return pts


def true2mathAngle(trueAngle, radians=False):
    """
    This function converts angles from the true system (i.e. the geographical
    system) to the math system (i.e. the trigonometry system).
    
    In the math convention an angle is measured CCW from East (i.e. East = 0 degrees,
    North = 90 degrees, etc.).
    In the true convention an angle is measured CW from North (i.e. North = 0 degrees,
    East = 90 degrees, etc.). This is also commonly referred to as compass angle.
    
    INPUT:
        trueAngle: numpy array containing true angles
        radians: flag indicating whether the angles are in degrees (default) or in radians
        
    OUTPUT:
        mathAngle: numpy array containing math angles
    """
    # Convert to degrees if angles are in radians
    if radians:
        trueAngle = np.rad2deg(trueAngle)
        
    mathAngle = 90 - trueAngle
    mathAngle = np.mod(mathAngle,360)
    
    return mathAngle


def dms2dd(dms):
    """
    This function converts angles from DMS (degrees-minutes-seconds) to DD 
    (decimal degrees).
    
    INPUT:
        dms: tuple containing degree, minute and second values
        
    OUTPUT:
        dd: decimal degrees angle
    """
    dd = float(dms[0]) + float(dms[1])/60 + float(dms[2])/(60*60)
    
    return dd


def evaluateGDOP(cell, siteLon, siteLat, g):
    """
    This function evaluates the GDOP value of a grid cell based on its coordinates 
    and on the coordinates of the radial sites.
    
    INPUT:
        cell: Series containing longitude and latitude of the grid cell for which the GDOP is evaluated
        siteLon: list containing the longitudes of the radial sites
        siteLat: list containing the latitudes of the radial sites
        g: Geod object with CRS.
        
    OUTPUT:
        gdop: GDOP value
    """
    # Convert grid cell Series to numpy arrays
    cell = cell.to_numpy()
    cellLon = cell[0]
    cellLat = cell[1]
    
    # Evaluate the radial angles from the radial sites to the grid cell
    radialAngles,az21,dist = g.inv(siteLon,siteLat,len(siteLon)*[cellLon],len(siteLat)*[cellLat])
    
    # Form the design matrix for GDOP evaluation
    Agdop = np.stack((np.array([np.cos(np.deg2rad(radialAngles))]),np.array([np.sin(np.deg2rad(radialAngles))])),axis=-1)[0,:,:]
    
    # Evaluate the covariance matrix Cgdop for GDOP evaluation
    Cgdop = np.linalg.inv(np.matmul(Agdop.T, Agdop))
    
    # Evaluate GDOP
    gdop = math.sqrt(Cgdop.trace())
    
    return gdop




