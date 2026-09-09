import logging
import numpy as np
import math
from hfradarpy.calc import scircle1, inpolygon, lonlat2km

class BaseObject:
    """the base object to provice some common functions"""

    def to_dict(self):
        """return the dict for all attributes in the object"""
        return dict(vars(self))

    def from_dict(self, **kwargs):
        """set the object attributes"""
        var = vars(self)
        for k, v in kwargs.items():
            if k in var:
                if isinstance(getattr(self, k), BaseObject):
                    getattr(self, k).from_dict(**v)
                else:
                    setattr(self, k, v)


class Grid(BaseObject):
    """
        The 'grid' class keeps info about the ocean velocities
    """

    def __init__(self):
        self.resolution_km = np.float64(0)
        self.projection = ""
        self.x_range = np.array([], dtype=np.float64)
        self.y_range = np.array([], dtype=np.float64)
        self.dx = np.float64(0)
        self.dy = np.float64(0)
        self.size = np.array([], dtype=np.float64)
        self.ocean_indices = np.array([], dtype=np.float64)
        self.ocean_xy = [[]]
        self.ocean_x_scircle = [[]]
        self.ocean_y_scircle = [[]]


class LandInfo:
    """
        Land information
    """

    def __init__(self, floatArray1, floatArray2, intVar):
        self.region = np.array(floatArray1)
        self.polygon = np.array(floatArray2)
        self.length = intVar


class RtvTotals(BaseObject):
    """
        Information needed for the radial totals
        All are arrays of variable length and grid of class Grid
    """

    def __init__(self):
        self.lat = np.array([], dtype=np.float64)
        self.lon = np.array([], dtype=np.float64)
        self.u_xvel = np.array([], dtype=np.float64)
        self.v_yvel = np.array([], dtype=np.float64)
        self.dopx = np.array([], dtype=np.float64)
        self.dopy = np.array([], dtype=np.float64)
        self.hdop = np.array([], dtype=np.float64)
        self.nRads = np.array([], dtype=np.float64)
        self.nRadSites = np.array([], dtype=np.float64)
        self.nSites = np.array([], dtype=np.float64)
        self.grid = Grid()
        self.land = LandInfo(self.lat, self.lon, 10)

    def to_dict(self):
        """return the dict for all attributes in the object"""
        d = super().to_dict()
        d['lat'] = np.array(d['lat'])
        d['lon'] = np.array(d['lon'])
        d['grid'] = d['grid'].to_dict()
        d.pop('land', None)
        return d


class SiteInfo:
    """
       Info about the data collection site
    """

    def __init__(self, network="", name="", beampattern="", useMinute=0):
        self.network = network
        self.name = name
        self.beampattern = beampattern
        self.useMinute = useMinute


class RtvInfo(BaseObject):
    """
       All the different variable parameters for RTV running
    """

    def __init__(self):
        self.grid_search_radius = np.float64(0)
        self.max_age = np.float64(0)
        self.max_rad_speed = np.float64(0)
        self.max_rtv_speed = np.float64(0)
        self.min_rad_sites = np.float64(0)
        self.min_radials = np.float64(0)
        self.uwls_max_hdop = np.float64(0)
        self.uwls_max_hdop_ascii = np.float64(0)
        self.uwls_max_hdop_nc = np.float64(0)
        self.oi_mdlvar = np.float64(0)
        self.oi_errvar = np.float64(0)
        self.oi_sx = np.float64(0)
        self.oi_sy = np.float64(0)
        self.oi_weighting = ''
        self.method = None
        self.new_state = None
        self.current_state = None


class UwlsTotalInfo:
    """
       Vector calculations info
    """

    def __init__(self):
        self.u = np.float64(np.nan)
        self.v = np.float64(np.nan)
        self.dopx = np.float64(np.nan)
        self.dopy = np.float64(np.nan)
        self.hdop = np.float64(np.nan)


class OITotalInfo:
    """
       Vector calculations info
    """

    def __init__(self):
        self.u = np.float64(np.nan)
        self.v = np.float64(np.nan)
        self.dopx = np.float64(np.nan)
        self.dopy = np.float64(np.nan)
        self.hdop = np.float64(np.nan)
        self.covr = np.float64(np.nan)



class RadialInfo(BaseObject):
    """
       Radial velocity info needed for calculations
    """

    def __init__(self):
        self.time = ""
        self.site = ""
        self.network = ""
        self.patterntype = ""
        self.manufacturer = ""
        self.file = ""
        self.dir = ""
        self.sitelatitude = np.float64(0)
        self.sitelongitude = 9999.0
        self.maxrange = 9999.0
        self.isNew = False
        self.longitude = []
        self.latitude = []
        self.speed = []
        self.heading = []


class compTotInfo(BaseObject):
    """
        Information needed for the radial totals
        All are arrays of variable length and grid of class Grid
    """

    def __init__(self):
        self.grid = Grid()


def uwlsTotal(rSpeed, rHeading) -> UwlsTotalInfo:
    """
    input of two vectors:
    rSpeed   - Column vector (n* x 1) of radial velocity magnitude
    rHeading - Column vector (n* x 1) of radial velocity heading
               in degres counterclockwise from +x (east)
    n* >= 2
    """

    retValUwlsTotal = UwlsTotalInfo()

    # copying the calculation
    X = np.zeros((np.size(rHeading), 2))
    X[:, 0] = np.cos(np.deg2rad(rHeading))
    X[:, 1] = np.sin(np.deg2rad(rHeading))

    X_transpose_X = np.matmul(np.transpose(X), X)
    det_transpose = np.linalg.det(X_transpose_X)

    # can't do inverse
    if abs(det_transpose) < 1e-5:
        return retValUwlsTotal

    C = np.linalg.inv(X_transpose_X)

    retValUwlsTotal.dopx = np.sqrt(C.item(0, 0))
    retValUwlsTotal.dopy = np.sqrt(C.item(1, 1))
    retValUwlsTotal.hdop = np.sqrt(C.item(0, 0) + C.item(1, 1))

    b = np.linalg.multi_dot([C, np.transpose(X), rSpeed])

    retValUwlsTotal.u = b[0]
    retValUwlsTotal.v = b[1]

    return retValUwlsTotal

def oiTotal(rSpeed,rHeading,rpLon,rpLat,gridloc,mdlvar,errvar,sx,sy,oi_weighting ) -> OITotalInfo:
    """
    input:
    rSpeed   - Column vector (n* x 1) of radial velocity magnitude
    rHeading - Column vector (n* x 1) of radial velocity heading
               in degres counterclockwise from +x (east)
    n* >= 2

        rloc: longitude and latitude associated with the radial location
		gridloc: vector grid location where uv vector components are estimated
        mdlvar: a priori model covariance of the surface currents
                     In the pointwise approach, user can set the a priori
                      model variance as a function of water depth, the
                      length from the coastline, or constant.
        errvar: observational error variance, which can be constant
						or hourly standard deviation of radial velocity (HSTD, so called
						temporal uncertainty). For example, if the observational
						uncertainty in HF radar observations is 3-5 cm/s, the corresponding
						error variance = 9-25 (cm/s)^2.
						As a standard error,
							HSTD^2/N can be used as the error variance.
						N is the number of cross spectra within a given time period.
						Although N is unknown, the bound of N is known. 1<= N <= 6.
		sx,sy: decorrelation length scales
		oi_weighting: option for correlation function; option == 1, Gaussian, option == 2, exponential function

    OUTPUT:
        u: U component of the total vector
        v: V component of the total vector
        xi: uncertainty normalized by the a priori model covariance. (2 by 2 matrix)
			xi(0,0) : normalized uncertainty of u = <(u_hat - u)^2>/<u^2> (good :0, poor: 1)
			xi(1,1) : normalized uncertainty of v = <(v_hat - v)^2>/<v^2> (good :0, poor: 1)
			xi(0,1) : directional information of u and v = <(u_hat -u)(v_hat- v)>/sqrt(<u^2><v^2>)

    REFERENCE: Kim, S.Y., Terrill, E.J. and Cornuelle, B.D., 2008. Mapping surface currents from
    HF radar radial velocity measurements using optimal interpolation.
    Journal of Geophysical Research: Oceans, 113(C10).
    """

    retValOITotal = OITotalInfo()

    # Form the design matrix (i.e. the angle matrix)
    X = np.zeros((np.size(rHeading), 2))
    X[:, 0] = np.cos(np.deg2rad(rHeading))
    X[:, 1] = np.sin(np.deg2rad(rHeading))

    X_transpose_X = np.matmul(np.transpose(X), X)
    det_transpose = np.linalg.det(X_transpose_X)

    # can't do inverse
    if abs(det_transpose) < 1e-5:
        return retValUwlsTotal

    nr = np.shape(rHeading)[0]

    # Matrix R based on the length of errvar
    if isinstance(errvar, int) or isinstance(errvar, float):
        R = np.eye(nr) * errvar   # For constant error variance
    else:
        R = np.diag(errvar)  # For HSTD (hourly temporal uncertainty)

    # P_ matrix (identity matrix scaled by mdlvar)
    P_ = np.eye(2) * mdlvar

    # lonlat2km equivalent
    #dx, dy = lonlat2km(rloc[:,0], rloc[:,1], gridloc[:, 0], gridloc[:, 1])
    dx, dy = lonlat2km(rpLon, rpLat, gridloc[:, 0], gridloc[:, 1])

    # Meshgrid operations for x and y locations
    #x1, x2 = np.meshgrid(rloc[:,0], rloc[:,0])
    #y1, y2 = np.meshgrid(rloc[:,1], rloc[:,1])
    x1, x2 = np.meshgrid(rpLon, rpLon)
    y1, y2 = np.meshgrid(rpLat, rpLat)
    ang1, ang2 = np.meshgrid(rHeading, rHeading)

    # lonlat2km for meshgrids
    drx_, dry_ = lonlat2km(x1, y1, x2, y2)

    if oi_weighting == 'gaussian':
        cmd = np.exp(-(dx ** 2 / sx ** 2 + dy ** 2 / sy ** 2))
        w_ = np.exp(-(drx_ ** 2 / sx ** 2 + dry_ ** 2 / sy ** 2)) * mdlvar
    elif oi_weighting == 'exponential':
        cmd = np.exp(-np.sqrt(dx ** 2 / sx ** 2 + dy ** 2 / sy ** 2))
        w_ = np.exp(-np.sqrt(drx_ ** 2 / sx ** 2 + dry_ ** 2 / sy ** 2)) * mdlvar
    else:
        warn()

    # Compute cmd in components (u and v)
    cmdu = cmd * np.cos(rHeading)
    cmdv = cmd * np.sin(rHeading)
    cmd = np.column_stack([cmdu, cmdv])

    # Calculate cdd
    cdd = w_ * (np.cos(ang1) * np.cos(ang2) + np.sin(ang1) * np.sin(ang2))

    # Compute cmdicdd and final values
    #cmdicdd = np.dot(P_, np.dot(cmd.T, np.linalg.inv(cdd + R))) #slower method
    A = cdd + R
    X = np.linalg.solve(A, cmd)
    cmdicdd = P_ @ X.T
    a = np.dot(cmdicdd, rSpeed)
    xi = np.linalg.inv(P_) @ (P_ - np.dot(np.dot(cmdicdd, cmd), P_))

    retValOITotal.dopx = xi[0, 0] # Uerr, normalized uncertainty of u (good: 0 poor: 1)
    retValOITotal.dopy = xi[1, 1]  # Verr, normalized uncertainty of v (good: 0 poor: 1)
    retValOITotal.hdop = math.sqrt(xi[0,0] ** 2 + xi[1,1] ** 2)  # OI Total Errors
    retValOITotal.covr = xi[0, 1]  # directional info of u and v (UV covariance)

    retValOITotal.u = a[0]
    retValOITotal.v = a[1]

    return retValOITotal



def rtvComputeTotals(compTotInput, radials, rtvInfoObj) -> RtvTotals:
    """
    Returns a structure containing total velocity solutions comptued from
    radial velocity measurements using input processing parameters.

    Inputs:
    Configuration structure, radial data structure

    Output:
    U - total solution structure with fields:
        lat, lon, u, v, dopx, dopy, hdop, nRads, nSites,
        grid: {resolution_km, projection, x_range, y_range, dx, dy, size,
               ocean_indices}

    Specific configuration fields required by rtvComputetotals are:
    grid, domain, resolution, rtv('min_rad_sites'), rtv('grid_search_radius'),
    rtv('min_radials'), rtv('max_rtv_speed'), rtv('uwls_max_hdop_mat')

    Small circles are used to narrow in on grid points were solutions are
    possible based on each site's origin and maximum radial data range and
    where the number of overlapping sites meet or exceed the
    rtv('min_rad_sites') value.

    Pre-computed small circles for each grid point corresponding to the
    rtv('grid_search_radius') value are then used to find radials that
    fall within the circle.  A total solution for a given gridpoint is
    computed if the radial data (1) contain data from a 'new' site, (2)
    have equal or greater than rtv('min_rad_sites') sites contributing,
    and (3) have equal to more than rtv('min_radials') radials.

    Finally, solutions are filtered for complex and infinate values as
    well as for speed exceeding rtv('max_rtv_speed') HDOP value exceeding
    rtv('uwls_max_hdop_mat').

    """

    U_totals = RtvTotals()

    # reduce total grid solution space
    nArrayLen = len(compTotInput.grid.ocean_indices)
    gridAllRadCount = np.zeros(nArrayLen)
    nRads = np.zeros(nArrayLen)
    nRadSites = np.zeros(nArrayLen)
    gridNewRadIndex = np.zeros(nArrayLen, dtype=bool)
    nSites = len(radials)

    # loop over each radial dataset (site)
    for iRadial in range(nSites):

        currRadial = radials[iRadial]

        # Compute small circlebased on maximum range of data, adding grid search radius, using WGS85
        scLat, scLon = scircle1(currRadial.sitelatitude,
                                currRadial.sitelongitude,
                                (currRadial.maxrange +
                                 rtvInfoObj.grid_search_radius))

        # Find total grid points inside the small circle
        grid_ocean_x = compTotInput.grid.ocean_xy[0]
        grid_ocean_y = compTotInput.grid.ocean_xy[1]
        inPoints = np.zeros(len(grid_ocean_x))

        if len(grid_ocean_x) != len(grid_ocean_y):
            logging.error("ERROR: Don't have same amount of x and y vals")
            return U_totals

        inPoints = inpolygon(grid_ocean_x, grid_ocean_y, scLon, scLat)

        # Increment grid count for points inside radial coverage
        gridAllRadCount = gridAllRadCount + inPoints

        # Index grid points covered by new data
        if currRadial.isNew:
            gridNewRadIndex = np.logical_or(gridNewRadIndex, inPoints)

    # Define grid points with potential solutions based on new radial coverage
    # and number of sites overlapping the grid point
    sPoint = np.argwhere(np.logical_and(gridNewRadIndex,
                                        (gridAllRadCount >= rtvInfoObj.min_rad_sites)))

    if len(sPoint) == 0:
        logging.info('No potential total solution points found')
        return U_totals

    # Define grid small circle field name
    scircle_xfield = ""
    scircle_yfield = ""
    if rtvInfoObj.grid_search_radius == np.floor(rtvInfoObj.grid_search_radius):
        scircle_xfield = f"ocean_x_scircle{rtvInfoObj.grid_search_radius:.0f}km"
        scircle_yfield = f"ocean_y_scircle{rtvInfoObj.grid_search_radius:.0f}km"

    elif rtvInfoObj.grid_search_radius * 1000 == np.floor(rtvInfoObj.grid_search_radius * 1000):
        scircle_xfield = f"ocean_x_scircle{rtvInfoObj.grid_search_radius * 1000:.0f}m"
        scircle_yfield = f"ocean_y_scircle{rtvInfoObj.grid_search_radius * 1000:.0f}m"

    else:
        msg = (f"Invalid grid search radius of {rtvInfoObj.grid_search_radius}"
               "km. Value must be a whole number when represented in meters")
        logging.error(msg)
        raise ValueError(msg)

    # Compute total solutions
    if rtvInfoObj.method == 'uwls':
        TotalComp = [UwlsTotalInfo() for i in range(nArrayLen)]
    elif rtvInfoObj.method == 'oi':
        TotalComp = [OITotalInfo() for i in range(nArrayLen)]
    scircle_xfield_arr = getattr(compTotInput.grid, f"{scircle_xfield}")
    scircle_yfield_arr = getattr(compTotInput.grid, f"{scircle_yfield}")

    nPoints = len(sPoint)

    # Loop over each potential solution grid point
    for iPoint in range(nPoints):

        # index of the solution point
        solutionIndex = int(sPoint[iPoint])
        scLat = scircle_yfield_arr[:, solutionIndex]
        scLon = scircle_xfield_arr[:, solutionIndex]

        nSitesContributing = 0
        containsNewData = False
        rpSpeed = []
        rpHeading = []
        rpLon = []
        rpLat = []
        gridloc = np.zeros((1, 2))
        gridloc[0, 0] = compTotInput.grid.ocean_xy[0][solutionIndex]
        gridloc[0, 1] = compTotInput.grid.ocean_xy[1][solutionIndex]

        # Loop over each site to find radials within the grid point's search radius
        for iSite in range(nSites):

            currRadial = radials[iSite]
            inPolPoints = inpolygon(currRadial.longitude, currRadial.latitude, scLon, scLat)

            if np.any(inPolPoints):
                rpSpeed.append(currRadial.speed[inPolPoints])
                rpHeading.append(currRadial.heading[inPolPoints])
                rpLon.append(currRadial.longitude[inPolPoints])
                rpLat.append(currRadial.latitude[inPolPoints])
                nSitesContributing = nSitesContributing + 1

                if (not containsNewData) and (currRadial.isNew):
                    containsNewData = True

        if len(rpSpeed) > 0:
            rpSpeed = np.concatenate(rpSpeed)
            rpHeading = np.concatenate(rpHeading)
            rpLon = np.concatenate(rpLon)
            rpLat = np.concatenate(rpLat)

        # See if we have:
        #  (1) New radial data with
        #  (2) enough contributing sites and
        #  (3) enough radials to compute a total
        haveEnoughSitesContributing = nSitesContributing >= rtvInfoObj.min_rad_sites
        haveEnoughRadials = len(rpSpeed) >= rtvInfoObj.min_radials

        if containsNewData and haveEnoughSitesContributing and haveEnoughRadials:
            # Compute total
            if rtvInfoObj.method == 'uwls':
                TotalComp[solutionIndex] = uwlsTotal(rpSpeed, rpHeading)
            elif rtvInfoObj.method == 'oi':
                TotalComp[solutionIndex] = oiTotal(rpSpeed, rpHeading, rpLon, rpLat, gridloc, rtvInfoObj.oi_mdlvar, rtvInfoObj.oi_errvar,rtvInfoObj.oi_sx,rtvInfoObj.oi_sy,rtvInfoObj.oi_weighting)
            else:
                TotalComp[solutionIndex] = []
            nRads[solutionIndex] = len(rpSpeed)
            nRadSites[solutionIndex] = nSitesContributing

    u_temp = np.array([tot.u for tot in TotalComp])
    v_temp = np.array([tot.v for tot in TotalComp])
    dopx_temp = np.array([tot.dopx for tot in TotalComp])
    dopy_temp = np.array([tot.dopy for tot in TotalComp])
    hdop_temp = np.array([tot.hdop for tot in TotalComp])
    covr_temp = np.array([tot.covr for tot in TotalComp])

    # Leave for now, but can call other hfradarpy routine to do the following QARTOD QC...
    # Filter total solutions: infinite, complex, speed threshold, HDOP threshold
    iInf = np.logical_or(np.isinf(u_temp), np.isinf(v_temp))
    iCpx = np.logical_or(np.logical_or(u_temp.imag > 0, v_temp.imag > 0),
                         np.logical_or(dopx_temp.imag > 0, dopy_temp.imag > 0))
    iSpd = np.sqrt(u_temp ** 2 + v_temp ** 2) > rtvInfoObj.max_rtv_speed
    iHdop = hdop_temp > rtvInfoObj.uwls_max_hdop
    mask = np.logical_or(np.logical_or(iInf, iCpx), np.logical_or(iSpd, iHdop))

    # Mask the sites
    u_temp[mask] = np.nan
    v_temp[mask] = np.nan
    dopx_temp[mask] = np.nan
    dopy_temp[mask] = np.nan
    hdop_temp[mask] = np.nan
    covr_temp[mask] = np.nan
    nRads[mask] = 0
    nRadSites[mask] = 0

    if np.any(mask):
        msg = (f"Masked {np.sum(mask)} total solutions\n"
               f"{np.sum(iInf)} inf, {np.sum(iCpx)} complex, "
               f"{np.sum(iSpd)} speed, {np.sum(iHdop)} hdop")
        logging.info(msg)
    else:
        logging.debug('No solutions eliminated by masking')

    # Structure data
    if np.any(nRads) > 0:
        U_totals.lat = compTotInput.grid.ocean_xy[1]
        U_totals.lon = compTotInput.grid.ocean_xy[0]
        U_totals.u_xvel = u_temp
        U_totals.v_yvel = v_temp
        U_totals.dopx = dopx_temp
        U_totals.dopy = dopy_temp
        U_totals.hdop = hdop_temp
        U_totals.covr = covr_temp
        U_totals.nRads = nRads
        U_totals.nSites = nRadSites
        U_totals.grid.resolution_km = compTotInput.grid.resolution_km
        U_totals.grid.projection = compTotInput.grid.projection
        U_totals.grid.x_range = compTotInput.grid.x_range
        U_totals.grid.y_range = compTotInput.grid.y_range
        U_totals.grid.dx = compTotInput.grid.dx
        U_totals.grid.dy = compTotInput.grid.dy
        U_totals.grid.size = compTotInput.grid.size
        U_totals.grid.ocean_indices = compTotInput.grid.ocean_indices
        U_totals.grid.ocean_xy = compTotInput.grid.ocean_xy

    return U_totals
