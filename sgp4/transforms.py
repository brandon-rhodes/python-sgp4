"""The coordinate transform functions from satellite.js.

https://github.com/shashwatak/satellite-js

"""

from math import atan2, cos, asin, pi, sin, sqrt

rad2deg = 180 / pi
twopi = 2.0 * pi

def _radiansToDegrees(radians):
    return radians * rad2deg

def degreesLat(radians):
    if (pi / 2) < radians < (-pi / 2):
        raise ValueError('Latitude radians must be in range [-pi/2 pi/2].')
    return _radiansToDegrees(radians)

def degreesLong(radians):
    if pi < radians < -pi:
        raise ValueError('Longitude radians must be in range [-pi pi].')
    return _radiansToDegrees(radians)

def geodeticToEcf(longitude, latitude, height=0):
    """Return ECF converted from Geodetic."""

    a = 6378.137
    b = 6356.7523142
    f = (a - b) / a
    e2 = ((2 * f) - (f * f))
    normal = a / sqrt(1 - (e2 * (sin(latitude) * sin(latitude))))

    x = (normal + height) * cos(latitude) * cos(longitude)
    y = (normal + height) * cos(latitude) * sin(longitude)
    z = ((normal * (1 - e2)) + height) * sin(latitude)

    return (x, y, z)

def eciToGeodetic(x, y, z, gmst):
    """Return Geodetic converted from ECI.
    
    http://www.celestrak.com/columns/v02n03/
    """
    
    a = 6378.137
    b = 6356.7523142
    R = sqrt((x * x) + (y * y))
    f = (a - b) / a
    e2 = ((2 * f) - (f * f))

    longitude = atan2(y, x) - gmst

    while longitude < -pi:
        longitude += twoPi
    while longitude > pi:
        longitude -= twoPi

    kmax = 20
    k = 0

    latitude = atan2(z, sqrt((x * x) + (y * y)))
    while k < kmax:
        C = 1 / sqrt(1 - (e2 * (sin(latitude) * sin(latitude))))
        latitude = atan2(z + (a * C * e2 * sin(latitude)), R)
        k += 1

    height = (R / cos(latitude)) - (a * C)
    return (longitude, latitude, height)

def ecfToEci(x, y, z, gmst):
    """Return ECI converted from ECF.
    
    ccar.colorado.edu/ASEN5070/handouts/coordsys.doc
    
    [X]     [C -S  0][X]
    [Y]  =  [S  C  0][Y]
    [Z]eci  [0  0  1][Z]ecf
    """
    
    X = (x * cos(gmst)) - (y * sin(gmst))
    Y = (x * (sin(gmst))) + (y * cos(gmst))
    Z = z

    return (X, Y, Z)

def eciToEcf(x, y, z, gmst):
    """Return ECF converted from ECI.
    
    ccar.colorado.edu/ASEN5070/handouts/coordsys.doc
    
    [X]     [C -S  0][X]
    [Y]  =  [S  C  0][Y]
    [Z]eci  [0  0  1][Z]ecf
    
    Inverse:
    [X]     [C  S  0][X]
    [Y]  =  [-S C  0][Y]
    [Z]ecf  [0  0  1][Z]eci
    """

    fx = (x * cos(gmst)) + (y * sin(gmst))
    fy = (x * (-sin(gmst))) + (y * cos(gmst))
    fz = z

    return (fx, fy, fz)

def _topocentric(longitude, latitude, x, y, z):
    """Return the topocentric from the observer to the satellite.
    
    http://www.celestrak.com/columns/v02n02/
    TS Kelso's method, except I'm using ECF frame and he uses ECI.
    """

    ox, oy, oz = geodeticToEcf(longitude, latitude)

    rx = x - ox
    ry = y - oy
    rz = z - oz

    topS = ((sin(latitude) * cos(longitude) * rx) \
            + (sin(latitude) * sin(longitude) * ry)) \
            - (cos(latitude) * rz)
    topE = (-sin(longitude) * rx) + (cos(longitude) * ry)

    topZ = (cos(latitude) * cos(longitude) * rx) \
        + (cos(latitude) * sin(longitude) * ry) \
        + (sin(latitude) * rz)

    return (topS, topE, topZ)

def _topocentricToLookAngles(topS, topE, topZ):
    """Return the look angles from topocentric."""

    rangeSat = sqrt((topS * topS) + (topE * topE) + (topZ * topZ))
    El = asin(topZ / rangeSat)
    Az = atan2(-topE, topS) + pi

    return (Az, El, rangeSat)

def ecfToLookAngles(longitude, latitude, x, y, z):
    """Return the look angles from the observer to the satellite."""
    
    topS, topE, topZ = _topocentric(longitude, latitude, x, y, z)
    return _topocentricToLookAngles(topS, topE, topZ)
