"""The Satellite class."""

from sgp4.earth_gravity import wgs72old, wgs72, wgs84
from sgp4.ext import jday
from sgp4.io import twoline2rv
from sgp4.propagation import sgp4

WGS72OLD = 0
WGS72 = 1
WGS84 = 2
gravity_constants = wgs72old, wgs72, wgs84  # indexed using enum values above
minutes_per_day = 1440.

class Satrec(object):
    """Slow Python-only version of the satellite object."""

    # Approximate the behavior of the C-accelerated class by locking
    # down attribute access, to avoid folks accidentally writing code
    # against this class and adding extra attributes, then moving to a
    # computer where the C-accelerated class is used and having their
    # code suddenly produce errors.
    __slots__ = (
        'Om', 'a', 'alta', 'altp', 'am', 'argpdot', 'argpo', 'atime', 'aycof',
        'bstar', 'cc1', 'cc4', 'cc5', 'classification', 'con41', 'd2', 'd2201',
        'd2211', 'd3', 'd3210', 'd3222', 'd4', 'd4410', 'd4422', 'd5220',
        'd5232', 'd5421', 'd5433', 'dedt', 'del1', 'del2', 'del3', 'delmo',
        'didt', 'dmdt', 'dnodt', 'domdt', 'e3', 'ecco', 'ee2', 'elnum', 'em',
        'ephtype', 'epoch', 'epochdays', 'epochyr', 'error', 'error_message',
        'eta', 'gsto', 'im', 'inclo', 'init', 'intldesg', 'irez', 'isimp',
        'j2', 'j3', 'j3oj2', 'j4', 'jdsatepoch', 'mdot', 'method', 'mm', 'mo',
        'mu', 'nddot', 'ndot', 'nm', 'no_kozai', 'no_unkozai', 'nodecf',
        'nodedot', 'nodeo', 'om', 'omgcof', 'operationmode', 'peo', 'pgho',
        'pho', 'pinco', 'plo', 'radiusearthkm', 'revnum', 'satnum', 'se2',
        'se3', 'sgh2', 'sgh3', 'sgh4', 'sh2', 'sh3', 'si2', 'si3', 'sinmao',
        'sl2', 'sl3', 'sl4', 't', 't2cof', 't3cof', 't4cof', 't5cof', 'tumin',
        'whichconst', 'x1mth2', 'x7thm1', 'xfact', 'xgh2', 'xgh3', 'xgh4',
        'xh2', 'xh3', 'xi2', 'xi3', 'xke', 'xl2', 'xl3', 'xl4', 'xlamo',
        'xlcof', 'xli', 'xmcof', 'xni', 'zmol', 'zmos',
    )

    array = None       # replaced, if needed, with NumPy array()
    jdsatepochF = 0.0  # for compatibility with accelerated version

    @property
    def no(self):
        return self.no_kozai

    @classmethod
    def twoline2rv(cls, line1, line2, whichconst=WGS72):
        whichconst = gravity_constants[whichconst]
        self = cls()
        twoline2rv(line1, line2, whichconst, 'i', self)
        self.epochyr %= 100  # undo my non-standard 4-digit year
        return self

    def sgp4(self, jd, fr):
        tsince = (jd - self.jdsatepoch + fr) * minutes_per_day
        r, v = sgp4(self, tsince)
        return self.error, r, v

    def sgp4_array(self, jd, fr):
        """Compute positions and velocities for the times in a NumPy array.

        Given NumPy arrays ``jd`` and ``fr`` of the same length that
        supply the whole part and the fractional part of one or more
        Julian dates, return a tuple ``(e, r, v)`` of three vectors:

        * ``e``: nonzero for any dates that produced errors, 0 otherwise.
        * ``r``: position vectors in kilometers.
        * ``v``: velocity vectors in kilometers per second.

        """
        # Import NumPy the first time sgp4_array() is called.
        array = self.array
        if array is None:
            from numpy import array
            Satrec.array = array

        results = []
        z = list(zip(jd, fr))
        for jd_i, fr_i in z:
            results.append(self.sgp4(jd_i, fr_i))
        elist, rlist, vlist = zip(*results)

        e = self.array(elist)
        r = self.array(rlist)
        v = self.array(vlist)

        r.shape = v.shape = len(jd), 3
        return e, r, v

class SatrecArray(object):
    """Slow Python-only version of the satellite array."""

    __slots__ = ('_satrecs',)

    array = None  # replaced with NumPy array(), if the user tries calling

    def __init__(self, satrecs):
        self._satrecs = satrecs
        # Import NumPy the first time a SatrecArray is instantiated.
        if self.array is None:
            from numpy import array
            SatrecArray.array = array

    def sgp4(self, jd, fr):
        """Compute positions and velocities for the satellites in this array.

        Given NumPy scalars or arrays ``jd`` and ``fr`` supplying the
        whole part and the fractional part of one or more Julian dates,
        return a tuple ``(e, r, v)`` of three vectors that are each as
        long as ``jd`` and ``fr``:

        * ``e``: nonzero for any dates that produced errors, 0 otherwise.
        * ``r``: (x,y,z) position vector in kilometers.
        * ``v``: (dx,dy,dz) velocity vector in kilometers per second.

        """
        results = []
        z = list(zip(jd, fr))
        for satrec in self._satrecs:
            for jd_i, fr_i in z:
                results.append(satrec.sgp4(jd_i, fr_i))
        elist, rlist, vlist = zip(*results)

        e = self.array(elist)
        r = self.array(rlist)
        v = self.array(vlist)

        jdlen = len(jd)
        mylen = len(self._satrecs)
        e.shape = (mylen, jdlen)
        r.shape = v.shape = (mylen, jdlen, 3)

        return e, r, v

class Satellite(object):
    """The old Satellite object for compatibility with sgp4 1.x.

    Most of this class's hundred-plus attributes are intermediate values
    of interest only to the propagation algorithm itself.  Here are the
    attributes set by ``sgp4.io.twoline2rv()`` in which users are likely
    to be interested:

    ``satnum``
        Unique satellite number given in the TLE file.
    ``epochyr``
        Full four-digit year of this element set's epoch moment.
    ``epochdays``
        Fractional days into the year of the epoch moment.
    ``jdsatepoch``
        Julian date of the epoch (computed from ``epochyr`` and ``epochdays``).
    ``ndot``
        First time derivative of the mean motion (ignored by SGP4).
    ``nddot``
        Second time derivative of the mean motion (ignored by SGP4).
    ``bstar``
        Ballistic drag coefficient B* in inverse earth radii.
    ``inclo``
        Inclination in radians.
    ``nodeo``
        Right ascension of ascending node in radians.
    ``ecco``
        Eccentricity.
    ``argpo``
        Argument of perigee in radians.
    ``mo``
        Mean anomaly in radians.
    ``no_kozai``
        Mean motion in radians per minute.

    """
    jdsatepochF = 0.0  # for compatibility with new Satrec; makes tests simpler

    # TODO: only offer this on legacy class we no longer document
    def propagate(self, year, month=1, day=1, hour=0, minute=0, second=0.0):
        """Return a position and velocity vector for a given date and time."""

        j = jday(year, month, day, hour, minute, second)
        m = (j - self.jdsatepoch) * minutes_per_day
        r, v = sgp4(self, m)
        return r, v

    @property
    def no(self):
        """Support renamed attribute for any code still using the old name."""
        return self.no_kozai
