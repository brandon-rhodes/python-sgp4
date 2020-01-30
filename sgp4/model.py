"""The Satellite class."""

from sgp4.earth_gravity import wgs72
from sgp4.ext import jday
from sgp4.io import twoline2rv
from sgp4.propagation import sgp4

minutes_per_day = 1440.

class Satrec(object):
    """Slow Python-only version of the satellite object."""

    jdsatepochF = 0.0  # for compatibility with accelerated version

    @classmethod
    def twoline2rv(cls, line1, line2):
        self = cls()
        twoline2rv(line1, line2, wgs72, 'i', self)
        return self

    def sgp4(self, jd, fr):
        tsince = (jd - self.jdsatepoch + fr) * minutes_per_day
        r, v = sgp4(self, tsince)
        return self.error, r, v

class SatrecArray(object):
    """Slow Python-only version of the satellite array."""

    def __init__(self, satrecs):
        self._satrecs = satrecs
        # Cache optional import that we now know we need.
        from numpy import array
        self.array = array

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
