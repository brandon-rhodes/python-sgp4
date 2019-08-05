"""The Satellite class."""

from sgp4.ext import jday
try:
    from sgp4.cpropagation import sgp4
except ImportError:
    from sgp4.propagation import sgp4

minutes_per_day = 1440.

class Satellite(object):
    """An earth-orbiting satellite as represented by the SGP4 model.

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
    ``no``
        Mean motion in radians per minute.

    """
    def propagate(self, year, month=1, day=1, hour=0, minute=0, second=0.0):
        """Return a position and velocity vector for a given date and time."""

        j = jday(year, month, day, hour, minute, second)
        m = (j - self.jdsatepoch) * minutes_per_day
        r, v = sgp4(self, m)
        return r, v
