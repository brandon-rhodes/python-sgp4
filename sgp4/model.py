"""The Satellite class."""

from sgp4.ext import jday
from sgp4.propagation import sgp4

minutes_per_day = 1440.

class Satellite(object):
    """An earth-orbiting satellite as represented by the SGP4 model."""

    # TODO: document commonly-used attributes here

    def propagate(self, year, month=1, day=1, hour=0, minute=0, second=0.0):
        """Return a position and velocity vector for a given date and time."""

        j = jday(year, month, day, hour, minute, second)
        m = (j - self.jdsatepoch) * minutes_per_day
        r, v = sgp4(self, m)
        return r, v
