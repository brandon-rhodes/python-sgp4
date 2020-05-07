"""Various conveniences.

Higher-level libraries like Skyfield that use this one usually have
their own date and time handling.  But for folks using this library by
itself, native Python datetime handling could be convenient.

"""
import datetime as dt
from .functions import jday

class _UTC(dt.tzinfo):
     'UTC'
     zero = dt.timedelta(0)
     def utcoffset(self, datetime):
          return self.zero
     def tzname(self, datetime):
          return 'UTC'
     def dst(self, datetime):
          return self.zero

UTC = _UTC()

def jday_datetime(datetime):
    """Return two floats that, when added, produce the specified Julian date.

    The first float returned gives the date, while the second float
    provides an additional offset for the particular hour, minute, and
    second of that date.  Because the second float is much smaller in
    magnitude it can, unlike the first float, be accurate down to very
    small fractions of a second.

    >>> jd, fr = jday(2020, 2, 11, 13, 57, 0)
    >>> jd
    2458890.5
    >>> fr
    0.58125

    Note that the first float, which gives the moment of midnight that
    commences the given calendar date, always carries the fraction
    ``.5`` because Julian dates begin and end at noon.  This made Julian
    dates more convenient for astronomers in Europe, by making the whole
    night belong to a single Julian date.

    The input is a native `datetime` object. Timezone of the input is
    converted internally to UTC.

    """
    u = datetime.astimezone(UTC)
    year = u.year
    mon = u.month
    day = u.day
    hr = u.hour
    minute = u.minute
    sec = u.second + u.microsecond * 1e-6

    return jday(year, mon, day, hr, minute, sec)
