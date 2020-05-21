"""Various conveniences.

Higher-level libraries like Skyfield that use this one usually have
their own date and time handling.  But for folks using this library by
itself, native Python datetime handling could be convenient.

"""
import datetime as dt
from .functions import days2mdhms, jday

class _UTC(dt.tzinfo):
    'UTC'
    zero = dt.timedelta(0)
    def __repr__(self):
        return 'UTC'
    def dst(self, datetime):
        return self.zero
    def tzname(self, datetime):
        return 'UTC'
    def utcoffset(self, datetime):
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

def sat_epoch_datetime(sat):
    """Return the epoch of the given satellite as a Python datetime."""
    year = sat.epochyr
    year += 1900 + (year < 57) * 100
    days = sat.epochdays
    month, day, hour, minute, second = days2mdhms(year, days)
    if month == 12 and day > 31:  # for that time the ISS epoch was "Dec 32"
        year += 1
        month = 1
        day -= 31
    second, fraction = divmod(second, 1.0)
    second = int(second)
    micro = int(fraction * 1e6)
    return dt.datetime(year, month, day, hour, minute, second, micro, UTC)
