"""General-purpose routines.

It seemed a shame for the ``api`` module to have to import large legacy
modules to offer simple date handling, so this small module holds the
routines instead.

"""
from math import trunc

def jday(year, mon, day, hr, minute, sec):
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
    ``.5`` because Julian dates begin and end at noon.  This made
    Julian dates more convenient for astronomers in Europe, by making
    the whole night belong to a single Julian date.

    This function is a simple translation to Python of the C++ routine
    ``jday()`` in Vallado's ``SGP4.cpp``.

    """
    jd = (367.0 * year
         - 7 * (year + ((mon + 9) // 12.0)) * 0.25 // 1.0
	   + 275 * mon / 9.0 // 1.0
	   + day
         + 1721013.5)
    fr = (sec + minute * 60.0 + hr * 3600.0) / 86400.0;
    return jd, fr

def days2mdhms(year, days, round_to_microsecond=6):
    """Convert a float point number of days into the year into date and time.

    Given the integer year plus the "day of the year" (where 1.0 means
    the beginning of January 1, 2.0 means the beginning of January 2,
    and so forth), returns the Gregorian calendar month, day, hour,
    minute, and floating point seconds.

    >>> days2mdhms(2000, 1.0)   # January 1
    (1, 1, 0, 0, 0.0)
    >>> days2mdhms(2000, 32.0)  # February 1
    (2, 1, 0, 0, 0.0)
    >>> days2mdhms(2000, 366.0)  # December 31, since 2000 was a leap year
    (12, 31, 0, 0, 0.0)

    """
    whole, fraction = divmod(days, 1.0)

    is_leap = year % 400 == 0 or (year % 4 == 0 and year % 100 != 0)
    month, day = _day_of_year_to_month_day(int(whole), is_leap)
    if month == 13:  # behave like the original in case of overflow
        month = 12
        day += 31

    # The 8 digits of floating point day specified in the TLE have a
    # resolution of exactly 1e-8 * 24 * 3600 * 1e6 = 864 microseconds,
    # so round off any floating-point noise beyond the microsecond.
    if round_to_microsecond:
        fraction += 0.5 / 86400e6

    second = fraction * 86400.0
    minute, second = divmod(second, 60.0)
    hour, minute = divmod(minute, 60.0)

    if round_to_microsecond:
        second = trunc(second * 1e6) / 1e6

    return month, day, int(hour), int(minute), second

def _day_of_year_to_month_day(day_of_year, is_leap):
    """Core logic for turning days into months, for easy testing."""
    february_bump = (2 - is_leap) * (day_of_year >= 60 + is_leap)
    august = day_of_year >= 215
    month, day = divmod(2 * (day_of_year - 1 + 30 * august + february_bump), 61)
    month += 1 - august
    day //= 2
    day += 1
    return month, day
