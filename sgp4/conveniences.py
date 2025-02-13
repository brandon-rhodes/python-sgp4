"""Various conveniences.

Higher-level libraries like Skyfield that use this one usually have
their own date and time handling.  But for folks using this library by
itself, native Python datetime handling could be convenient.

"""
import datetime as dt
import sgp4
from math import pi
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

_ATTRIBUTES = []
_ATTR_MAXES = {}
_MAX_VALUES = {'2pi': 2*pi, 'pi': pi}

def _load_attributes():
    for line in sgp4.__doc__.splitlines():
        if line.endswith('*'):
            title = line.strip('*')
            _ATTRIBUTES.append(title)
        elif line.startswith('| ``'):
            pieces = line.split('``')
            name = pieces[1]
            _ATTRIBUTES.append(name)
            i = 2
            while pieces[i] == ', ':
                another_name = pieces[i+1]
                _ATTRIBUTES.append(another_name)
                i += 2
            if '<' in line:
                _ATTR_MAXES[name] = '2pi' if ('2pi' in line) else 'pi'

def check_satrec(sat):
    """Check whether satellite orbital elements are within range."""

    if not _ATTRIBUTES:
        _load_attributes()

    e = []

    for name, max_name in sorted(_ATTR_MAXES.items()):
        value = getattr(sat, name)
        if 0.0 <= value < _MAX_VALUES[max_name]:
            continue
        e.append('  {0} = {1:f} is outside the range 0 <= {0} < {2}\n'
                 .format(name, value, max_name))

    if e:
        raise ValueError('satellite parameters out of range:\n' + '\n'.join(e))


def dump_satrec(sat, sat2=None):
    """Yield lines that list the attributes of one or two satellites."""

    if not _ATTRIBUTES:
        _load_attributes()

    for item in _ATTRIBUTES:
        if item[0].isupper():
            title = item
            yield '\n'
            yield '# -------- {0} --------\n'.format(title)
        else:
            name = item
            value = getattr(sat, item, '(not set)')
            line = '{0} = {1!r}\n'.format(item, value)
            if sat2 is not None:
                value2 = getattr(sat2, name, '(not set)')
                verdict = '==' if (value == value2) else '!='
                line = '{0:39} {1} {2!r}\n'.format(line[:-1], verdict, value2)
            yield line
