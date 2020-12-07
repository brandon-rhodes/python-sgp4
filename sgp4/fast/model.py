from numba import typeof
from numba.core.types import float64, int32, int64, string
from numba.experimental import jitclass

from sgp4.earth_gravity import wgs72old, wgs72, wgs84
from sgp4.ext import invjday, jday
from sgp4.io import twoline2rv as io_twoline2rv
from sgp4.fast.propagation import sgp4, sgp4init

WGS72OLD = 0
WGS72 = 1
WGS84 = 2
gravity_constants = wgs72old, wgs72, wgs84  # indexed using enum values above
minutes_per_day = 1440.

spec = [
    ('Om', float64),
    ('a', float64),
    ('alta', float64),
    ('altp', float64),
    ('am', float64),
    ('argpdot', float64),
    ('argpo', float64),
    ('atime', float64),
    ('aycof', float64),
    ('bstar', float64),
    ('cc1', float64),
    ('cc4', float64),
    ('cc5', float64),
    ('classification', string),
    ('con41', float64),
    ('d2', float64),
    ('d2201', float64),
    ('d2211', float64),
    ('d3', float64),
    ('d3210', float64),
    ('d3222', float64),
    ('d4', float64),
    ('d4410', float64),
    ('d4422', float64),
    ('d5220', float64),
    ('d5232', float64),
    ('d5421', float64),
    ('d5433', float64),
    ('dedt', float64),
    ('del1', float64),
    ('del2', float64),
    ('del3', float64),
    ('delmo', float64),
    ('didt', float64),
    ('dmdt', float64),
    ('dnodt', float64),
    ('domdt', float64),
    ('e3', float64),
    ('ecco', float64),
    ('ee2', float64),
    ('elnum', float64),
    ('em', float64),
    ('ephtype', string),  # TODO: ephtype is int in C++, but string in Python!
    ('epochdays', float64),
    ('epochyr', int32),
    ('error', int32),
    ('error_message', string),
    ('eta', float64),
    ('gsto', float64),
    ('im', float64),
    ('inclo', float64),
    ('init', string),
    ('intldesg', string),
    ('irez', int32),
    ('isimp', int32),
    ('j2', float64),
    ('j3', float64),
    ('j3oj2', float64),
    ('j4', float64),
    ('jdsatepoch', float64),
    ('mdot', float64),
    ('method', string),
    ('mm', float64),
    ('mo', float64),
    ('mu', float64),
    ('nddot', float64),
    ('ndot', float64),
    ('nm', float64),
    ('no_kozai', float64),
    ('no_unkozai', float64),
    ('nodecf', float64),
    ('nodedot', float64),
    ('nodeo', float64),
    ('om', float64),
    ('omgcof', float64),
    ('operationmode', string),
    ('peo', float64),
    ('pgho', float64),
    ('pho', float64),
    ('pinco', float64),
    ('plo', float64),
    ('radiusearthkm', float64),
    ('revnum', string),  # TODO: revnum is long in C++, but string in Python!
    ('satnum', int64),
    ('se2', float64),
    ('se3', float64),
    ('sgh2', float64),
    ('sgh3', float64),
    ('sgh4', float64),
    ('sh2', float64),
    ('sh3', float64),
    ('si2', float64),
    ('si3', float64),
    ('sinmao', float64),
    ('sl2', float64),
    ('sl3', float64),
    ('sl4', float64),
    ('t', float64),
    ('t2cof', float64),
    ('t3cof', float64),
    ('t4cof', float64),
    ('t5cof', float64),
    ('tumin', float64),
    ('whichconst', typeof(wgs84)),
    ('x1mth2', float64),
    ('x7thm1', float64),
    ('xfact', float64),
    ('xgh2', float64),
    ('xgh3', float64),
    ('xgh4', float64),
    ('xh2', float64),
    ('xh3', float64),
    ('xi2', float64),
    ('xi3', float64),
    ('xke', float64),
    ('xl2', float64),
    ('xl3', float64),
    ('xl4', float64),
    ('xlamo', float64),
    ('xlcof', float64),
    ('xli', float64),
    ('xmcof', float64),
    ('xni', float64),
    ('zmol', float64),
    ('zmos', float64),
    ('jdsatepochF', float64)
]


@jitclass(spec)
class Satrec:

    @property
    def no(self):
        return self.no_kozai

    def __init__(self):
        # Never called directly, but numba requires it
        # Therefore, we do not do anything
        pass

    def sgp4init(self, whichconst, opsmode, satnum, epoch, bstar,
                 ndot, nddot, ecco, argpo, inclo, mo, no_kozai, nodeo):
        whichconst = gravity_constants[whichconst]

        y, m, d, H, M, S = invjday(epoch + 2433281.5)
        jan0epoch = jday(y, 1, 0, 0, 0, 0.0) - 2433281.5

        self.epochyr = y % 1000
        self.epochdays = epoch - jan0epoch
        self.jdsatepoch, self.jdsatepochF = divmod(epoch, 1.0)
        self.jdsatepoch += 2433281.5

        sgp4init(whichconst, opsmode, satnum, epoch, bstar, ndot, nddot,
                 ecco, argpo, inclo, mo, no_kozai, nodeo, self)

    def sgp4(self, jd, fr):
        tsince = ((jd - self.jdsatepoch) * minutes_per_day +
                  (fr - self.jdsatepochF) * minutes_per_day)
        r, v = sgp4(self, tsince)
        return self.error, r, v

    def sgp4_tsince(self, tsince):
        r, v = sgp4(self, tsince)
        return self.error, r, v


# Separate function instead of @classmethod
# because @jitclass does not support classmethods,
# see https://github.com/numba/numba/issues/4688
def twoline2rv(satrec, line1, line2, whichconst=WGS72):
    whichconst = gravity_constants[whichconst]
    io_twoline2rv(line1, line2, whichconst, 'i', satrec, init=False)
    sgp4init(whichconst, 'i', satrec.satnum, satrec.jdsatepoch - 2433281.5, satrec.bstar,
             satrec.ndot, satrec.nddot, satrec.ecco, satrec.argpo, satrec.inclo, satrec.mo,
             satrec.no_kozai, satrec.nodeo, satrec)

    # Install a fancy split JD of the kind the C++ natively supports.
    # We rebuild it from the TLE year and day to maintain precision.
    year = satrec.epochyr
    days, fraction = divmod(satrec.epochdays, 1.0)
    satrec.jdsatepoch = year * 365 + (year - 1) // 4 + days + 1721044.5
    satrec.jdsatepochF = round(fraction, 8)  # exact number of digits in TLE

    # Undo my non-standard 4-digit year
    satrec.epochyr %= 100
    return satrec
