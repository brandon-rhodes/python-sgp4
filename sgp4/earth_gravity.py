"""Three earth-gravity models for use with SGP4."""

from collections import namedtuple
from sgp4.propagation import getgravconst

EarthGravity = namedtuple(
    'EarthGravity',
    'tumin mu radiusearthkm xke j2 j3 j4 j3oj2',
    )

wgs72old = EarthGravity(*getgravconst('wgs72old'))
wgs72 = EarthGravity(*getgravconst('wgs72'))
wgs84 = EarthGravity(*getgravconst('wgs84'))
