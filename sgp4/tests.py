"""Test suite for SPG4."""

import os
from unittest import TestCase
from sgp4.io import twoline2rv

class SatRec(object):
    pass


class Tests(TestCase):

    def test_tle(self):
        tlepath = os.path.join(os.path.dirname(__file__), 'SGP4-VER.TLE')
        with open(tlepath) as tlefile:
            lines = tlefile.readlines()

        for i in range(len(lines)):

            if not lines[i].startswith('1'):
                continue

            line1 = lines[i]
            line2 = lines[i + 1]

            satrec = SatRec()

            twoline2rv(line1, line2, 'v', None, 'i', 'wgs84', satrec)
            break
