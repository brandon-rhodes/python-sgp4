"""Test suite for SPG4."""

import os
from unittest import TestCase

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

            print line1
            break
