"""Test suite for SPG4."""

import os
from unittest import TestCase
from sgp4.io import twoline2rv
from sgp4.propagation import sgp4

class SatRec(object):
    pass


class Tests(TestCase):

    def test_tle(self):
        dirpath = os.path.dirname(__file__)
        tlepath = os.path.join(dirpath, 'SGP4-VER.TLE')
        with open(tlepath) as tlefile:
            lines = tlefile.readlines()

        whichconst = 'wgs72'

        for i in range(len(lines)):

            if not lines[i].startswith('1'):
                continue

            line1 = lines[i]
            line2 = lines[i + 1]

            satrec = SatRec()

            twoline2rv(line1, line2, 'c', None, 'i', whichconst, satrec)

            # Call the propagator to get the initial state vector value.
            ro = [0.0, 0.0, 0.0]
            vo = [0.0, 0.0, 0.0]
            sgp4(whichconst, satrec,  0.0, ro,  vo);
            print satrec.t, ro[0], ro[1], ro[2], vo[0], vo[1], vo[2]

            tcppath = os.path.join(dirpath, 'tcppver.out')
            with open(tcppath) as tcpfile:
                tcplines = tcpfile.readlines()

            print ' '.join(repr(float(field)) for field
                           in tcplines[1].split()[:7])

            sgp4(whichconst, satrec,  360.0, ro,  vo);
            print satrec.t, ro[0], ro[1], ro[2], vo[0], vo[1], vo[2]
            print ' '.join(repr(float(field)) for field
                           in tcplines[2].split()[:7])

            break
