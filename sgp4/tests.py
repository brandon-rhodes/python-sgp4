"""Test suite for SGP4."""

import os
import re
import sys
from doctest import DocTestSuite, ELLIPSIS
from unittest import TestCase
from math import pi, isnan

from sgp4.earth_gravity import wgs72
from sgp4.ext import invjday, newtonnu, rv2coe
from sgp4.io import twoline2rv
from sgp4.propagation import sgp4

thisdir = os.path.dirname(__file__)
error = 2e-7
rad = 180.0 / pi


class Tests(TestCase):

    def test_tle_verify(self):
        # Check whether a test run produces the output in tcppver.out

        whichconst = 'wgs72'
        error_list = []
        actual = generate_test_output(whichconst, error_list)
        previous_data_line = None

        # Iterate across "tcppver.out", making sure that we ourselves
        # produce a line that looks very much like the corresponding
        # line in that file.

        tcppath = os.path.join(thisdir, 'tcppver.out')
        with open(tcppath) as tcpfile:
            for i, expected_line in enumerate(tcpfile, start = 1):

                try:
                    actual_line = next(actual)
                except StopIteration:
                    raise ValueError(
                        'WARNING: our output ended early, on line %d' % (i,))

                if actual_line == '(Use previous data line)':
                    actual_line = ('       0.00000000' +
                                   previous_data_line[17:107])

                # Compare the lines.  The first seven fields are printed
                # to very high precision, so we allow a small error due
                # to rounding differences; the rest are printed to lower
                # precision, and so can be compared textually.

                if 'xx' in actual_line:
                    similar = (actual_line == expected_line)
                else:
                    afields = actual_line.split()
                    efields = expected_line.split()
                    actual7 = [ float(a) for a in afields[:7] ]
                    expected7 = [ float(e) for e in efields[:7] ]
                    similar = (
                        len(actual7) == len(expected7)
                        and
                        all(
                            -error < (a - e) < error
                             for a, e in zip(actual7, expected7)
                             )
                        and
                        afields[7:] == efields[7:]  # just compare text
                        )

                if not similar:
                    raise ValueError(
                        'Line %d of output does not match:\n'
                        '\n'
                        'Expect: %s'
                        'Actual: %s'
                        % (i, expected_line, actual_line))

                if 'xx' not in actual_line:
                    previous_data_line = actual_line

        # Make sure the test file is not missing lines.

        missing_count = 0
        for actual_line in actual:
            missing_count += 1

        if missing_count > 0:
            raise ValueError('we produced %d extra lines' % (missing_count,))

        self.assertEqual(error_list, [
            (1, 'mean eccentricity -0.001329'
             ' not within range 0.0 <= e < 1.0'),
            (1, 'mean eccentricity -0.001208'
             ' not within range 0.0 <= e < 1.0'),
            (6, 'mrt 0.996159 is less than 1.0'
             ' indicating the satellite has decayed'),
            (6, 'mrt 0.996252 is less than 1.0'
             ' indicating the satellite has decayed'),
            (4, 'semilatus rectum -0.103223 is less than zero'),
            (3, 'perturbed eccentricity -122.217193'
             ' not within range 0.0 <= e <= 1.0'),
            (6, 'mrt 0.830534 is less than 1.0'
             ' indicating the satellite has decayed'),
            ])

    def test_hyperbolic_orbit(self):
        # Exercise the newtonnu() code path with asinh() to see whether
        # we can replace it with the one from Python's math module.

        self.assertEqual(newtonnu(1.0, 2.9),   # parabolic
                         (8.238092752965605, 194.60069989482898))

        self.assertEqual(newtonnu(1.1, 2.7),   # hyperbolic
                         (4.262200676156417, 34.76134082028372))

    def test_bad_first_line(self):
        with self.assertRaisesRegexp(ValueError, re.escape("""TLE format error

The Two-Line Element (TLE) format was designed for punch cards and is
therefore very strict about the position of every space and digit in a
TLE line.  Your line does not quite match.  Here is the official format
for line 1 followed by the line you provided:

1 NNNNNC NNNNNAAA NNNNN.NNNNNNNN +.NNNNNNNN +NNNNN-N +NNNNN-N N NNNNN
1 00005U 58002B   00179.78495062  .000000234 00000-0  28098-4 0  4753""")):
            twoline2rv(good1.replace('23 ', '234'), good2, wgs72)

    def test_bad_second_line(self):
        with self.assertRaisesRegexp(ValueError, re.escape("""TLE format error

The Two-Line Element (TLE) format was designed for punch cards and is
therefore very strict about the position of every space and digit in a
TLE line.  Your line does not quite match.  Here is the official format
for line 2 followed by the line you provided:

2 NNNNN NNN.NNNN NNN.NNNN NNNNNNN NNN.NNNN NNN.NNNN NN.NNNNNNNNNNNNNN
2 00005 34 .268234 8.7242 1859667 331.7664  19.3264 10.82419157413667""")):
            twoline2rv(good1, good2.replace(' 34', '34 '), wgs72)


good1 = '1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753'
good2 = '2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667'


def generate_test_output(whichconst, error_list):
    """Generate lines like those in the test file tcppver.out.

    This iterates through the satellites in "SGP4-VER.TLE", which are
    each supplemented with a time start/stop/step over which we are
    supposed to print results.

    """
    whichconst = wgs72
    tlepath = os.path.join(thisdir, 'SGP4-VER.TLE')
    with open(tlepath) as tlefile:
        tlelines = iter(tlefile.readlines())

    for line1 in tlelines:

        if not line1.startswith('1'):
            continue

        line2 = next(tlelines)
        satrec = twoline2rv(line1, line2, whichconst)

        yield '%ld xx\n' % (satrec.satnum,)

        for line in generate_satellite_output(satrec, line2, error_list):
            yield line


def generate_satellite_output(satrec, line2, error_list):
    """Print a data line for each time in line2's start/stop/step field."""

    mu = satrec.whichconst.mu

    r, v = sgp4(satrec, 0.0)
    if isnan(r[0]) and isnan(r[1]) and isnan(r[2]):
        error_list.append((satrec.error, satrec.error_message))
        yield '(Use previous data line)'
        return
    yield format_short_line(satrec, r, v)

    tstart, tend, tstep = (float(field) for field in line2[69:].split())

    tsince = tstart
    while tsince <= tend:
        if tsince == tstart == 0.0:
            tsince += tstep
            continue  # avoid duplicating the first line

        r, v = sgp4(satrec, tsince)

        if isnan(r[0]) and isnan(r[1]) and isnan(r[2]):
            error_list.append((satrec.error, satrec.error_message))
            return
        yield format_long_line(satrec, mu, r, v)

        tsince += tstep

    if tsince - tend < tstep - 1e-6:  # do not miss last line!
        r, v = sgp4(satrec, tend)
        if isnan(r[0]) and isnan(r[1]) and isnan(r[2]):
            error_list.append((satrec.error, satrec.error_message))
            return
        yield format_long_line(satrec, mu, r, v)


def format_short_line(satrec, r, v):
    """Short line, using the same format string that testcpp.cpp uses."""

    return ' %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f\n' % (
        satrec.t, r[0], r[1], r[2], v[0], v[1], v[2])


def format_long_line(satrec, mu, r, v):
    """Long line, using the same format string that testcpp.cpp uses."""

    short = format_short_line(satrec, r, v).strip('\n')

    jd = satrec.jdsatepoch + satrec.t / 1440.0
    year, mon, day, hr, minute, sec = invjday(jd)

    (p, a, ecc, incl, node, argp, nu, m, arglat, truelon, lonper
     ) = rv2coe(r, v, mu)

    return short + (' %14.6f %8.6f %10.5f %10.5f %10.5f %10.5f %10.5f'
                    ' %5i%3i%3i %2i:%2i:%9.6f\n') % (
        a, ecc, incl*rad, node*rad, argp*rad, nu*rad,
        m*rad, year, mon, day, hr, minute, sec)


def load_tests(loader, tests, ignore):
    """Run our main documentation as a test."""

    # Python 2.6 formats floating-point numbers a bit differently and
    # breaks the doctest.
    if sys.version_info >= (2, 7):
        tests.addTests(DocTestSuite('sgp4', optionflags=ELLIPSIS))

    return tests
