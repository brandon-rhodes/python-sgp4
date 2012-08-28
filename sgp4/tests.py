"""Test suite for SGP4."""

import os
from doctest import DocTestSuite
from unittest import TestCase
from math import pi

from sgp4.earth_gravity import wgs72
from sgp4.ext import invjday, rv2coe
from sgp4.io import twoline2rv
from sgp4.propagation import sgp4

thisdir = os.path.dirname(__file__)
error = 2e-7
rad = 180.0 / pi


class Tests(TestCase):

    def test_tle_verify(self):
        # Check whether a test run produces the output in tcppver.out

        whichconst = 'wgs72'
        actual = generate_test_output(whichconst)
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


def generate_test_output(whichconst):
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

        for line in generate_satellite_output(satrec, line2):
            yield line


def generate_satellite_output(satrec, line2):
    """Print a data line for each time in line2's start/stop/step field."""

    mu = satrec.whichconst.mu

    r, v = sgp4(satrec, 0.0)
    if r is None:
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

        if r is None:
            return
        yield format_long_line(satrec, mu, r, v)

        tsince += tstep

    if tsince - tend < tstep - 1e-6:  # do not miss last line!
        r, v = sgp4(satrec, tend)
        if r is None:
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
    tests.addTests(DocTestSuite('sgp4'))
    return tests
