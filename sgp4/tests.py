"""Test suite for SPG4."""

import os
from unittest import TestCase
from sgp4.io import twoline2rv
from sgp4.propagation import sgp4

thisdir = os.path.dirname(__file__)
error = 2e-7


class SatRec(object):
    pass


class Tests(TestCase):

    def test_tle_verify(self):
        whichconst = 'wgs72'
        actual = generate_test_output(whichconst)
        previous_data_line = None

        tcppath = os.path.join(thisdir, 'tcppver.out')
        with open(tcppath) as tcpfile:
            for i, expected_line in enumerate(tcpfile, start = 1):

                # TODO: what are we supposed to do with the extra fields
                # that lie past character 107?
                if len(expected_line) > 107:
                    expected_line = expected_line[:107] + expected_line[-1:]

                try:
                    actual_line = next(actual)
                except StopIteration:
                    print 'WARNING: our output ended early, on line %d' % (i,)
                    break

                if actual_line == '(Use previous data line)':
                    actual_line = '       0.00000000' + previous_data_line[17:]

                # Compare the lines.

                if 'xx' in actual_line:
                    similar = (actual_line == expected_line)
                else:
                    actuals = (float(a) for a in actual_line.split())
                    expecteds = (float(e) for e in expected_line.split())
                    similar = all(-error < (a - e) < error
                                   for a, e in zip(actuals, expecteds))

                if not similar:
                    print(
                        'Line %d of output does not match:\n'
                        '\n'
                        'Expect: %s'
                        'Actual: %s'
                        % (i, expected_line, actual_line))

                if 'xx' not in actual_line:
                    previous_data_line = actual_line

        missing_count = 0
        for actual_line in actual:
            missing_count += 1

        if missing_count > 0:
            raise ValueError('we produced %d extra lines' % (missing_count,))


def generate_test_output(whichconst):
    """Generate lines like those in the test file tcppver.out."""

    tlepath = os.path.join(thisdir, 'SGP4-VER.TLE')
    with open(tlepath) as tlefile:
        tlelines = iter(tlefile.readlines())

    for line1 in tlelines:

        if not line1.startswith('1'):
            continue

        line2 = next(tlelines)

        satrec = SatRec()
        twoline2rv(line1, line2, 'c', None, 'i', whichconst, satrec)
        yield '%ld xx\n' % (satrec.satnum,)

        for line in generate_satellite_output(whichconst, satrec, line2):
            yield line


def generate_satellite_output(whichconst, satrec, line2):

    ro, vo = sgp4(whichconst, satrec, 0.0)
    if ro is None:
        yield '(Use previous data line)'
        return
    yield format_test_line(satrec, ro, vo)

    tstart, tend, tstep = (float(field) for field in line2[69:].split())

    tsince = tstart
    while tsince <= tend:
        if tsince == tstart == 0.0:
            tsince += tstep
            continue  # avoid duplicating the first line

        ro, vo = sgp4(whichconst, satrec, tsince)
        if ro is None:
            return
        yield format_test_line(satrec, ro, vo)

        tsince += tstep

    if tsince - tend < tstep - 1e-6:  # do not miss last line!
        ro, vo = sgp4(whichconst, satrec, tend)
        if ro is None:
            return
        yield format_test_line(satrec, ro, vo)


def format_test_line(satrec, ro, vo):
    return ' %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f\n' % (
        satrec.t, ro[0], ro[1], ro[2], vo[0], vo[1], vo[2])
