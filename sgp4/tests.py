"""Test suite for SGP4."""

try:
    from unittest2 import TestCase, main
except:
    from unittest import TestCase, main

import datetime as dt
import re
import sys
from doctest import DocTestSuite, ELLIPSIS
from math import pi, isnan
from pkgutil import get_data

from sgp4.api import WGS72OLD, WGS72, WGS84, Satrec, jday, accelerated
from sgp4.earth_gravity import wgs72
from sgp4.ext import invjday, newtonnu, rv2coe
from sgp4.functions import _day_of_year_to_month_day
from sgp4.propagation import sgp4, sgp4init
from sgp4 import conveniences, io
from sgp4.exporter import export_tle
import sgp4.model as model

_testcase = TestCase('setUp')
assertEqual = _testcase.assertEqual
assertAlmostEqual = _testcase.assertAlmostEqual
assertRaises = _testcase.assertRaises
assertRaisesRegex = getattr(_testcase, 'assertRaisesRegex',
                            _testcase.assertRaisesRegexp)

error = 2e-7
rad = 180.0 / pi
LINE1 = '1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753'
LINE2 = '2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667'
BAD2  = '2 00007  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413669'
VANGUARD_ATTRS = {
    # Identity
    'satnum': 5,
    'operationmode': 'i',
    # Time
    'jdsatepoch': 2451722.5,
    # Orbit
    'bstar': 2.8098e-05,
    'ndot': 6.96919666594958e-13,
    'nddot': 0.0,
    'ecco': 0.1859667,
    'argpo': 5.790416027488515,
    'inclo': 0.5980929187319208,
    'mo': 0.3373093125574321,
    'no_kozai': 0.04722944544077857,
    'nodeo': 6.08638547138321,
}
VANGUARD_EPOCH = 18441.7849506199999894
VANGUARD_EPOCHDAYS = 179.78495062

# Handle deprecated assertRaisesRegexp, but allow its use Python 2.6 and 2.7
if sys.version_info[:2] == (2, 7) or sys.version_info[:2] == (2, 6):
    TestCase.assertRaisesRegex = TestCase.assertRaisesRegexp

# ------------------------------------------------------------------------
#                           Core Attributes
#

def test_satrec_built_with_twoline2rv():
    sat = Satrec.twoline2rv(LINE1, LINE2)
    verify_vanguard_1(sat)
    assertEqual(sat.epochdays, VANGUARD_EPOCHDAYS)

def test_legacy_built_with_twoline2rv():
    sat = io.twoline2rv(LINE1, LINE2, wgs72)
    verify_vanguard_1(sat, legacy=True)
    assertEqual(sat.epochdays, VANGUARD_EPOCHDAYS)

def test_satrec_initialized_with_sgp4init():
    # epochyr and epochdays are not set by sgp4init
    sat = Satrec()
    sat.sgp4init(
        WGS72,
        'i',
        VANGUARD_ATTRS['satnum'],
        VANGUARD_EPOCH,
        *sgp4init_args(VANGUARD_ATTRS)
    )
    verify_vanguard_1(sat)

def test_satrec_initialized_with_sgp4init_in_afspc_mode():
    sat = Satrec()
    sat.sgp4init(
        WGS72,
        'a',
        VANGUARD_ATTRS['satnum'],
        VANGUARD_EPOCH,
        *sgp4init_args(VANGUARD_ATTRS)
    )
    assertEqual(sat.operationmode, 'a')

def test_legacy_initialized_with_sgp4init():
    sat = model.Satellite()
    sgp4init(
        wgs72, 'i', VANGUARD_ATTRS['satnum'], VANGUARD_EPOCH,
        *sgp4init_args(VANGUARD_ATTRS) + (sat,)
    )
    verify_vanguard_1(sat, legacy=True)

# ------------------------------------------------------------------------
#                 Other Officially Supported Routines
#

def test_jday2():
    jd, fr = jday(2019, 10, 9, 16, 57, 15)
    assertEqual(jd, 2458765.5)
    assertAlmostEqual(fr, 0.7064236111111111)

def test_jday_datetime():
    # define local time
    # UTC equivalent: 2011-11-03 20:05:23+00:00

    class UTC_plus_4(dt.tzinfo):
        'UTC'
        offset = dt.timedelta(hours=4)
        def utcoffset(self, datetime):
            return self.offset
        def tzname(self, datetime):
            return 'UTC plus 4'
        def dst(self, datetime):
            return self.offset

    datetime_local = dt.datetime(2011, 11, 4, 0, 5, 23, 0, UTC_plus_4())
    jd, fr = conveniences.jday_datetime(datetime_local)

    # jd of this date is 2455868.5 + 0.8370717592592593
    assertEqual(jd, 2455868.5)
    assertAlmostEqual(fr, 0.8370717592592593)

def test_sat_epoch_datetime():
    sat = Satrec.twoline2rv(LINE1, LINE2)
    datetime = conveniences.sat_epoch_datetime(sat)
    zone = conveniences.UTC
    assertEqual(datetime, dt.datetime(2000, 6, 27, 18, 50, 19, 733568, zone))

def test_good_tle_checksum():
    for line in LINE1, LINE2:
        checksum = int(line[-1])
        assertEqual(io.compute_checksum(line), checksum)
        assertEqual(io.fix_checksum(line[:68]), line)
        io.verify_checksum(line)

def test_bad_tle_checksum():
    checksum = LINE1[-1]
    assertEqual(checksum, '3')
    bad = LINE1[:68] + '7'
    assertRaises(ValueError, io.verify_checksum, bad)
    assertEqual(io.fix_checksum(bad), LINE1)

def test_tle_export():
    """Check `export_tle()` round-trip using all the TLEs in the test file.

    This iterates through the satellites in "SGP4-VER.TLE",
    generates `Satrec` objects and exports the TLEs.  These exported
    TLEs are then compared to the original TLE, closing the loop (or
    the round-trip).

    """
    data = get_data(__name__, 'SGP4-VER.TLE')
    tle_lines = iter(data.decode('ascii').splitlines())

    # Skip these lines, known errors
    # Resulting TLEs are equivalent (same values in the Satrec object), but they are not the same
    # 25954: BSTAR = 0 results in a negative exp, not positive
    # 29141: BSTAR = 0.13519 results in a negative exp, not positive
    # 33333: Checksum error as expected on both lines
    # 33334: Checksum error as expected on line 1
    # 33335: Checksum error as expected on line 1
    expected_errs_line1 = set([25954, 29141, 33333, 33334, 33335])
    expected_errs_line2 = set([33333, 33335])

    if accelerated:
        # Non-standard: omits the ephemeris type integer.
        expected_errs_line1.add(11801)

    for line1 in tle_lines:

        if not line1.startswith('1'):
            continue

        line2 = next(tle_lines)

        # trim lines to normal TLE string size
        line1 = line1[:69]
        line2 = line2[:69]
        satrec = Satrec.twoline2rv(line1, line2)

        # Generate TLE from satrec
        out_line1, out_line2 = export_tle(satrec)

        if satrec.satnum not in expected_errs_line1:
            assertEqual(out_line1, line1)
        if satrec.satnum not in expected_errs_line2:
            assertEqual(out_line2, line2)

def test_all_three_gravity_models_with_twoline2rv():
    # The numbers below are those produced by Vallado's C++ code.
    # (Why does the Python version not produce the same values to
    # high accuracy, instead of agreeing to only 4 places?)

    assert_wgs72old(Satrec.twoline2rv(LINE1, LINE2, WGS72OLD))
    assert_wgs72(Satrec.twoline2rv(LINE1, LINE2, WGS72))
    assert_wgs84(Satrec.twoline2rv(LINE1, LINE2, WGS84))

    # Not specifying a gravity model should select WGS72.

    assert_wgs72(Satrec.twoline2rv(LINE1, LINE2))

def test_all_three_gravity_models_with_sgp4init():
    # Gravity models specified with sgp4init() should also change the
    # positions generated.

    sat = Satrec()
    args = sgp4init_args(VANGUARD_ATTRS)

    sat.sgp4init(WGS72OLD, 'i', VANGUARD_ATTRS['satnum'], VANGUARD_EPOCH, *args)
    assert_wgs72old(sat)

    sat.sgp4init(WGS72, 'i', VANGUARD_ATTRS['satnum'], VANGUARD_EPOCH, *args)
    assert_wgs72(sat)

    sat.sgp4init(WGS84, 'i', VANGUARD_ATTRS['satnum'], VANGUARD_EPOCH, *args)
    assert_wgs84(sat)

GRAVITY_DIGITS = 12 if accelerated else 4

def assert_wgs72old(sat):
    e, r, v = sat.sgp4_tsince(309.67110720001529)
    assertAlmostEqual(r[0], -3754.251473242793, GRAVITY_DIGITS)
    assertAlmostEqual(r[1], 7876.346815095482, GRAVITY_DIGITS)
    assertAlmostEqual(r[2], 4719.220855042922, GRAVITY_DIGITS)

def assert_wgs72(sat):
    e, r, v = sat.sgp4_tsince(309.67110720001529)
    assertAlmostEqual(r[0], -3754.2514743216166, GRAVITY_DIGITS)
    assertAlmostEqual(r[1], 7876.346817439062, GRAVITY_DIGITS)
    assertAlmostEqual(r[2], 4719.220856478582, GRAVITY_DIGITS)

def assert_wgs84(sat):
    e, r, v = sat.sgp4_tsince(309.67110720001529)
    assertAlmostEqual(r[0], -3754.2437675772426, GRAVITY_DIGITS)
    assertAlmostEqual(r[1], 7876.3549956188945, GRAVITY_DIGITS)
    assertAlmostEqual(r[2], 4719.227897029576, GRAVITY_DIGITS)

# ------------------------------------------------------------------------
#                            Special Cases
#

def test_intldesg_with_6_characters():
    sat = Satrec.twoline2rv(LINE1, LINE2)
    assertEqual(sat.intldesg, '58002B')

def test_intldesg_with_7_characters():
    sat = Satrec.twoline2rv(
        '1 39444U 13066AE  20110.89708219  .00000236  00000-0'
        '  35029-4 0  9992',
        '2 39444  97.5597 114.3769 0059573 102.0933 258.6965 '
        '14.82098949344697',
    )
    assertEqual(sat.intldesg, '13066AE')

def test_hyperbolic_orbit():
    # Exercise the newtonnu() code path with asinh() to see whether
    # we can replace it with the one from Python's math module.

    e0, m = newtonnu(1.0, 2.9)  # parabolic
    assertAlmostEqual(e0, 8.238092752965605, places=12)
    assertAlmostEqual(m, 194.60069989482898, places=12)

    e0, m = newtonnu(1.1, 2.7)   # hyperbolic
    assertAlmostEqual(e0, 4.262200676156417, places=12)
    assertAlmostEqual(m, 34.76134082028372, places=12)

def test_correct_epochyr():
    # Make sure that the non-standard four-digit epochyr I switched
    # to in the Python version of SGP4 is reverted back to the
    # official behavior when that code is used behind Satrec.
    sat = Satrec.twoline2rv(LINE1, LINE2)
    assertEqual(sat.epochyr, 0)

def test_legacy_epochyr():
    # Apparently I saw fit to change the meaning of this attribute
    # in the Python version of SGP4.
    sat = io.twoline2rv(LINE1, LINE2, wgs72)
    assertEqual(sat.epochyr, 2000)

def test_support_for_old_no_attribute():
    s = io.twoline2rv(LINE1, LINE2, wgs72)
    assert s.no == s.no_kozai

def test_months_and_days():
    # Make sure our hand-written months-and-days routine is perfect.

    month_lengths = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    day_of_year = 1
    for month, length in enumerate(month_lengths, 1):
        for day in range(1, length + 1):
            tup = _day_of_year_to_month_day(day_of_year, False)
            assertEqual((month, day), tup)
            day_of_year += 1

    month_lengths[1] = 29  # February, during a leap year
    day_of_year = 1
    for month, length in enumerate(month_lengths, 1):
        for day in range(1, length + 1):
            tup = _day_of_year_to_month_day(day_of_year, True)
            assertEqual((month, day), tup)
            day_of_year += 1


def test_december_32():
    # ISS [Orbit 606], whose date is 2019 plus 366.82137887 days.
    # The core SGP4 routines handled this fine, but my hamfisted
    # attempt to provide a Python datetime for "convenience" ran
    # into an overflow.
    a = '1 25544U 98067A   19366.82137887  .00016717  00000-0  10270-3 0  9129'
    b = '2 25544  51.6392  96.6358 0005156  88.7140 271.4601 15.49497216  6061'
    correct_epoch = dt.datetime(2020, 1, 1, 19, 42, 47, 134368)

    # Legacy API.
    sat = io.twoline2rv(a, b, wgs72)
    assertEqual(sat.epoch, correct_epoch)

    correct_epoch = correct_epoch.replace(tzinfo=conveniences.UTC)

    # Modern API.
    sat = Satrec.twoline2rv(a, b)
    assertEqual(conveniences.sat_epoch_datetime(sat), correct_epoch)


def test_bad_first_line():
    with assertRaisesRegex(ValueError, re.escape("""TLE format error

The Two-Line Element (TLE) format was designed for punch cards, and so
is very strict about the position of every period, space, and digit.
Your line does not quite match.  Here is the official format for line 1
with an N where each digit should go, followed by the line you provided:

1 NNNNNC NNNNNAAA NNNNN.NNNNNNNN +.NNNNNNNN +NNNNN-N +NNNNN-N N NNNNN
1 00005U 58002B   00179.78495062  .000000234 00000-0  28098-4 0  4753""")):
        io.twoline2rv(LINE1.replace('23 ', '234'), LINE2, wgs72)

def test_bad_second_line():
    with assertRaisesRegex(ValueError, re.escape("""TLE format error

The Two-Line Element (TLE) format was designed for punch cards, and so
is very strict about the position of every period, space, and digit.
Your line does not quite match.  Here is the official format for line 2
with an N where each digit should go, followed by the line you provided:

2 NNNNN NNN.NNNN NNN.NNNN NNNNNNN NNN.NNNN NNN.NNNN NN.NNNNNNNNNNNNNN
2 00005 34 .268234 8.7242 1859667 331.7664  19.3264 10.82419157413667""")):
        io.twoline2rv(LINE1, LINE2.replace(' 34', '34 '), wgs72)

def test_mismatched_lines():
    msg = "Object numbers in lines 1 and 2 do not match"
    with assertRaisesRegex(ValueError, re.escape(msg)):
        io.twoline2rv(LINE1, BAD2, wgs72)

# ------------------------------------------------------------------------
#                           Helper routines
#

def verify_vanguard_1(sat, legacy=False):
    attrs = VANGUARD_ATTRS

    if legacy:
        attrs = attrs.copy()
        del attrs['jdsatepoch']

    for name, value in attrs.items():
        try:
            assertEqual(getattr(sat, name), value)
        except AssertionError as e:
            message, = e.args
            e.args = ('for attribute %s, %s' % (name, message),)
            raise e

    if not legacy:
        assertAlmostEqual(sat.jdsatepochF, 0.78495062, delta=1e-13)

def sgp4init_args(d):
    """Given a dict of orbital parameters, return them in sgp4init order."""
    return (d['bstar'], d['ndot'], d['nddot'], d['ecco'], d['argpo'],
            d['inclo'], d['mo'], d['no_kozai'], d['nodeo'])

# ----------------------------------------------------------------------
#                           INTEGRATION TEST
#
# This runs both new and old satellite objects against every example
# computation in the official `tcppver.out` test case file.  Instead of
# trying to parse the file, it instead re-generates it using Python
# satellite objects, then compares the resulting text with the file.

def test_satrec_against_tcppver_using_julian_dates():

    def invoke(satrec, tsince):
        whole, fraction = divmod(tsince / 1440.0, 1.0)
        jd = satrec.jdsatepoch + whole
        fr = satrec.jdsatepochF + fraction
        e, r, v = satrec.sgp4(jd, fr)
        return e, r, v

    run_satellite_against_tcppver(Satrec.twoline2rv, invoke, [1,1,6,6,4,3,6])

def test_satrec_against_tcppver_using_tsince():

    def invoke(satrec, tsince):
        e, r, v = satrec.sgp4_tsince(tsince)
        return e, r, v

    run_satellite_against_tcppver(Satrec.twoline2rv, invoke, [1,1,6,6,4,3,6])

def test_legacy_against_tcppver():

    def make_legacy_satellite(line1, line2):
        sat = io.twoline2rv(line1, line2, wgs72)
        return sat

    def run_legacy_sgp4(satrec, tsince):
        r, v = sgp4(satrec, tsince)
        return (satrec.error, satrec.error_message), r, v

    errs = [
        (1, 'mean eccentricity -0.001329 not within range 0.0 <= e < 1.0'),
        (1, 'mean eccentricity -0.001208 not within range 0.0 <= e < 1.0'),
        (6, 'mrt 0.996159 is less than 1.0'
         ' indicating the satellite has decayed'),
        (6, 'mrt 0.996252 is less than 1.0'
         ' indicating the satellite has decayed'),
        (4, 'semilatus rectum -0.103223 is less than zero'),
        (3, 'perturbed eccentricity -122.217193'
         ' not within range 0.0 <= e <= 1.0'),
        (6, 'mrt 0.830534 is less than 1.0'
         ' indicating the satellite has decayed'),
    ]

    run_satellite_against_tcppver(make_legacy_satellite, run_legacy_sgp4, errs)

def run_satellite_against_tcppver(twoline2rv, invoke, expected_errors):
    # Check whether a test run produces the output in tcppver.out

    error_list = []
    actual = generate_test_output(twoline2rv, invoke, error_list)
    previous_data_line = None

    # Iterate across "tcppver.out", making sure that we ourselves
    # produce a line that looks very much like the corresponding
    # line in that file.

    data = get_data(__name__, 'tcppver.out')
    tcppver_lines = data.decode('ascii').splitlines(True)
    for i, expected_line in enumerate(tcppver_lines, start = 1):

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
                'Expected: %r\n'
                'Got back: %r'
                % (i, expected_line, actual_line))

        if 'xx' not in actual_line:
            previous_data_line = actual_line

    # Make sure the test file is not missing lines.

    missing_count = 0
    for actual_line in actual:
        missing_count += 1

    if missing_count > 0:
        raise ValueError('we produced %d extra lines' % (missing_count,))

    assertEqual(error_list, expected_errors)

def generate_test_output(twoline2rv, invoke, error_list):
    """Generate lines like those in the test file tcppver.out.

    This iterates through the satellites in "SGP4-VER.TLE", which are
    each supplemented with a time start/stop/step over which we are
    supposed to print results.

    """
    data = get_data(__name__, 'SGP4-VER.TLE')
    tle_lines = iter(data.decode('ascii').splitlines())

    for line1 in tle_lines:

        if not line1.startswith('1'):
            continue

        line2 = next(tle_lines)
        satrec = twoline2rv(line1, line2)

        yield '%ld xx\n' % (satrec.satnum,)

        for line in generate_satellite_output(
                satrec, invoke, line2, error_list):
            yield line

def generate_satellite_output(satrec, invoke, line2, error_list):
    """Print a data line for each time in line2's start/stop/step field."""

    mu = wgs72.mu

    e, r, v = invoke(satrec, 0.0)
    if isnan(r[0]) and isnan(r[1]) and isnan(r[2]):
        error_list.append(e)
        yield '(Use previous data line)'
        return
    yield format_short_line(0.0, r, v)

    tstart, tend, tstep = (float(field) for field in line2[69:].split())

    tsince = tstart
    while tsince <= tend:
        if tsince == tstart == 0.0:
            tsince += tstep
            continue  # avoid duplicating the first line

        e, r, v = invoke(satrec, tsince)

        if e != 0 and e != (0, None):
            error_list.append(e)
            return
        yield format_long_line(satrec, tsince, mu, r, v)

        tsince += tstep

    if tsince - tend < tstep - 1e-6:  # do not miss last line!
        e, r, v = invoke(satrec, tend)
        if e != 0 and e != (0, None):
            error_list.append(e)
            return
        yield format_long_line(satrec, tend, mu, r, v)

def format_short_line(tsince, r, v):
    """Short line, using the same format string that testcpp.cpp uses."""

    return ' %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f\n' % (
        tsince, r[0], r[1], r[2], v[0], v[1], v[2])

def format_long_line(satrec, tsince, mu, r, v):
    """Long line, using the same format string that testcpp.cpp uses."""

    short = format_short_line(tsince, r, v).strip('\n')

    jd = satrec.jdsatepoch + satrec.jdsatepochF + tsince / 1440.0
    year, mon, day, hr, minute, sec = invjday(jd)

    (p, a, ecc, incl, node, argp, nu, m, arglat, truelon, lonper
     ) = rv2coe(r, v, mu)

    return short + (
        ' %14.6f %8.6f %10.5f %10.5f %10.5f %10.5f %10.5f'
        ' %5i%3i%3i %2i:%2i:%9.6f\n'
    ) % (
        a, ecc, incl*rad, node*rad, argp*rad, nu*rad,
        m*rad, year, mon, day, hr, minute, sec,
    )

# ----------------------------------------------------------------------

def load_tests(loader, tests, ignore):
    """Run our main documentation as a test, plus all test functions."""

    from sgp4.wulfgar import add_test_functions
    add_test_functions(loader, tests, __name__)

    # Python 2.6 formats floating-point numbers a bit differently and
    # breaks the doctest, so we only run the doctest on later versions.
    if sys.version_info >= (2, 7):
        tests.addTests(DocTestSuite('sgp4', optionflags=ELLIPSIS))
        tests.addTests(DocTestSuite('sgp4.conveniences', optionflags=ELLIPSIS))
        tests.addTests(DocTestSuite('sgp4.functions', optionflags=ELLIPSIS))

    return tests

if __name__ == '__main__':
    main()
