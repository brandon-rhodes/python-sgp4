# -*- coding: utf-8 -*-
"""Track earth satellite TLE orbits using up-to-date 2020 version of SGP4

This Python package computes the position and velocity of an
earth-orbiting satellite, given the satellite's TLE orbital elements
from a source like `Celestrak <https://celestrak.com/>`_.  It implements
the most recent version of SGP4, and is regularly run against the SGP4
test suite to make sure that its satellite position predictions **agree
to within 0.1 mm** with the predictions of the standard distribution of
the algorithm.  This error is far less than the 1–3 km/day by which
satellites themselves deviate from the ideal orbits described in TLE
files.

* If your platform supports it, this package compiles the verbatim
  source code from the official C++ version of SGP4.  You can call the
  routine directly, or through an array API that loops over arrays of
  satellites and arrays of times with machine code instead of Python.

* Otherwise, a slower but reliable Python implementation of SGP4 is used
  instead.

Note that the SGP4 propagator returns raw *x,y,z* Cartesian coordinates
in a “True Equator Mean Equinox” (TEME) reference frame that’s centered
on the Earth but does not rotate with it — an “Earth centered inertial”
(ECI) reference frame.  The SGP4 propagator itself does not implement
the math to convert these positions into more official ECI frames like
J2000 or the ICRF; nor to convert positions into any Earth-centered
Earth-fixed (ECEF) frames like the ITRS; nor to convert them to
latitudes and longitudes through an Earth ellipsoid like WGS84.

For conversions into other coordinate frames, look for a comprehensive
astronomy library that is built atop this one, like the `Skyfield
<https://rhodesmill.org/skyfield/>`_ library:

https://rhodesmill.org/skyfield/earth-satellites.html

Usage
-----

This library uses the same function names as the official C++ code, to
help users who may already be familiar with SGP4 in other languages.
Here is how to compute the x,y,z position and velocity for the
International Space Station at 12:50:19 on 29 June 2000:

>>> from sgp4.api import Satrec
>>>
>>> s = '1 25544U 98067A   19343.69339541  .00001764  00000-0  38792-4 0  9991'
>>> t = '2 25544  51.6439 211.2001 0007417  17.6667  85.6398 15.50103472202482'
>>> satellite = Satrec.twoline2rv(s, t)
>>>
>>> jd, fr = 2458827, 0.362605
>>> e, r, v = satellite.sgp4(jd, fr)
>>> e
0
>>> print(r)  # True Equator Mean Equinox position (km)
(-6102.44..., -986.33..., -2820.31...)
>>> print(v)  # True Equator Mean Equinox velocity (km/s)
(-1.45..., -5.52..., 5.10...)

As input, you can provide either:

* A simple floating-point Julian Date for ``jd`` and the value 0.0 for
  ``fr``, if you are happy with the precision of a 64-bit floating point
  number.  Note that modern Julian Dates are greater than 2,450,000
  which means that nearly half of the precision of a 64-bit float will
  be consumed by the whole part that specifies the day.  The remaining
  digits will provide a precision for the fraction of around 20.1 µs.
  This should be no problem for the accuracy of your result — satellite
  positions usually off by a few kilometers anyway, far less than a
  satellite moves in 20.1 µs — but if you run a solver that dives down
  into the microseconds while searching for a rising or setting time,
  the solver might be bothered by the 20.1 µs plateau between each jump
  in the satellite’s position.

* Or, you can provide a coarse date ``jd`` plus a very precise fraction
  ``fr`` that supplies the rest of the value.  The Julian Date for which
  the satellite position is computed is the sum of the two values.  One
  common practice is to provide the whole number as ``jd`` and the
  fraction as ``fr``; another is to have ``jd`` carry the fraction 0.5
  since UTC midnight occurs halfway through each Julian Date.  Either
  way, splitting the value allows a solver to run all the way down into
  the nanoseconds and still see SGP4 respond smoothly to tiny date
  adjustments with tiny changes in the resulting satellite position.

Here is how to intrepret the results:

* ``e`` will be a non-zero error code if the satellite position could
  not be computed for the given date.  You can ``from sgp4.api import
  SGP4_ERRORS`` to access a dictionary mapping error codes to error
  messages explaining what each code means.

* ``r`` measures the satellite position in **kilometers** from the
  center of the earth in the idiosyncratic True Equator Mean Equinox
  coordinate frame used by SGP4.

* ``v`` velocity is the rate at which the position is changing,
  expressed in **kilometers per second**.

If your application does not natively handle Julian dates, you can
compute ``jd`` and ``fr`` from calendar dates using ``jday()``.

>>> from sgp4.api import jday
>>> jd, fr = jday(2019, 12, 9, 12, 0, 0)
>>> jd
2458826.5
>>> fr
0.5

OMM
---

The industry is making adjustments because the fixed-width TLE format
will soon run out of satellite numbers.

* Some TLE files now use a new “Alpha-5” convention that expands the
  range of satellite numbers by using an initial letter; for example,
  “E8493” means satellite 148493.  This library now supports the Alpha-5
  convention and should return the correct integer in Python.

* Some authorities are now distributing satellite elements in an “OMM”
  Orbit Mean Elements Message format that replaces the TLE format.

You can already try out experimental support for OMM:

>>> from sgp4 import omm

Reading OMM data takes two steps, because OMM supports several different
text formats.  First, parse the input text to recover the field names
and values that it stores; second, build a Python satellite object from
those field values.  For example, to load OMM from XML:

>>> with open('sample_omm.xml') as f:
...     fields = next(omm.parse_xml(f))
>>> sat = Satrec()
>>> omm.initialize(sat, fields)

Or, to load OMM from CSV:

>>> with open('sample_omm.csv') as f:
...     fields = next(omm.parse_csv(f))
>>> sat = Satrec()
>>> omm.initialize(sat, fields)

Either way, the satellite object should wind up properly initialized and
ready to start producing positions.

Epoch
-----

Over a given satellite’s lifetime, dozens or hundreds of different TLE
records will be produced as its orbit evolves.  Each TLE record
specifies the “epoch date” for which it is most accurate.  Typically a
TLE is only useful for a couple of weeks to either side of its epoch
date, beyond which its predictions become unreliable.

Satellite objects natively provide their epoch as a two-digit year and
then a fractional number of days into the year:

>>> satellite.epochyr
19
>>> satellite.epochdays
343.69339541

Because Sputnik was launched in 1957, satellite element sets will never
refer to an earlier year, so years 57 through 99 mean 1957–1999 while 0
through 56 mean 2000–2056.  The TLE format will presumably be obsolete
in 2057 and have to be upgraded to 4-digit years.

To turn the number of days and its fraction into a calendar date and
time, use the ``days2mdhms()`` function.

>>> from sgp4.api import days2mdhms
>>> month, day, hour, minute, second = days2mdhms(19, 343.69339541)
>>> month
12
>>> day
9
>>> hour
16
>>> minute
38
>>> second
29.363424

The SGP4 library also translates those two numbers into a Julian date
and fractional Julian date, since Julian dates are more commonly used in
astronomy.

>>> satellite.jdsatepoch
2458826.5
>>> satellite.jdsatepochF
0.69339541

Finally, a convenience function is available in the library if you need
the epoch date and time as Python ``datetime``.

>>> from sgp4.conveniences import sat_epoch_datetime
>>> sat_epoch_datetime(satellite)
datetime.datetime(2019, 12, 9, 16, 38, 29, 363423, tzinfo=UTC)

Array Acceleration
------------------

To avoid the expense of Python loops when you have many dates, you can
pass them as arrays to another method that understands NumPy:

>>> import numpy as np
>>> np.set_printoptions(precision=2)

>>> jd = np.array((2458826, 2458826, 2458826, 2458826))
>>> fr = np.array((0.0001, 0.0002, 0.0003, 0.0004))

>>> e, r, v = satellite.sgp4_array(jd, fr)

>>> print(e)
[0 0 0 0]
>>> print(r)
[[-3431.31  2620.15 -5252.97]
 [-3478.86  2575.14 -5243.87]
 [-3526.09  2529.89 -5234.28]
 [-3572.98  2484.41 -5224.19]]
>>> print(v)
[[-5.52 -5.19  1.02]
 [-5.49 -5.22  1.08]
 [-5.45 -5.25  1.14]
 [-5.41 -5.28  1.2 ]]

To avoid the expense of Python loops when you have many satellites and
dates, build a ``SatrecArray`` from several individual satellites.  Its
``sgp4()`` method will expect both ``jd`` and ``fr`` to be NumPy arrays,
so if you only have one date, be sure to provide NumPy arrays of length
one.  Here is a sample computation for 2 satellites and 4 dates:

>>> s = '1 20580U 90037B   19342.88042116  .00000361  00000-0  11007-4 0  9996'
>>> t = '2 20580  28.4682 146.6676 0002639 185.9222 322.7238 15.09309432427086'
>>> satellite2 = Satrec.twoline2rv(s, t)

>>> from sgp4.api import SatrecArray
>>> a = SatrecArray([satellite, satellite2])
>>> e, r, v = a.sgp4(jd, fr)

>>> np.set_printoptions(precision=2)
>>> print(e)
[[0 0 0 0]
 [0 0 0 0]]
>>> print(r)
[[[-3431.31  2620.15 -5252.97]
  [-3478.86  2575.14 -5243.87]
  [-3526.09  2529.89 -5234.28]
  [-3572.98  2484.41 -5224.19]]
<BLANKLINE>
 [[ 5781.85  2564.   -2798.22]
  [ 5749.36  2618.59 -2814.63]
  [ 5716.35  2672.94 -2830.78]
  [ 5682.83  2727.05 -2846.68]]]
>>> print(v)
[[[-5.52 -5.19  1.02]
  [-5.49 -5.22  1.08]
  [-5.45 -5.25  1.14]
  [-5.41 -5.28  1.2 ]]
<BLANKLINE>
 [[-3.73  6.33 -1.91]
  [-3.79  6.3  -1.88]
  [-3.85  6.28 -1.85]
  [-3.91  6.25 -1.83]]]

Attributes
----------

The attributes of a ``Satrec`` object carry the data loaded from the TLE
entry.  Look at the class's documentation for details.

Export
------

If you have a ``Satrec`` you want to share with friends or persist to a
file, there’s an export routine that will turn it back into a TLE:

>>> from sgp4.exporter import export_tle
>>> line1, line2 = export_tle(satellite)
>>> line1
'1 25544U 98067A   19343.69339541  .00001764  00000-0  38792-4 0  9991'
>>> line2
'2 25544  51.6439 211.2001 0007417  17.6667  85.6398 15.50103472202482'

Gravity
-------

The SGP4 algorithm operates atop a set of constants specifying how
strong the Earth’s gravity is.  The most recent official paper on SGP4
(see below) specifies that “We use WGS-72 as the default value”, so this
Python module uses the same default.  But in case you want to use either
the old legacy version of the WGS-72 constants, or else the non-standard
but more modern WGS-84 constants, the ``twoline2rv()`` constructor takes
an optional argument:

>>> from sgp4.api import WGS72OLD, WGS72, WGS84
>>> satellite3 = Satrec.twoline2rv(s, t, WGS84)

You will in general get less accurate results if you choose WGS-84.
Even though it reflects more recent and accurate measures of the Earth,
satellite TLEs across the industry are most likely generated with WGS-72
as their basis.  The positions you generate will better agree with the
real positions of each satellite if you use the same underlying gravity
constants as were used to generate the TLE.

Providing your own elements
---------------------------

If instead of parsing a TLE you want to provide your own orbital
elements, you can call the ``sgp4init()`` method of any existing
satellite object to reset it to those new elements.

>>> sat = Satrec()
>>> sat.sgp4init(
...     WGS72,           # gravity model
...     'i',             # 'a' = old AFSPC mode, 'i' = improved mode
...     5,               # satnum: Satellite number
...     18441.785,       # epoch: days since 1949 December 31 00:00 UT
...     2.8098e-05,      # bstar: drag coefficient (/earth radii)
...     6.969196665e-13, # ndot: ballistic coefficient (revs/day)
...     0.0,             # nddot: second derivative of mean motion (revs/day^3)
...     0.1859667,       # ecco: eccentricity
...     5.7904160274885, # argpo: argument of perigee (radians)
...     0.5980929187319, # inclo: inclination (radians)
...     0.3373093125574, # mo: mean anomaly (radians)
...     0.0472294454407, # no_kozai: mean motion (radians/minute)
...     6.0863854713832, # nodeo: right ascension of ascending node (radians)
... )

To compute the “epoch” value, simply take a normal Julian date and
subtract ``2433281.5`` days.

In addition to setting the attributes natively set by the underlying
``sgp4init()`` routine, this library also goes ahead and sets the date
fields ``epochyr``, ``epochdays``, ``jdsatepoch``, and ``jdsatepochF``.

The character provided as the second argument can be ``'a'`` to run the
computations so that they are compatible with the old Air Force Space
Command edition of the library, or ``'i'`` to run the new and improved
version of the SGP4 algorithm.

You can also directly access a satellite’s orbital parameters by asking
for the attributes ``sat.epoch``, ``sat.bstar``, and so forth, using the
names given in the comments above.

Validation against the official algorithm
-----------------------------------------

This implementation passes all of the automated tests in the August 2010
release of the reference implementation of SGP4 by Vallado et al., who
originally published their revision of SGP4 in 2006:

    Vallado, David A., Paul Crawford, Richard Hujsak, and T.S. Kelso,
    “Revisiting Spacetrack Report #3,” presented at the AIAA/AAS
    Astrodynamics Specialist Conference, Keystone, CO, 2006 August
    21–24.

If you would like to review the paper, it is `available online
<https://www.celestrak.com/publications/AIAA/2006-6753/>`_.  You can
always download the latest version of their code for comparison against
this Python module (or other implementations) at `AIAA-2006-6753.zip
<https://www.celestrak.com/publications/AIAA/2006-6753/AIAA-2006-6753.zip>`_.

For developers
--------------

Developers can check out this full project from GitHub:

https://github.com/brandon-rhodes/python-sgp4

To run its unit tests, install Python 2, Python 3, and the ``tox``
testing tool.  The tests runing in Python 2 will exercise the fallback
pure-Python version of the routines, while Python 3 exercises the fast
new C++ accelerated code::

    cd python-sgp4
    tox

Legacy API
----------

Before this library pivoted to wrapping Vallado's official C++ code and
was operating in pure Python only, it had a slightly quirkier API, which
is still supported for compatibility with older clients.  You can learn
about it by reading the documentation from version 1.4 or earlier:

https://pypi.org/project/sgp4/1.4/

Changelog
---------

| 2020-12-16 — 2.14 — New data formats: added OMM message support for both XML and CSV, and added support for the new Alpha-5 extension to TLE files.
| 2020-10-14 — 2.13 — Enhanced ``sgp4init()`` with custom code that also sets the ``epochdays`` and ``epochyr`` satellite attributes.
| 2020-05-28 — 2.12 — Moved the decision of whether to set the locale during ``twoline2rv()`` from import time to runtime, for users who change locales after their application is up and running.
| 2020-05-24 — 2.11 — Fixed a regression in how dates are split into hours, minutes, and seconds that would sometimes produce a time whose second=60, crashing the pure-Python version of the library.
| 2020-05-22 — 2.10 — Switch the locale temporarily to ``C`` during the C++ accelerated ``twoline2rv()``, since it does not protect its ``sscanf()`` calls from locales that, like German, expect comma decimal points instead of the period decimal points always used in a TLE.
| 2020-05-21 — 2.9 — Added ``sat_epoch_datetime()``, expanded documentation around converting a satellite epoch to a date and time, and started rounding the epoch to exactly the digits provided in the TLE; and removed the ``Satrec.epoch`` attribute from Python fallback code to better match the C++ version.
| 2020-05-07 — 2.8 — New function ``jday_datetime()`` is now available in the ``sgp4.conveniences`` module, thanks to Egemen Imre.
| 2020-04-24 — 2.7 — New method ``sgp4init()`` (thank you, Chris Lewicki!) is available.
| 2020-04-20 — 2.6 — New routine ``export_tle()`` (thank you, Egemen Imre!) is available. Improved how the accelerated C++ backend parses the ``intldesg`` string and the ``revnum`` integer.
| 2020-03-22 — 2.5 — Gave the new accelerated ``twoline2rv()`` an optional argument that lets the user choose a non-standard set of gravity constants.
| 2020-02-25 — 2.4 — Improved the ``jday()`` docstring; made the old legacy Python resilient if the day of the month is out-of-range (past the end of the month) in a TLE; and Mark Rutten fixed the C++ so it compiles on Windows!
| 2020-02-04 — 2.3 — Removed experimental code that caused performance problems for users with Numba installed.
| 2020-02-02 — 2.2 — A second release on Palindrome Day: fix the Satrec ``.epochyr`` attribute so it behaves the same way in Python as it does in the official C library, where it is only the last 2 digits of the year; and make ``.no`` available in the Python fallback case as well.
| 2020-02-02 — 2.1 — Add vectorized array method to Satrec object; add ``.no`` attribute to new Satrec object to support old code that has not migrated to the new name ``.no_kozai``; gave Python wrapper classes ``__slots__`` to avoid the expense of a per-object attribute dictionary.
| 2020-01-30 — 2.0 — Rewrite API to use genuine Vallado C++ code on those systems where it can be compiled; add accelerated vectorized array interface; make ``gstime()`` a public function; clarify format error message.
| 2015-01-15 — 1.4 — Display detailed help when TLE input does not match format.
| 2014-06-26 — 1.3 — Return ``(NaN,NaN,NaN)`` vectors on error and set ``.error_message``
| 2013-11-29 — 1.2 — Made ``epochyr`` 4 digits; add ``datetime`` for ``.epoch``
| 2012-11-22 — 1.1 — Python 3 compatibility; more documentation
| 2012-08-27 — 1.0 — Initial release

"""
