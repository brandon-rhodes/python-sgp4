# -*- coding: utf-8 -*-
"""Track earth satellite TLE orbits using up-to-date 2010 version of SGP4

This Python package computes the position and velocity of an
earth-orbiting satellite, given the satellite's TLE orbital elements
from a source like `Celestrak <http://celestrak.com/>`_.  It implements
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

Note that this package produces raw Earth-centered cartesian
coordinates.  It does not implement all the steps necessary to convert
satellite positions into geographic coordinates.  For that, look for a
comprehensive astronomy library that is built atop this one, like the
`Skyfield <http://rhodesmill.org/skyfield/>`_ library:

http://rhodesmill.org/skyfield/earth-satellites.html

To run the test suite for this module, clone its repository from GitHub:

https://github.com/brandon-rhodes/python-sgp4

Then invoke the tests using the Python Standard Library::

    python -m unittest discover sgp4

The C++ function names have been retained, since users may already be
familiar with this library in other languages.  Here is how to compute
the x,y,z position and velocity for the International Space Station at
12:50:19 on 29 June 2000:

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
>>> print(r)
(-6102.44..., -986.33..., -2820.31...)
>>> print(v)
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

* Or, you can provide a coarse date ``jd`` that is within a few weeks or
  months of the satellite’s epoch, and a very precise fraction ``fr``
  that supplies the rest of the value.  The Julian Date for which the
  satellite position is computed is the sum of the two values.  One
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

The attributes of a ``Satrec`` object carry the data loaded from the TLE
entry.  Look at the class's documentation for details.

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

This implementation passes all of the automated tests in the August 2010
release of the reference implementation of SGP4 by Vallado et al., who
originally published their revision of SGP4 in 2006:

    Vallado, David A., Paul Crawford, Richard Hujsak, and T.S. Kelso,
    “Revisiting Spacetrack Report #3,” presented at the AIAA/AAS
    Astrodynamics Specialist Conference, Keystone, CO, 2006 August
    21–24.

If you would like to review the paper, it is `available online
<http://www.celestrak.com/publications/AIAA/2006-6753/>`_.  You can
always download the latest version of their code for comparison against
this Python module (or other implementations) at `AIAA-2006-6753.zip
<http://www.celestrak.com/publications/AIAA/2006-6753/AIAA-2006-6753.zip>`_.

Legacy API
----------

Before this library pivoted to wrapping Vallado's official C++ code and
was operating in pure Python only, it had a slightly quirkier API, which
is still supported for compatibility with older clients.  You can learn
about it by reading the documentation from version 1.4 or earlier:

https://pypi.org/project/sgp4/1.4/

Changelog
---------

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
