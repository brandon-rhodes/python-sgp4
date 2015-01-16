# -*- coding: utf-8 -*-
"""Track earth satellite TLE orbits using up-to-date 2010 version of SGP4

This Python package computes the position and velocity of an
earth-orbiting satellite, given the satellite's TLE orbital elements
from a source like `Celestrak <http://celestrak.com/>`_.  It implements
the most recent version of SGP4, and is regularly run against the SGP4
test suite to make sure that its satellite position predictions **agree
to within 0.1 mm** of the predictions of the standard C++ implementation
of the algorithm.  This error is far less than the 1–3 km/day by which
satellites themselves deviate from the ideal orbits described in TLE
files.

The C++ function names have been retained, since users may already be
familiar with this library in other languages.  Here is how to compute
the x,y,z position and velocity for Vanguard 1 at 12:50:19 on 29
June 2000:

>>> from sgp4.earth_gravity import wgs72
>>> from sgp4.io import twoline2rv
>>>
>>> line1 = ('1 00005U 58002B   00179.78495062  '
...          '.00000023  00000-0  28098-4 0  4753')
>>> line2 = ('2 00005  34.2682 348.7242 1859667 '
...          '331.7664  19.3264 10.82419157413667')
>>>
>>> satellite = twoline2rv(line1, line2, wgs72)
>>> position, velocity = satellite.propagate(
...     2000, 6, 29, 12, 50, 19)
>>>
>>> print(satellite.error)    # nonzero on error
0
>>> print(satellite.error_message)
None
>>> print(position)
(5576.056952..., -3999.371134..., -1521.957159...)
>>> print(velocity)
(4.772627..., 5.119817..., 4.275553...)

The position vector measures the satellite position in **kilometers**
from the center of the earth.  The velocity is the rate at which those
three parameters are changing, expressed in **kilometers per second**.

There are three gravity models available that you can import from the
``earth_gravity`` module:

* ``wgs72``
* ``wgs72old``
* ``wgs84``

The ``wgs72`` model seems to be the most commonly used in the satellite
tracking community, and is probably the model behind most TLE elements
that are available for download.

The ``twoline2rv()`` function returns a ``Satellite`` object whose
attributes carry the data loaded from the TLE entry:

* Unique satellite number, as given in the TLE file.

  >>> satellite.satnum
  5

* The epoch of the element set, expressed three ways:
  as the integer year plus the floating point number of days into the year;
  as a floating-point Julian date; and as Python ``datetime`` object.

  >>> satellite.epochyr
  2000
  >>> satellite.epochdays
  179.78495062
  >>> satellite.jdsatepoch
  2451723.28495062
  >>> satellite.epoch
  datetime.datetime(2000, 6, 27, 18, 50, 19, 733567)

This implementation passes all of the automated tests in the August 2010
release of the reference implementation of SGP4 by Vallado et al., who
originally published their revision of SGP4 in 2006:

    Vallado, David A., Paul Crawford, Richard Hujsak, and T.S. Kelso, “Revisiting Spacetrack Report #3,” presented at the AIAA/AAS Astrodynamics Specialist Conference, Keystone, CO, 2006 August 21–24.

If you would like to review the paper, it is `available online
<http://www.celestrak.com/publications/AIAA/2006-6753/>`_.  You can
always download the latest version of their code for comparison against
this Python module (or other implementations) at `AIAA-2006-6753.zip
<http://www.celestrak.com/publications/AIAA/2006-6753/AIAA-2006-6753.zip>`_.

This module was adapted from Vallado's C++ code since its revision date
was the most recently updated SGP4 implementation in their zip file:

* C++, August 2010
* Fortran, August 2008
* Pascal, August 2008
* Matlab, May 2008
* Java, July 2005

Changelog
---------

| 2015-01-15 — 1.4 — Display detailed help when TLE input does not match format.
| 2014-06-26 — 1.3 — Return ``(NaN,NaN,NaN)`` vectors on error and set ``.error_message``
| 2013-11-29 — 1.2 — Made ``epochyr`` 4 digits; add ``datetime`` for ``.epoch``
| 2012-11-22 — 1.1 — Python 3 compatibility; more documentation
| 2012-08-27 — 1.0 — Initial release

"""
