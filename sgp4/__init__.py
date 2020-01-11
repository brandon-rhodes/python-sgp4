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

Note that these official SGP4 routines do not implement all the steps
necessary to convert satellite positions into geographic coordinates;
they need to be integrated into a more comprehensive astronomy library.
One example is the Python `Skyfield <http://rhodesmill.org/skyfield/>`_
package, which uses this SGP4 library when you ask it to turn satellite
elements into Earth positions as described in its documentation:

http://rhodesmill.org/skyfield/earth-satellites.html

These SGP4 routines, by contrast, produce only raw spatial coordinates.

To run the test suite for this module, clone its repository from GitHub:

https://github.com/brandon-rhodes/python-sgp4

Then invoke the tests using the Python Standard Library::

    python -m unittest discover sgp4

The C++ function names have been retained, since users may already be
familiar with this library in other languages.  Here is how to compute
the x,y,z position and velocity for Vanguard 1 at 12:50:19 on 29
June 2000:

>>> from sgp4.api import Satrec
>>>
>>> s = '1 25544U 98067A   19343.69339541  .00001764  00000-0  38792-4 0  9991'
>>> t = '2 25544  51.6439 211.2001 0007417  17.6667  85.6398 15.50103472202482'
>>>
>>> satellite = Satrec.twoline2rv(s, t)
>>> jd, fr = 2458827, 0.362605
>>> error, position, velocity = satellite.sgp4(jd, fr)
>>> error
0
>>> print(position)
(-6102.44..., -986.33..., -2820.31...)
>>> print(velocity)
(-1.45..., -5.52..., 5.10...)

>>> from sgp4.api import jday
>>> jd, fr = jday(2019, 12, 9, 12, 0, 0)
>>> jd
2458826.5
>>> fr
0.5

>>> s = '1 20580U 90037B   19342.88042116  .00000361  00000-0  11007-4 0  9996'
>>> t = '2 20580  28.4682 146.6676 0002639 185.9222 322.7238 15.09309432427086'
>>> satellite2 = Satrec.twoline2rv(s, t)


>>> import numpy as np
>>> jd = np.array((2458826, 2458826, 2458826, 2458826))
>>> fr = np.array((0.0001, 0.0002, 0.0003, 0.0004))

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

| 2019-08-?? — 1.5 — Make ``gstime()`` a public function; clarify format error message.
| 2015-01-15 — 1.4 — Display detailed help when TLE input does not match format.
| 2014-06-26 — 1.3 — Return ``(NaN,NaN,NaN)`` vectors on error and set ``.error_message``
| 2013-11-29 — 1.2 — Made ``epochyr`` 4 digits; add ``datetime`` for ``.epoch``
| 2012-11-22 — 1.1 — Python 3 compatibility; more documentation
| 2012-08-27 — 1.0 — Initial release

"""
