# -*- coding: utf-8 -*-
"""Track earth satellite TLE orbits using up-to-date 2010 version of SGP4

This Python package computes the position and velocity of an
earth-orbiting satellite, given the satellite's TLE orbital elements
from a source like `Celestrak <http://celestrak.com/>`_.  It implements
the most recent version of SGP4, and has been tested to make sure that
its satellite position predictions **agree to within 1 µm** of the
standard C++ implementation of the algorithm.

The C++ function names have been retained, since users may already be
familiar with this library from other languages.  The result looks
somewhat clunky in Python, but should be easy to use.  Here are the
x,y,z position and velocity for the Vanguard 1 satellite in 1958:

>>> from sgp4.io import twoline2rv
>>> from sgp4.propagation import sgp4, wgs84
>>>
>>> line1 = ('1 00005U 58002B   00179.78495062  '
...          '.00000023  00000-0  28098-4 0  4753')
>>> line2 = ('2 00005  34.2682 348.7242 1859667 '
...          '331.7664  19.3264 10.82419157413667')
>>>
>>> satellite = twoline2rv(line1, line2, wgs84)
>>> position, velocity = sgp4(wgs84, satellite, 0.0)
>>> position
[7022.46647249137, -1400.0665618178339, 0.05106558274635007]
>>> velocity
[1.8938310806788694, 6.405894872518269, 4.534806700953066]

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
was the most recently updated SPG4 implementation in their zip file:

* C++, August 2010
* Fortran, August 2008
* Pascal, August 2008
* Matlab, May 2008
* Java, July 2005

"""
