# -*- coding: utf-8 -*-
"""Track earth satellite TLE orbits using up-to-date 2010 version of SGP4

This package can compute the position of an earth-orbiting satellite
from the TLE orbital elements that you can download from `Celestrak
<http://celestrak.com/>`_.  This implementation passes every one of the
automated tests in the August 2010 release of SGP4 by Vallado et al,
who originally published their revision of the algorithm in 2006:

    Vallado, David A., Paul Crawford, Richard Hujsak, and T.S. Kelso, “Revisiting Spacetrack Report #3,” presented at the AIAA/AAS Astrodynamics Specialist Conference, Keystone, CO, 2006 August 21–24.

If you would like to review the paper, it is `available online
<http://www.celestrak.com/publications/AIAA/2006-6753/>`_.  You can
always download the latest version of their code for comparison against
this Python module (or other implementations) as `AIAA-2006-6753.zip
<http://www.celestrak.com/publications/AIAA/2006-6753/AIAA-2006-6753.zip>`_.

"""
