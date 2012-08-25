# -*- coding: utf-8 -*-
"""Computes earth satellite positions using up-to-date 2006 algorithm

This package can compute the position of an earth-orbiting satellite
from the TLE orbital elements that can be downloaded from `Celestrak
<http://celestrak.com/>`_.  This implementation passes every one of the
automated tests crafted by Vallado et al for their 2006 paper updating
the traditional SGP4 algorithm:

    Vallado, David A., Paul Crawford, Richard Hujsak, and T.S. Kelso, “Revisiting Spacetrack Report #3,” presented at the AIAA/AAS Astrodynamics Specialist Conference, Keystone, CO, 2006 August 21–24.

If you would like to review the paper, it is `available online
<http://www.celestrak.com/publications/AIAA/2006-6753/>`_.

"""
