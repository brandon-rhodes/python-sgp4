"""General-purpose routines.

It seemed a shame for the ``api`` module to have to import large legacy
modules to offer simple date handling, so this small module holds the
routines instead.

"""
def jday(year, mon, day, hr, minute, sec):
     """Return two floats that, when added, produce the specified Julian date.

     The first float returned gives the date, while the second float
     provides an additional offset for the particular hour, minute, and
     second of that date.  Because the second float is much smaller in
     magnitude it can, unlike the first float, be accurate down to very
     small fractions of a second.

     >>> jd, fr = jday(2020, 2, 11, 13, 57, 0)
     >>> jd
     2458890.5
     >>> fr
     0.58125

     Note that the first float, which gives the moment of midnight that
     commences the given calendar date, always carries the fraction
     ``.5`` because Julian dates begin and end at noon.  This made
     Julian dates more convenient for astronomers in Europe, by making
     the whole night belong to a single Julian date.

     This function is a simple translation to Python of the C++ routine
     ``jday()`` in Vallado's ``SGP4.cpp``.

     """
     jd = (367.0 * year
           - 7 * (year + ((mon + 9) // 12.0)) * 0.25 // 1.0
	   + 275 * mon / 9.0 // 1.0
	   + day
           + 1721013.5)
     fr = (sec + minute * 60.0 + hr * 3600.0) / 86400.0;
     return jd, fr
