"""General-purpose routines.

It seemed a shame for the ``api`` module to have to import large legacy
modules to offer simple date handling, so this small module holds the
routines instead.

"""
def jday(year, mon, day, hr, minute, sec):
     """Return two floats that, when added, produce the specified Julian date.

     The first float specifies the day, while the second float specifies
     an additional offset for the hour, minute, and second.  Because the
     second float is much smaller in magnitude it can, unlike the first
     float, be accurate down to very small fractions of a second.

     """
     jd = (367.0 * year
           - 7 * (year + ((mon + 9) // 12.0)) * 0.25 // 1.0
	   + 275 * mon / 9.0 // 1.0
	   + day
           + 1721013.5)
     fr = (sec + minute * 60.0 + hr * 3600.0) / 86400.0;
     return jd, fr
