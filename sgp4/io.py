"""Read the TLE earth satellite file format.

This is a minimally-edited copy of "sgp4io.cpp".

"""
import re
from math import pi, pow
from sgp4.ext import days2mdhms, jday
from sgp4.model import Satellite
from sgp4.propagation import sgp4init

SPACE = ord(' ')
INT_RE = re.compile(r'[+-]?\d*')
FLOAT_RE = re.compile(r'[+-]?\d*(\.\d*)?')


"""
/*     ----------------------------------------------------------------
*
*                               sgp4io.cpp
*
*    this file contains a function to read two line element sets. while 
*    not formerly part of the sgp4 mathematical theory, it is 
*    required for practical implemenation.
*
*                            companion code for
*               fundamentals of astrodynamics and applications
*                                    2007
*                              by david vallado
*
*       (w) 719-573-2600, email dvallado@agi.com
*
*    current :
*              27 Aug 10  david vallado
*                           fix input format and delete unused variables in twoline2rv
*    changes :
*               3 sep 08  david vallado
*                           add operationmode for afspc (a) or improved (i)
*               9 may 07  david vallado
*                           fix year correction to 57
*              27 mar 07  david vallado
*                           misc fixes to manual inputs
*              14 aug 06  david vallado
*                           original baseline
*       ----------------------------------------------------------------      */
"""

"""
/* -----------------------------------------------------------------------------
*
*                           function twoline2rv
*
*  this function converts the two line element set character string data to
*    variables and initializes the sgp4 variables. several intermediate varaibles
*    and quantities are determined. note that the result is a structure so multiple
*    satellites can be processed simultaneously without having to reinitialize. the
*    verification mode is an important option that permits quick checks of any
*    changes to the underlying technical theory. this option works using a
*    modified tle file in which the start, stop, and delta time values are
*    included at the end of the second line of data. this only works with the
*    verification mode. the catalog mode simply propagates from -1440 to 1440 min
*    from epoch and is useful when performing entire catalog runs.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs        :
*    longstr1    - first line of the tle
*    longstr2    - second line of the tle
*    typerun     - type of run                    verification 'v', catalog 'c', 
*                                                 manual 'm'
*    typeinput   - type of manual input           mfe 'm', epoch 'e', dayofyr 'd'
*    opsmode     - mode of operation afspc or improved 'a', 'i'
*    whichconst  - which set of constants to use  72, 84
*
*  outputs       :
*    satrec      - structure containing all the sgp4 satellite information
*
*  coupling      :
*    getgravconst-
*    days2mdhms  - conversion of days to month, day, hour, minute, second
*    jday        - convert day month year hour minute second into julian date
*    sgp4init    - initialize the sgp4 variables
*
*  references    :
*    norad spacetrack report #3
*    vallado, crawford, hujsak, kelso  2006
  --------------------------------------------------------------------------- */
"""

def twoline2rv(longstr1, longstr2, whichconst, afspc_mode=False):
       """Return a Satellite imported from two lines of TLE data.

       Provide the two TLE lines as strings `longstr1` and `longstr2`,
       and select which standard set of gravitational constants you want
       by providing `gravity_constants`:

       `sgp4.propagation.wgs72` - Standard WGS 72 model
       `sgp4.propagation.wgs84` - More recent WGS 84 model
       `sgp4.propagation.wgs72old` - Legacy support for old SGP4 behavior

       Normally, computations are made using various recent improvements
       to the algorithm.  If you want to turn some of these off and go
       back into "afspc" mode, then set `afspc_mode` to `True`.

       """
       opsmode = 'a' if afspc_mode else 'i'

       deg2rad  =   pi / 180.0;         #    0.0174532925199433
       xpdotp   =  1440.0 / (2.0 *pi);  #  229.1831180523293

       revnum = 0;
       elnum = 0;

       year = 0;

       tumin = whichconst.tumin

       satrec = Satellite()
       satrec.error = 0;
       satrec.whichconst = whichconst  # Python extension: remembers its consts

       # This is Python, so we make the strings mutable before setting
       # the C++ code loose on them.
       longstr1 = bytearray(longstr1)
       longstr2 = bytearray(longstr2)

       #  set the implied decimal points since doing a formated read
       #  fixes for bad input data values (missing, ...)
       for j in xrange(10, 16):
           if longstr1[j] == SPACE:
               longstr1[j] = '_';

       if longstr1[44] != SPACE:
           longstr1[43] = longstr1[44];
       longstr1[44] = '.';
       if longstr1[7] == SPACE:
           longstr1[7] = 'U';
       if longstr1[9] == SPACE:
           longstr1[9] = '.';
       for j in xrange(45, 50):
           if longstr1[j] == SPACE:
               longstr1[j] = '0';
       if longstr1[51] == SPACE:
           longstr1[51] = '0';
       if longstr1[53] != SPACE:
           longstr1[52] = longstr1[53];
       longstr1[53] = '.';
       longstr2[25] = '.';
       for j in xrange(26, 33):
           if longstr2[j] == SPACE:
               longstr2[j] = '0';
       if longstr1[62] == SPACE:
           longstr1[62] = '0';
       if longstr1[68] == SPACE:
           longstr1[68] = '0';

       # Convert mutable but inconvenient bytearrays into real strings.
       longstr1 = longstr1.decode('ascii')
       longstr2 = longstr2.decode('ascii')

       (cardnumb,satrec.satnum,classification, intldesg, satrec.epochyr,
        satrec.epochdays,satrec.ndot, satrec.nddot, nexp, satrec.bstar,
        ibexp, numb, elnum) = \
       sscanf(longstr1,"%2d %5ld %1c %10s %2d %12lf %11lf %7lf %2d %7lf %2d %2d %6ld ",
           )

       if longstr2[52] == ' ':
           (cardnumb,satrec.satnum, satrec.inclo,
            satrec.nodeo,satrec.ecco, satrec.argpo, satrec.mo, satrec.no,
            revnum) = \
               sscanf(longstr2,"%2d %5ld %9lf %9lf %8lf %9lf %9lf %10lf %6ld \n",
                      )
       else:
           (cardnumb,satrec.satnum, satrec.inclo,
            satrec.nodeo,satrec.ecco, satrec.argpo, satrec.mo, satrec.no,
            revnum) = \
               sscanf(longstr2,"%2d %5ld %9lf %9lf %8lf %9lf %9lf %11lf %6ld \n",
                      )

       #  ---- find no, ndot, nddot ----
       satrec.no   = satrec.no / xpdotp; #   rad/min
       satrec.nddot= satrec.nddot * pow(10.0, nexp);
       satrec.bstar= satrec.bstar * pow(10.0, ibexp);

       #  ---- convert to sgp4 units ----
       satrec.a    = pow( satrec.no*tumin , (-2.0/3.0) );
       satrec.ndot = satrec.ndot  / (xpdotp*1440.0);  #   ? * minperday
       satrec.nddot= satrec.nddot / (xpdotp*1440.0*1440);

       #  ---- find standard orbital elements ----
       satrec.inclo = satrec.inclo  * deg2rad;
       satrec.nodeo = satrec.nodeo  * deg2rad;
       satrec.argpo = satrec.argpo  * deg2rad;
       satrec.mo    = satrec.mo     * deg2rad;

       satrec.alta = satrec.a*(1.0 + satrec.ecco) - 1.0;
       satrec.altp = satrec.a*(1.0 - satrec.ecco) - 1.0;

       """
       // ----------------------------------------------------------------
       // find sgp4epoch time of element set
       // remember that sgp4 uses units of days from 0 jan 1950 (sgp4epoch)
       // and minutes from the epoch (time)
       // ----------------------------------------------------------------

       // ---------------- temp fix for years from 1957-2056 -------------------
       // --------- correct fix will occur when year is 4-digit in tle ---------
       """
       if satrec.epochyr < 57:
           year= satrec.epochyr + 2000;
       else:
           year= satrec.epochyr + 1900;

       mon,day,hr,minute,sec = days2mdhms(year, satrec.epochdays);
       satrec.jdsatepoch = jday(year,mon,day,hr,minute,sec);

       #  ---------------- initialize the orbit at sgp4epoch -------------------
       sgp4init( whichconst, opsmode, satrec.satnum, satrec.jdsatepoch-2433281.5, satrec.bstar,
                 satrec.ecco, satrec.argpo, satrec.inclo, satrec.mo, satrec.no,
                 satrec.nodeo, satrec);

       return satrec


def sscanf(data, format):
    """Yes: a bootleg sscanf(), instead of tediously rewriting the above!"""

    directives = format.split()
    values = []
    start = 0

    for directive in directives:
        conversion = directive[-1]
        fieldwidthstr = directive[1:-1].strip('l')
        fieldwidth = int(fieldwidthstr) if fieldwidthstr else 999

        # scanf(3) skips "any amount of white space, including none"
        while start < len(data) and data[start] == ' ':  # space
            start += 1

        if start == len(data):
            break

        # Field will end early if whitespace, or an illegal character,
        # is encountered.
        if conversion == 'd':
            source = INT_RE.match(data[start:]).group(0)
            end = start + min(len(source), fieldwidth)
        elif conversion == 'f':
            source = FLOAT_RE.match(data[start:]).group(0)
            end = start + min(len(source), fieldwidth)
        else:
            for end in range(start + 1, start + fieldwidth + 1):
                if data[end].isspace():
                    break

        source = data[start:end]

        # Convert!
        if conversion in ('c', 's'):
            values.append(source)
        elif conversion == 'd':
            values.append(int(source))
        elif conversion == 'f':
            values.append(float(source))
        else:
            raise ValueError('unknown format specifier %r' % (conversion,))

        # Start over with the next token.
        start = end

    return values
