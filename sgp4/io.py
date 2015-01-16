"""Read the TLE earth satellite file format.

This is a minimally-edited copy of "sgp4io.cpp".

"""
import re
from datetime import datetime
from math import pi, pow
from sgp4.ext import days2mdhms, jday
from sgp4.model import Satellite
from sgp4.propagation import sgp4init

INT_RE = re.compile(r'[+-]?\d*')
FLOAT_RE = re.compile(r'[+-]?\d*(\.\d*)?')

LINE1 = '1 NNNNNC NNNNNAAA NNNNN.NNNNNNNN +.NNNNNNNN +NNNNN-N +NNNNN-N N NNNNN'
LINE2 = '2 NNNNN NNN.NNNN NNN.NNNN NNNNNNN NNN.NNNN NNN.NNNN NN.NNNNNNNNNNNNNN'

error_message = """TLE format error

The Two-Line Element (TLE) format was designed for punch cards and is
therefore very strict about the position of every space and digit in a
TLE line.  Your line does not quite match.  Here is the official format
for line {} followed by the line you provided:

{}
{}"""

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

    `sgp4.earth_gravity.wgs72` - Standard WGS 72 model
    `sgp4.earth_gravity.wgs84` - More recent WGS 84 model
    `sgp4.earth_gravity.wgs72old` - Legacy support for old SGP4 behavior

    Normally, computations are made using various recent improvements
    to the algorithm.  If you want to turn some of these off and go
    back into "afspc" mode, then set `afspc_mode` to `True`.

    """
    opsmode = 'a' if afspc_mode else 'i'

    deg2rad  =   pi / 180.0;         #    0.0174532925199433
    xpdotp   =  1440.0 / (2.0 *pi);  #  229.1831180523293

    tumin = whichconst.tumin

    satrec = Satellite()
    satrec.error = 0;
    satrec.whichconst = whichconst  # Python extension: remembers its consts

    line = longstr1.rstrip()
    try:
        assert line.startswith('1 ')
        satrec.satnum = int(line[2:7])
        # classification = line[7] or 'U'
        assert line[8] == ' '
        # intldesg = line[9:17]
        two_digit_year = int(line[18:20])
        assert line[23] == '.'
        satrec.epochdays = float(line[20:32])
        assert line[32] == ' '
        assert line[34] == '.'
        satrec.ndot = float(line[33:43])
        assert line[43] == ' '
        satrec.nddot = float(line[44] + '.' + line[45:50])
        nexp = int(line[50:52])
        assert line[52] == ' '
        satrec.bstar = float(line[53] + '.' + line[54:59])
        ibexp = int(line[59:61])
        assert line[61] == ' '
        assert line[63] == ' '
        # numb = int(line[62])
        # elnum = int(line[64:68])
    except (AssertionError, IndexError, ValueError):
        raise ValueError(error_message.format(1, LINE1, line))

    line = longstr2.rstrip()
    try:
        assert line.startswith('2 ')
        satrec.satnum = int(line[2:7])  # TODO: check it matches line 1?
        assert line[7] == ' '
        assert line[11] == '.'
        satrec.inclo = float(line[8:16])
        assert line[16] == ' '
        assert line[20] == '.'
        satrec.nodeo = float(line[17:25])
        assert line[25] == ' '
        satrec.ecco = float('0.' + line[26:33].replace(' ', '0'))
        assert line[33] == ' '
        assert line[37] == '.'
        satrec.argpo = float(line[34:42])
        assert line[42] == ' '
        assert line[46] == '.'
        satrec.mo = float(line[43:51])
        assert line[51] == ' '
        satrec.no = float(line[52:63])
        #revnum = line[63:68]
    except (AssertionError, IndexError, ValueError):
        raise ValueError(error_message.format(2, LINE2, line))

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
    if two_digit_year < 57:
        year = two_digit_year + 2000;
    else:
        year = two_digit_year + 1900;

    mon,day,hr,minute,sec = days2mdhms(year, satrec.epochdays);
    sec_whole, sec_fraction = divmod(sec, 1.0)

    satrec.epochyr = year
    satrec.jdsatepoch = jday(year,mon,day,hr,minute,sec);
    satrec.epoch = datetime(year, mon, day, hr, minute, int(sec_whole),
                            int(sec_fraction * 1000000.0 // 1.0))

    #  ---------------- initialize the orbit at sgp4epoch -------------------
    sgp4init(whichconst, opsmode, satrec.satnum, satrec.jdsatepoch-2433281.5, satrec.bstar,
             satrec.ecco, satrec.argpo, satrec.inclo, satrec.mo, satrec.no,
             satrec.nodeo, satrec)

    return satrec
