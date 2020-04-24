"""The sgp4 procedures for analytical propagation of a satellite.

I have made the rather unorthodox decision to leave as much of this C++
code alone as possible: if a line of code would run without change in
Python, then I refused to re-indent it or remove its terminal semicolon,
so that in the future it will be easier to keep updating this file as
the original author's C++ continues to improve.  Thus, 5-space
indentation (!) prevails in this file.

I have even kept all of the C++ block comments (by turning them into
Python string constants) to make this easier to navigate and maintain,
as well as to make it more informative for people who encounter this
code for the first time here in its Python form.

| - Brandon Rhodes
|   Common Grounds Coffee House, Bluffton, Ohio
|   On a very hot August day in 2012

"""
from math import atan2, cos, fabs, pi, sin, sqrt

deg2rad = pi / 180.0;
_nan = float('NaN')
false = (_nan, _nan, _nan)
true = True
twopi = 2.0 * pi

"""
/*     ----------------------------------------------------------------
*
*                               sgp4unit.cpp
*
*    this file contains the sgp4 procedures for analytical propagation
*    of a satellite. the code was originally released in the 1980 and 1986
*    spacetrack papers. a detailed discussion of the theory and history
*    may be found in the 2006 aiaa paper by vallado, crawford, hujsak,
*    and kelso.
*
*                            companion code for
*               fundamentals of astrodynamics and applications
*                                    2013
*                              by david vallado
*
*     (w) 719-573-2600, email dvallado@agi.com, davallado@gmail.com
*
*    current :
*               7 dec 15  david vallado
*                           fix jd, jdfrac
*    changes :
*               3 nov 14  david vallado
*                           update to msvs2013 c++
*              30 aug 10  david vallado
*                           delete unused variables in initl
*                           replace pow integer 2, 3 with multiplies for speed
*               3 nov 08  david vallado
*                           put returns in for error codes
*              29 sep 08  david vallado
*                           fix atime for faster operation in dspace
*                           add operationmode for afspc (a) or improved (i)
*                           performance mode
*              16 jun 08  david vallado
*                           update small eccentricity check
*              16 nov 07  david vallado
*                           misc fixes for better compliance
*              20 apr 07  david vallado
*                           misc fixes for constants
*              11 aug 06  david vallado
*                           chg lyddane choice back to strn3, constants, misc doc
*              15 dec 05  david vallado
*                           misc fixes
*              26 jul 05  david vallado
*                           fixes for paper
*                           note that each fix is preceded by a
*                           comment with "sgp4fix" and an explanation of
*                           what was changed
*              10 aug 04  david vallado
*                           2nd printing baseline working
*              14 may 01  david vallado
*                           2nd edition baseline
*                     80  norad
*                           original baseline
*       ----------------------------------------------------------------      */
"""

"""
/* -----------------------------------------------------------------------------
*
*                           procedure dpper
*
*  this procedure provides deep space long period periodic contributions
*    to the mean elements.  by design, these periodics are zero at epoch.
*    this used to be dscom which included initialization, but it's really a
*    recurring function.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    e3          -
*    ee2         -
*    peo         -
*    pgho        -
*    pho         -
*    pinco       -
*    plo         -
*    se2 , se3 , sgh2, sgh3, sgh4, sh2, sh3, si2, si3, sl2, sl3, sl4 -
*    t           -
*    xh2, xh3, xi2, xi3, xl2, xl3, xl4 -
*    zmol        -
*    zmos        -
*    ep          - eccentricity                           0.0 - 1.0
*    inclo       - inclination - needed for lyddane modification
*    nodep       - right ascension of ascending node
*    argpp       - argument of perigee
*    mp          - mean anomaly
*
*  outputs       :
*    ep          - eccentricity                           0.0 - 1.0
*    inclp       - inclination
*    nodep        - right ascension of ascending node
*    argpp       - argument of perigee
*    mp          - mean anomaly
*
*  locals        :
*    alfdp       -
*    betdp       -
*    cosip  , sinip  , cosop  , sinop  ,
*    dalf        -
*    dbet        -
*    dls         -
*    f2, f3      -
*    pe          -
*    pgh         -
*    ph          -
*    pinc        -
*    pl          -
*    sel   , ses   , sghl  , sghs  , shl   , shs   , sil   , sinzf , sis   ,
*    sll   , sls
*    xls         -
*    xnoh        -
*    zf          -
*    zm          -
*
*  coupling      :
*    none.
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/
"""

def _dpper(satrec, inclo, init, ep, inclp, nodep, argpp, mp, opsmode):

     # Copy satellite attributes into local variables for convenience
     # and symmetry in writing formulae.

     e3 = satrec.e3
     ee2 = satrec.ee2
     peo = satrec.peo
     pgho = satrec.pgho
     pho = satrec.pho
     pinco = satrec.pinco
     plo = satrec.plo
     se2 = satrec.se2
     se3 = satrec.se3
     sgh2 = satrec.sgh2
     sgh3 = satrec.sgh3
     sgh4 = satrec.sgh4
     sh2 = satrec.sh2
     sh3 = satrec.sh3
     si2 = satrec.si2
     si3 = satrec.si3
     sl2 = satrec.sl2
     sl3 = satrec.sl3
     sl4 = satrec.sl4
     t = satrec.t
     xgh2 = satrec.xgh2
     xgh3 = satrec.xgh3
     xgh4 = satrec.xgh4
     xh2 = satrec.xh2
     xh3 = satrec.xh3
     xi2 = satrec.xi2
     xi3 = satrec.xi3
     xl2 = satrec.xl2
     xl3 = satrec.xl3
     xl4 = satrec.xl4
     zmol = satrec.zmol
     zmos = satrec.zmos

     #  ---------------------- constants -----------------------------
     zns   = 1.19459e-5;
     zes   = 0.01675;
     znl   = 1.5835218e-4;
     zel   = 0.05490;

     #  --------------- calculate time varying periodics -----------
     zm    = zmos + zns * t;
     # be sure that the initial call has time set to zero
     if init == 'y':
         zm = zmos;
     zf    = zm + 2.0 * zes * sin(zm);
     sinzf = sin(zf);
     f2    =  0.5 * sinzf * sinzf - 0.25;
     f3    = -0.5 * sinzf * cos(zf);
     ses   = se2* f2 + se3 * f3;
     sis   = si2 * f2 + si3 * f3;
     sls   = sl2 * f2 + sl3 * f3 + sl4 * sinzf;
     sghs  = sgh2 * f2 + sgh3 * f3 + sgh4 * sinzf;
     shs   = sh2 * f2 + sh3 * f3;
     zm    = zmol + znl * t;
     if init == 'y':
         zm = zmol;
     zf    = zm + 2.0 * zel * sin(zm);
     sinzf = sin(zf);
     f2    =  0.5 * sinzf * sinzf - 0.25;
     f3    = -0.5 * sinzf * cos(zf);
     sel   = ee2 * f2 + e3 * f3;
     sil   = xi2 * f2 + xi3 * f3;
     sll   = xl2 * f2 + xl3 * f3 + xl4 * sinzf;
     sghl  = xgh2 * f2 + xgh3 * f3 + xgh4 * sinzf;
     shll  = xh2 * f2 + xh3 * f3;
     pe    = ses + sel;
     pinc  = sis + sil;
     pl    = sls + sll;
     pgh   = sghs + sghl;
     ph    = shs + shll;

     if init == 'n':

       pe    = pe - peo;
       pinc  = pinc - pinco;
       pl    = pl - plo;
       pgh   = pgh - pgho;
       ph    = ph - pho;
       inclp = inclp + pinc;
       ep    = ep + pe;
       sinip = sin(inclp);
       cosip = cos(inclp);

       """
       /* ----------------- apply periodics directly ------------ */
       //  sgp4fix for lyddane choice
       //  strn3 used original inclination - this is technically feasible
       //  gsfc used perturbed inclination - also technically feasible
       //  probably best to readjust the 0.2 limit value and limit discontinuity
       //  0.2 rad = 11.45916 deg
       //  use next line for original strn3 approach and original inclination
       //  if (inclo >= 0.2)
       //  use next line for gsfc version and perturbed inclination
       """

       if inclp >= 0.2:

           ph /= sinip
           pgh -= cosip * ph
           argpp += pgh
           nodep += ph
           mp += pl

       else:

           #  ---- apply periodics with lyddane modification ----
           sinop  = sin(nodep);
           cosop  = cos(nodep);
           alfdp  = sinip * sinop;
           betdp  = sinip * cosop;
           dalf   =  ph * cosop + pinc * cosip * sinop;
           dbet   = -ph * sinop + pinc * cosip * cosop;
           alfdp  = alfdp + dalf;
           betdp  = betdp + dbet;
           nodep  = nodep % twopi if nodep >= 0.0 else -(-nodep % twopi)
           #   sgp4fix for afspc written intrinsic functions
           #  nodep used without a trigonometric function ahead
           if nodep < 0.0 and opsmode == 'a':
               nodep = nodep + twopi;
           xls = mp + argpp + pl + pgh + (cosip - pinc * sinip) * nodep
           xnoh   = nodep;
           nodep  = atan2(alfdp, betdp);
           #   sgp4fix for afspc written intrinsic functions
           #  nodep used without a trigonometric function ahead
           if nodep < 0.0 and opsmode == 'a':
               nodep = nodep + twopi;
           if fabs(xnoh - nodep) > pi:
             if nodep < xnoh:
                nodep = nodep + twopi;
             else:
                nodep = nodep - twopi;
           mp += pl
           argpp = xls - mp - cosip * nodep;

     return ep, inclp, nodep, argpp, mp

"""
/*-----------------------------------------------------------------------------
*
*                           procedure dscom
*
*  this procedure provides deep space common items used by both the secular
*    and periodics subroutines.  input is provided as shown. this routine
*    used to be called dpper, but the functions inside weren't well organized.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    epoch       -
*    ep          - eccentricity
*    argpp       - argument of perigee
*    tc          -
*    inclp       - inclination
*    nodep       - right ascension of ascending node
*    np          - mean motion
*
*  outputs       :
*    sinim  , cosim  , sinomm , cosomm , snodm  , cnodm
*    day         -
*    e3          -
*    ee2         -
*    em          - eccentricity
*    emsq        - eccentricity squared
*    gam         -
*    peo         -
*    pgho        -
*    pho         -
*    pinco       -
*    plo         -
*    rtemsq      -
*    se2, se3         -
*    sgh2, sgh3, sgh4        -
*    sh2, sh3, si2, si3, sl2, sl3, sl4         -
*    s1, s2, s3, s4, s5, s6, s7          -
*    ss1, ss2, ss3, ss4, ss5, ss6, ss7, sz1, sz2, sz3         -
*    sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33        -
*    xgh2, xgh3, xgh4, xh2, xh3, xi2, xi3, xl2, xl3, xl4         -
*    nm          - mean motion
*    z1, z2, z3, z11, z12, z13, z21, z22, z23, z31, z32, z33         -
*    zmol        -
*    zmos        -
*
*  locals        :
*    a1, a2, a3, a4, a5, a6, a7, a8, a9, a10         -
*    betasq      -
*    cc          -
*    ctem, stem        -
*    x1, x2, x3, x4, x5, x6, x7, x8          -
*    xnodce      -
*    xnoi        -
*    zcosg  , zsing  , zcosgl , zsingl , zcosh  , zsinh  , zcoshl , zsinhl ,
*    zcosi  , zsini  , zcosil , zsinil ,
*    zx          -
*    zy          -
*
*  coupling      :
*    none.
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/
"""

def _dscom(
       epoch,  ep,     argpp,   tc,     inclp,
       nodep,  np,
       e3,     ee2,
       peo,    pgho,  pho,
       pinco, plo, se2,   se3,
       sgh2,  sgh3,  sgh4,   sh2,   sh3,
       si2,   si3,   sl2,    sl3,   sl4,
       xgh2,  xgh3,   xgh4,  xh2,
       xh3,   xi2,   xi3,    xl2,   xl3,
       xl4,   zmol,  zmos,
     ):

     #  -------------------------- constants -------------------------
     zes     =  0.01675;
     zel     =  0.05490;
     c1ss    =  2.9864797e-6;
     c1l     =  4.7968065e-7;
     zsinis  =  0.39785416;
     zcosis  =  0.91744867;
     zcosgs  =  0.1945905;
     zsings  = -0.98088458;

     #  --------------------- local variables ------------------------
     nm     = np;
     em     = ep;
     snodm  = sin(nodep);
     cnodm  = cos(nodep);
     sinomm = sin(argpp);
     cosomm = cos(argpp);
     sinim  = sin(inclp);
     cosim  = cos(inclp);
     emsq   = em * em;
     betasq = 1.0 - emsq;
     rtemsq = sqrt(betasq);

     #  ----------------- initialize lunar solar terms ---------------
     peo    = 0.0;
     pinco  = 0.0;
     plo    = 0.0;
     pgho   = 0.0;
     pho    = 0.0;
     day    = epoch + 18261.5 + tc / 1440.0;
     xnodce = (4.5236020 - 9.2422029e-4 * day) % twopi
     stem   = sin(xnodce);
     ctem   = cos(xnodce);
     zcosil = 0.91375164 - 0.03568096 * ctem;
     zsinil = sqrt(1.0 - zcosil * zcosil);
     zsinhl = 0.089683511 * stem / zsinil;
     zcoshl = sqrt(1.0 - zsinhl * zsinhl);
     gam    = 5.8351514 + 0.0019443680 * day;
     zx     = 0.39785416 * stem / zsinil;
     zy     = zcoshl * ctem + 0.91744867 * zsinhl * stem;
     zx     = atan2(zx, zy);
     zx     = gam + zx - xnodce;
     zcosgl = cos(zx);
     zsingl = sin(zx);

     #  ------------------------- do solar terms ---------------------
     zcosg = zcosgs;
     zsing = zsings;
     zcosi = zcosis;
     zsini = zsinis;
     zcosh = cnodm;
     zsinh = snodm;
     cc    = c1ss;
     xnoi  = 1.0 / nm;

     for lsflg in 1, 2:

         a1  =   zcosg * zcosh + zsing * zcosi * zsinh;
         a3  =  -zsing * zcosh + zcosg * zcosi * zsinh;
         a7  =  -zcosg * zsinh + zsing * zcosi * zcosh;
         a8  =   zsing * zsini;
         a9  =   zsing * zsinh + zcosg * zcosi * zcosh;
         a10 =   zcosg * zsini;
         a2  =   cosim * a7 + sinim * a8;
         a4  =   cosim * a9 + sinim * a10;
         a5  =  -sinim * a7 + cosim * a8;
         a6  =  -sinim * a9 + cosim * a10;

         x1  =  a1 * cosomm + a2 * sinomm;
         x2  =  a3 * cosomm + a4 * sinomm;
         x3  = -a1 * sinomm + a2 * cosomm;
         x4  = -a3 * sinomm + a4 * cosomm;
         x5  =  a5 * sinomm;
         x6  =  a6 * sinomm;
         x7  =  a5 * cosomm;
         x8  =  a6 * cosomm;

         z31 = 12.0 * x1 * x1 - 3.0 * x3 * x3;
         z32 = 24.0 * x1 * x2 - 6.0 * x3 * x4;
         z33 = 12.0 * x2 * x2 - 3.0 * x4 * x4;
         z1  =  3.0 *  (a1 * a1 + a2 * a2) + z31 * emsq;
         z2  =  6.0 *  (a1 * a3 + a2 * a4) + z32 * emsq;
         z3  =  3.0 *  (a3 * a3 + a4 * a4) + z33 * emsq;
         z11 = -6.0 * a1 * a5 + emsq *  (-24.0 * x1 * x7-6.0 * x3 * x5);
         z12 = -6.0 *  (a1 * a6 + a3 * a5) + emsq * \
                (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5));
         z13 = -6.0 * a3 * a6 + emsq * (-24.0 * x2 * x8 - 6.0 * x4 * x6);
         z21 =  6.0 * a2 * a5 + emsq * (24.0 * x1 * x5 - 6.0 * x3 * x7);
         z22 =  6.0 *  (a4 * a5 + a2 * a6) + emsq * \
                (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8));
         z23 =  6.0 * a4 * a6 + emsq * (24.0 * x2 * x6 - 6.0 * x4 * x8);
         z1  = z1 + z1 + betasq * z31;
         z2  = z2 + z2 + betasq * z32;
         z3  = z3 + z3 + betasq * z33;
         s3  = cc * xnoi;
         s2  = -0.5 * s3 / rtemsq;
         s4  = s3 * rtemsq;
         s1  = -15.0 * em * s4;
         s5  = x1 * x3 + x2 * x4;
         s6  = x2 * x3 + x1 * x4;
         s7  = x2 * x4 - x1 * x3;

         #  ----------------------- do lunar terms -------------------
         if lsflg == 1:

             ss1   = s1;
             ss2   = s2;
             ss3   = s3;
             ss4   = s4;
             ss5   = s5;
             ss6   = s6;
             ss7   = s7;
             sz1   = z1;
             sz2   = z2;
             sz3   = z3;
             sz11  = z11;
             sz12  = z12;
             sz13  = z13;
             sz21  = z21;
             sz22  = z22;
             sz23  = z23;
             sz31  = z31;
             sz32  = z32;
             sz33  = z33;
             zcosg = zcosgl;
             zsing = zsingl;
             zcosi = zcosil;
             zsini = zsinil;
             zcosh = zcoshl * cnodm + zsinhl * snodm;
             zsinh = snodm * zcoshl - cnodm * zsinhl;
             cc    = c1l;

     zmol = (4.7199672 + 0.22997150  * day - gam) % twopi
     zmos = (6.2565837 + 0.017201977 * day) % twopi

     #  ------------------------ do solar terms ----------------------
     se2  =   2.0 * ss1 * ss6;
     se3  =   2.0 * ss1 * ss7;
     si2  =   2.0 * ss2 * sz12;
     si3  =   2.0 * ss2 * (sz13 - sz11);
     sl2  =  -2.0 * ss3 * sz2;
     sl3  =  -2.0 * ss3 * (sz3 - sz1);
     sl4  =  -2.0 * ss3 * (-21.0 - 9.0 * emsq) * zes;
     sgh2 =   2.0 * ss4 * sz32;
     sgh3 =   2.0 * ss4 * (sz33 - sz31);
     sgh4 = -18.0 * ss4 * zes;
     sh2  =  -2.0 * ss2 * sz22;
     sh3  =  -2.0 * ss2 * (sz23 - sz21);

     #  ------------------------ do lunar terms ----------------------
     ee2  =   2.0 * s1 * s6;
     e3   =   2.0 * s1 * s7;
     xi2  =   2.0 * s2 * z12;
     xi3  =   2.0 * s2 * (z13 - z11);
     xl2  =  -2.0 * s3 * z2;
     xl3  =  -2.0 * s3 * (z3 - z1);
     xl4  =  -2.0 * s3 * (-21.0 - 9.0 * emsq) * zel;
     xgh2 =   2.0 * s4 * z32;
     xgh3 =   2.0 * s4 * (z33 - z31);
     xgh4 = -18.0 * s4 * zel;
     xh2  =  -2.0 * s2 * z22;
     xh3  =  -2.0 * s2 * (z23 - z21);

     return (
       snodm, cnodm, sinim,  cosim, sinomm,
       cosomm,day,   e3,     ee2,   em,
       emsq,  gam,   peo,    pgho,  pho,
       pinco, plo,   rtemsq, se2,   se3,
       sgh2,  sgh3,  sgh4,   sh2,   sh3,
       si2,   si3,   sl2,    sl3,   sl4,
       s1,    s2,    s3,     s4,    s5,
       s6,    s7,    ss1,    ss2,   ss3,
       ss4,   ss5,   ss6,    ss7,   sz1,
       sz2,   sz3,   sz11,   sz12,  sz13,
       sz21,  sz22,  sz23,   sz31,  sz32,
       sz33,  xgh2,  xgh3,   xgh4,  xh2,
       xh3,   xi2,   xi3,    xl2,   xl3,
       xl4,   nm,    z1,     z2,    z3,
       z11,   z12,   z13,    z21,   z22,
       z23,   z31,   z32,    z33,   zmol,
       zmos
       )

"""
/*-----------------------------------------------------------------------------
*
*                           procedure dsinit
*
*  this procedure provides deep space contributions to mean motion dot due
*    to geopotential resonance with half day and one day orbits.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    cosim, sinim-
*    emsq        - eccentricity squared
*    argpo       - argument of perigee
*    s1, s2, s3, s4, s5      -
*    ss1, ss2, ss3, ss4, ss5 -
*    sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33 -
*    t           - time
*    tc          -
*    gsto        - greenwich sidereal time                   rad
*    mo          - mean anomaly
*    mdot        - mean anomaly dot (rate)
*    no          - mean motion
*    nodeo       - right ascension of ascending node
*    nodedot     - right ascension of ascending node dot (rate)
*    xpidot      -
*    z1, z3, z11, z13, z21, z23, z31, z33 -
*    eccm        - eccentricity
*    argpm       - argument of perigee
*    inclm       - inclination
*    mm          - mean anomaly
*    xn          - mean motion
*    nodem       - right ascension of ascending node
*
*  outputs       :
*    em          - eccentricity
*    argpm       - argument of perigee
*    inclm       - inclination
*    mm          - mean anomaly
*    nm          - mean motion
*    nodem       - right ascension of ascending node
*    irez        - flag for resonance           0-none, 1-one day, 2-half day
*    atime       -
*    d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433    -
*    dedt        -
*    didt        -
*    dmdt        -
*    dndt        -
*    dnodt       -
*    domdt       -
*    del1, del2, del3        -
*    ses  , sghl , sghs , sgs  , shl  , shs  , sis  , sls
*    theta       -
*    xfact       -
*    xlamo       -
*    xli         -
*    xni
*
*  locals        :
*    ainv2       -
*    aonv        -
*    cosisq      -
*    eoc         -
*    f220, f221, f311, f321, f322, f330, f441, f442, f522, f523, f542, f543  -
*    g200, g201, g211, g300, g310, g322, g410, g422, g520, g521, g532, g533  -
*    sini2       -
*    temp        -
*    temp1       -
*    theta       -
*    xno2        -
*
*  coupling      :
*    getgravconst
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/
"""

def _dsinit(
       # sgp4fix no longer needed pass in xke
       # whichconst,
       xke,
       cosim,  emsq,   argpo,   s1,     s2,
       s3,     s4,     s5,      sinim,  ss1,
       ss2,    ss3,    ss4,     ss5,    sz1,
       sz3,    sz11,   sz13,    sz21,   sz23,
       sz31,   sz33,   t,       tc,     gsto,
       mo,     mdot,   no,      nodeo,  nodedot,
       xpidot, z1,     z3,      z11,    z13,
       z21,    z23,    z31,     z33,    ecco,
       eccsq,  em,    argpm,  inclm, mm,
       nm,    nodem,
       irez,
       atime, d2201, d2211,  d3210, d3222,
       d4410, d4422, d5220,  d5232, d5421,
       d5433, dedt,  didt,   dmdt,
       dnodt, domdt, del1,   del2,  del3,
       xfact, xlamo, xli,    xni,
     ):

     q22    = 1.7891679e-6;
     q31    = 2.1460748e-6;
     q33    = 2.2123015e-7;
     root22 = 1.7891679e-6;
     root44 = 7.3636953e-9;
     root54 = 2.1765803e-9;
     rptim  = 4.37526908801129966e-3; # equates to 7.29211514668855e-5 rad/sec
     root32 = 3.7393792e-7;
     root52 = 1.1428639e-7;
     x2o3   = 2.0 / 3.0;
     znl    = 1.5835218e-4;
     zns    = 1.19459e-5;

	 # sgp4fix identify constants and allow alternate values
	 # just xke is used here so pass it in rather than have multiple calls
     # xke = whichconst.xke

     #  -------------------- deep space initialization ------------
     irez = 0;
     if 0.0034906585 < nm < 0.0052359877:
         irez = 1;
     if 8.26e-3 <= nm <= 9.24e-3 and em >= 0.5:
         irez = 2;

     #  ------------------------ do solar terms -------------------
     ses  =  ss1 * zns * ss5;
     sis  =  ss2 * zns * (sz11 + sz13);
     sls  = -zns * ss3 * (sz1 + sz3 - 14.0 - 6.0 * emsq);
     sghs =  ss4 * zns * (sz31 + sz33 - 6.0);
     shs  = -zns * ss2 * (sz21 + sz23);
     #  sgp4fix for 180 deg incl
     if inclm < 5.2359877e-2 or inclm > pi - 5.2359877e-2:
       shs = 0.0;
     if sinim != 0.0:
       shs = shs / sinim;
     sgs  = sghs - cosim * shs;

     #  ------------------------- do lunar terms ------------------
     dedt = ses + s1 * znl * s5;
     didt = sis + s2 * znl * (z11 + z13);
     dmdt = sls - znl * s3 * (z1 + z3 - 14.0 - 6.0 * emsq);
     sghl = s4 * znl * (z31 + z33 - 6.0);
     shll = -znl * s2 * (z21 + z23);
     #  sgp4fix for 180 deg incl
     if inclm < 5.2359877e-2 or inclm > pi - 5.2359877e-2:
         shll = 0.0;
     domdt = sgs + sghl;
     dnodt = shs;
     if sinim != 0.0:

         domdt = domdt - cosim / sinim * shll;
         dnodt = dnodt + shll / sinim;


     #  ----------- calculate deep space resonance effects --------
     dndt   = 0.0;
     theta  = (gsto + tc * rptim) % twopi
     em     = em + dedt * t;
     inclm  = inclm + didt * t;
     argpm  = argpm + domdt * t;
     nodem  = nodem + dnodt * t;
     mm     = mm + dmdt * t;
     """
     //   sgp4fix for negative inclinations
     //   the following if statement should be commented out
     //if (inclm < 0.0)
     //  {
     //    inclm  = -inclm;
     //    argpm  = argpm - pi;
     //    nodem = nodem + pi;
     //  }
     """

     #  -------------- initialize the resonance terms -------------
     if irez != 0:

         aonv = pow(nm / xke, x2o3);

         #  ---------- geopotential resonance for 12 hour orbits ------
         if irez == 2:

             cosisq = cosim * cosim;
             emo    = em;
             em     = ecco;
             emsqo  = emsq;
             emsq   = eccsq;
             eoc    = em * emsq;
             g201   = -0.306 - (em - 0.64) * 0.440;

             if em <= 0.65:

                 g211 =    3.616  -  13.2470 * em +  16.2900 * emsq;
                 g310 =  -19.302  + 117.3900 * em - 228.4190 * emsq +  156.5910 * eoc;
                 g322 =  -18.9068 + 109.7927 * em - 214.6334 * emsq +  146.5816 * eoc;
                 g410 =  -41.122  + 242.6940 * em - 471.0940 * emsq +  313.9530 * eoc;
                 g422 = -146.407  + 841.8800 * em - 1629.014 * emsq + 1083.4350 * eoc;
                 g520 = -532.114  + 3017.977 * em - 5740.032 * emsq + 3708.2760 * eoc;

             else:

                 g211 =   -72.099 +   331.819 * em -   508.738 * emsq +   266.724 * eoc;
                 g310 =  -346.844 +  1582.851 * em -  2415.925 * emsq +  1246.113 * eoc;
                 g322 =  -342.585 +  1554.908 * em -  2366.899 * emsq +  1215.972 * eoc;
                 g410 = -1052.797 +  4758.686 * em -  7193.992 * emsq +  3651.957 * eoc;
                 g422 = -3581.690 + 16178.110 * em - 24462.770 * emsq + 12422.520 * eoc;
                 if em > 0.715:
                     g520 =-5149.66 + 29936.92 * em - 54087.36 * emsq + 31324.56 * eoc;
                 else:
                     g520 = 1464.74 -  4664.75 * em +  3763.64 * emsq;

             if em < 0.7:

                 g533 = -919.22770 + 4988.6100 * em - 9064.7700 * emsq + 5542.21  * eoc;
                 g521 = -822.71072 + 4568.6173 * em - 8491.4146 * emsq + 5337.524 * eoc;
                 g532 = -853.66600 + 4690.2500 * em - 8624.7700 * emsq + 5341.4  * eoc;

             else:

                 g533 =-37995.780 + 161616.52 * em - 229838.20 * emsq + 109377.94 * eoc;
                 g521 =-51752.104 + 218913.95 * em - 309468.16 * emsq + 146349.42 * eoc;
                 g532 =-40023.880 + 170470.89 * em - 242699.48 * emsq + 115605.82 * eoc;

             sini2=  sinim * sinim;
             f220 =  0.75 * (1.0 + 2.0 * cosim+cosisq);
             f221 =  1.5 * sini2;
             f321 =  1.875 * sinim  *  (1.0 - 2.0 * cosim - 3.0 * cosisq);
             f322 = -1.875 * sinim  *  (1.0 + 2.0 * cosim - 3.0 * cosisq);
             f441 = 35.0 * sini2 * f220;
             f442 = 39.3750 * sini2 * sini2;
             f522 =  9.84375 * sinim * (sini2 * (1.0 - 2.0 * cosim- 5.0 * cosisq) +
                     0.33333333 * (-2.0 + 4.0 * cosim + 6.0 * cosisq) );
             f523 = sinim * (4.92187512 * sini2 * (-2.0 - 4.0 * cosim +
                    10.0 * cosisq) + 6.56250012 * (1.0+2.0 * cosim - 3.0 * cosisq));
             f542 = 29.53125 * sinim * (2.0 - 8.0 * cosim+cosisq *
                    (-12.0 + 8.0 * cosim + 10.0 * cosisq));
             f543 = 29.53125 * sinim * (-2.0 - 8.0 * cosim+cosisq *
                    (12.0 + 8.0 * cosim - 10.0 * cosisq));
             xno2  =  nm * nm;
             ainv2 =  aonv * aonv;
             temp1 =  3.0 * xno2 * ainv2;
             temp  =  temp1 * root22;
             d2201 =  temp * f220 * g201;
             d2211 =  temp * f221 * g211;
             temp1 =  temp1 * aonv;
             temp  =  temp1 * root32;
             d3210 =  temp * f321 * g310;
             d3222 =  temp * f322 * g322;
             temp1 =  temp1 * aonv;
             temp  =  2.0 * temp1 * root44;
             d4410 =  temp * f441 * g410;
             d4422 =  temp * f442 * g422;
             temp1 =  temp1 * aonv;
             temp  =  temp1 * root52;
             d5220 =  temp * f522 * g520;
             d5232 =  temp * f523 * g532;
             temp  =  2.0 * temp1 * root54;
             d5421 =  temp * f542 * g521;
             d5433 =  temp * f543 * g533;
             xlamo =  (mo + nodeo + nodeo-theta - theta) % twopi
             xfact =  mdot + dmdt + 2.0 * (nodedot + dnodt - rptim) - no;
             em    = emo;
             emsq  = emsqo;

         #  ---------------- synchronous resonance terms --------------
         if irez == 1:

             g200  = 1.0 + emsq * (-2.5 + 0.8125 * emsq);
             g310  = 1.0 + 2.0 * emsq;
             g300  = 1.0 + emsq * (-6.0 + 6.60937 * emsq);
             f220  = 0.75 * (1.0 + cosim) * (1.0 + cosim);
             f311  = 0.9375 * sinim * sinim * (1.0 + 3.0 * cosim) - 0.75 * (1.0 + cosim);
             f330  = 1.0 + cosim;
             f330  = 1.875 * f330 * f330 * f330;
             del1  = 3.0 * nm * nm * aonv * aonv;
             del2  = 2.0 * del1 * f220 * g200 * q22;
             del3  = 3.0 * del1 * f330 * g300 * q33 * aonv;
             del1  = del1 * f311 * g310 * q31 * aonv;
             xlamo = (mo + nodeo + argpo - theta) % twopi
             xfact = mdot + xpidot - rptim + dmdt + domdt + dnodt - no;

         #  ------------ for sgp4, initialize the integrator ----------
         xli   = xlamo;
         xni   = no;
         atime = 0.0;
         nm    = no + dndt;

     return (
       em,    argpm,  inclm, mm,
       nm,    nodem,
       irez, atime,
       d2201, d2211,  d3210, d3222,
       d4410, d4422, d5220,  d5232,
       d5421, d5433, dedt,  didt,
       dmdt,  dndt, dnodt, domdt,
       del1,   del2,  del3, xfact,
       xlamo, xli,    xni,
       )

"""
/*-----------------------------------------------------------------------------
*
*                           procedure dspace
*
*  this procedure provides deep space contributions to mean elements for
*    perturbing third body.  these effects have been averaged over one
*    revolution of the sun and moon.  for earth resonance effects, the
*    effects have been averaged over no revolutions of the satellite.
*    (mean motion)
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433 -
*    dedt        -
*    del1, del2, del3  -
*    didt        -
*    dmdt        -
*    dnodt       -
*    domdt       -
*    irez        - flag for resonance           0-none, 1-one day, 2-half day
*    argpo       - argument of perigee
*    argpdot     - argument of perigee dot (rate)
*    t           - time
*    tc          -
*    gsto        - gst
*    xfact       -
*    xlamo       -
*    no          - mean motion
*    atime       -
*    em          - eccentricity
*    ft          -
*    argpm       - argument of perigee
*    inclm       - inclination
*    xli         -
*    mm          - mean anomaly
*    xni         - mean motion
*    nodem       - right ascension of ascending node
*
*  outputs       :
*    atime       -
*    em          - eccentricity
*    argpm       - argument of perigee
*    inclm       - inclination
*    xli         -
*    mm          - mean anomaly
*    xni         -
*    nodem       - right ascension of ascending node
*    dndt        -
*    nm          - mean motion
*
*  locals        :
*    delt        -
*    ft          -
*    theta       -
*    x2li        -
*    x2omi       -
*    xl          -
*    xldot       -
*    xnddt       -
*    xndt        -
*    xomi        -
*
*  coupling      :
*    none        -
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/
"""

def _dspace(
       irez,
       d2201,  d2211,  d3210,   d3222,  d4410,
       d4422,  d5220,  d5232,   d5421,  d5433,
       dedt,   del1,   del2,    del3,   didt,
       dmdt,   dnodt,  domdt,   argpo,  argpdot,
       t,      tc,     gsto,    xfact,  xlamo,
       no,
       atime, em,    argpm,  inclm, xli,
       mm,    xni,   nodem,  nm,
       ):

     fasx2 = 0.13130908;
     fasx4 = 2.8843198;
     fasx6 = 0.37448087;
     g22   = 5.7686396;
     g32   = 0.95240898;
     g44   = 1.8014998;
     g52   = 1.0508330;
     g54   = 4.4108898;
     rptim = 4.37526908801129966e-3; # equates to 7.29211514668855e-5 rad/sec
     stepp =    720.0;
     stepn =   -720.0;
     step2 = 259200.0;

     #  ----------- calculate deep space resonance effects -----------
     dndt   = 0.0;
     theta  = (gsto + tc * rptim) % twopi
     em     = em + dedt * t;

     inclm  = inclm + didt * t;
     argpm  = argpm + domdt * t;
     nodem  = nodem + dnodt * t;
     mm     = mm + dmdt * t;

     """
     //   sgp4fix for negative inclinations
     //   the following if statement should be commented out
     //  if (inclm < 0.0)
     // {
     //    inclm = -inclm;
     //    argpm = argpm - pi;
     //    nodem = nodem + pi;
     //  }

     /* - update resonances : numerical (euler-maclaurin) integration - */
     /* ------------------------- epoch restart ----------------------  */
     //   sgp4fix for propagator problems
     //   the following integration works for negative time steps and periods
     //   the specific changes are unknown because the original code was so convoluted

     // sgp4fix take out atime = 0.0 and fix for faster operation
     """
     ft    = 0.0;
     if irez != 0:

         #  sgp4fix streamline check
         if atime == 0.0 or t * atime <= 0.0 or fabs(t) < fabs(atime):

             atime  = 0.0;
             xni    = no;
             xli    = xlamo;

         # sgp4fix move check outside loop
         if t > 0.0:
               delt = stepp;
         else:
               delt = stepn;

         iretn = 381; # added for do loop
         # iret  =   0; # added for loop
         while iretn == 381:

             #  ------------------- dot terms calculated -------------
             #  ----------- near - synchronous resonance terms -------
             if irez != 2:

                 xndt  = del1 * sin(xli - fasx2) + del2 * sin(2.0 * (xli - fasx4)) + \
                         del3 * sin(3.0 * (xli - fasx6));
                 xldot = xni + xfact;
                 xnddt = del1 * cos(xli - fasx2) + \
                         2.0 * del2 * cos(2.0 * (xli - fasx4)) + \
                         3.0 * del3 * cos(3.0 * (xli - fasx6));
                 xnddt = xnddt * xldot;

             else:

                 # --------- near - half-day resonance terms --------
                 xomi  = argpo + argpdot * atime;
                 x2omi = xomi + xomi;
                 x2li  = xli + xli;
                 xndt  = (d2201 * sin(x2omi + xli - g22) + d2211 * sin(xli - g22) +
                       d3210 * sin(xomi + xli - g32)  + d3222 * sin(-xomi + xli - g32)+
                       d4410 * sin(x2omi + x2li - g44)+ d4422 * sin(x2li - g44) +
                       d5220 * sin(xomi + xli - g52)  + d5232 * sin(-xomi + xli - g52)+
                       d5421 * sin(xomi + x2li - g54) + d5433 * sin(-xomi + x2li - g54));
                 xldot = xni + xfact;
                 xnddt = (d2201 * cos(x2omi + xli - g22) + d2211 * cos(xli - g22) +
                       d3210 * cos(xomi + xli - g32) + d3222 * cos(-xomi + xli - g32) +
                       d5220 * cos(xomi + xli - g52) + d5232 * cos(-xomi + xli - g52) +
                       2.0 * (d4410 * cos(x2omi + x2li - g44) +
                       d4422 * cos(x2li - g44) + d5421 * cos(xomi + x2li - g54) +
                       d5433 * cos(-xomi + x2li - g54)));
                 xnddt = xnddt * xldot;

             #  ----------------------- integrator -------------------
             #  sgp4fix move end checks to end of routine
             if fabs(t - atime) >= stepp:
                 # iret  = 0;
                 iretn = 381;

             else:
                 ft    = t - atime;
                 iretn = 0;

             if iretn == 381:

                 xli   = xli + xldot * delt + xndt * step2;
                 xni   = xni + xndt * delt + xnddt * step2;
                 atime = atime + delt;

         nm = xni + xndt * ft + xnddt * ft * ft * 0.5;
         xl = xli + xldot * ft + xndt * ft * ft * 0.5;
         if irez != 1:
             mm   = xl - 2.0 * nodem + 2.0 * theta;
             dndt = nm - no;

         else:
             mm   = xl - nodem - argpm + theta;
             dndt = nm - no;

         nm = no + dndt;

     return (
       atime, em,    argpm,  inclm, xli,
       mm,    xni,   nodem,  dndt,  nm,
       )

"""
/*-----------------------------------------------------------------------------
*
*                           procedure initl
*
*  this procedure initializes the spg4 propagator. all the initialization is
*    consolidated here instead of having multiple loops inside other routines.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    satn        - satellite number - not needed, placed in satrec
*    xke         - reciprocal of tumin
*    j2          - j2 zonal harmonic
*    ecco        - eccentricity                           0.0 - 1.0
*    epoch       - epoch time in days from jan 0, 1950. 0 hr
*    inclo       - inclination of satellite
*    no          - mean motion of satellite
*
*  outputs       :
*    ainv        - 1.0 / a
*    ao          - semi major axis
*    con41       -
*    con42       - 1.0 - 5.0 cos(i)
*    cosio       - cosine of inclination
*    cosio2      - cosio squared
*    eccsq       - eccentricity squared
*    method      - flag for deep space                    'd', 'n'
*    omeosq      - 1.0 - ecco * ecco
*    posq        - semi-parameter squared
*    rp          - radius of perigee
*    rteosq      - square root of (1.0 - ecco*ecco)
*    sinio       - sine of inclination
*    gsto        - gst at time of observation               rad
*    no          - mean motion of satellite
*
*  locals        :
*    ak          -
*    d1          -
*    del         -
*    adel        -
*    po          -
*
*  coupling      :
*    getgravconst- no longer used
*    gstime      - find greenwich sidereal time from the julian date
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/
"""

def _initl(
       # not needeed. included in satrec if needed later
       # satn,
       # sgp4fix assin xke and j2
       # whichconst,
       xke, j2,
       ecco,   epoch,  inclo,   no,
       method,
       opsmode,
       ):

     # sgp4fix use old way of finding gst

     #  ----------------------- earth constants ----------------------
     #  sgp4fix identify constants and allow alternate values
	 #  only xke and j2 are used here so pass them in directly
     #  tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 = whichconst
     x2o3   = 2.0 / 3.0;

     #  ------------- calculate auxillary epoch quantities ----------
     eccsq  = ecco * ecco;
     omeosq = 1.0 - eccsq;
     rteosq = sqrt(omeosq);
     cosio  = cos(inclo);
     cosio2 = cosio * cosio;

     #  ------------------ un-kozai the mean motion -----------------
     ak    = pow(xke / no, x2o3);
     d1    = 0.75 * j2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq);
     del_  = d1 / (ak * ak);
     adel  = ak * (1.0 - del_ * del_ - del_ *
             (1.0 / 3.0 + 134.0 * del_ * del_ / 81.0));
     del_  = d1/(adel * adel);
     no    = no / (1.0 + del_);

     ao    = pow(xke / no, x2o3);
     sinio = sin(inclo);
     po    = ao * omeosq;
     con42 = 1.0 - 5.0 * cosio2;
     con41 = -con42-cosio2-cosio2;
     ainv  = 1.0 / ao;
     posq  = po * po;
     rp    = ao * (1.0 - ecco);
     method = 'n';

     #  sgp4fix modern approach to finding sidereal time
     if opsmode == 'a':

         #  sgp4fix use old way of finding gst
         #  count integer number of days from 0 jan 1970
         ts70  = epoch - 7305.0;
         ds70 = (ts70 + 1.0e-8) // 1.0;
         tfrac = ts70 - ds70;
         #  find greenwich location at epoch
         c1    = 1.72027916940703639e-2;
         thgr70= 1.7321343856509374;
         fk5r  = 5.07551419432269442e-15;
         c1p2p = c1 + twopi;
         gsto  = (thgr70 + c1*ds70 + c1p2p*tfrac + ts70*ts70*fk5r) % twopi
         if gsto < 0.0:
             gsto = gsto + twopi;

     else:
        gsto = _gstime(epoch + 2433281.5);

     return (
       no,
       method,
       ainv,  ao,    con41,  con42, cosio,
       cosio2,eccsq, omeosq, posq,
       rp,    rteosq,sinio , gsto,
       )

"""
/*-----------------------------------------------------------------------------
*
*                             procedure sgp4init
*
*  this procedure initializes variables for sgp4.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    opsmode     - mode of operation afspc or improved 'a', 'i'
*    whichconst  - which set of constants to use  72, 84
*    satn        - satellite number
*    bstar       - sgp4 type drag coefficient              kg/m2er
*    ecco        - eccentricity
*    epoch       - epoch time in days from jan 0, 1950. 0 hr
*    argpo       - argument of perigee (output if ds)
*    inclo       - inclination
*    mo          - mean anomaly (output if ds)
*    no          - mean motion
*    nodeo       - right ascension of ascending node
*
*  outputs       :
*    satrec      - common values for subsequent calls
*    return code - non-zero on error.
*                   1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
*                   2 - mean motion less than 0.0
*                   3 - pert elements, ecc < 0.0  or  ecc > 1.0
*                   4 - semi-latus rectum < 0.0
*                   5 - epoch elements are sub-orbital
*                   6 - satellite has decayed
*
*  locals        :
*    cnodm  , snodm  , cosim  , sinim  , cosomm , sinomm
*    cc1sq  , cc2    , cc3
*    coef   , coef1
*    cosio4      -
*    day         -
*    dndt        -
*    em          - eccentricity
*    emsq        - eccentricity squared
*    eeta        -
*    etasq       -
*    gam         -
*    argpm       - argument of perigee
*    nodem       -
*    inclm       - inclination
*    mm          - mean anomaly
*    nm          - mean motion
*    perige      - perigee
*    pinvsq      -
*    psisq       -
*    qzms24      -
*    rtemsq      -
*    s1, s2, s3, s4, s5, s6, s7          -
*    sfour       -
*    ss1, ss2, ss3, ss4, ss5, ss6, ss7         -
*    sz1, sz2, sz3
*    sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33        -
*    tc          -
*    temp        -
*    temp1, temp2, temp3       -
*    tsi         -
*    xpidot      -
*    xhdot1      -
*    z1, z2, z3          -
*    z11, z12, z13, z21, z22, z23, z31, z32, z33         -
*
*  coupling      :
*    getgravconst-
*    initl       -
*    dscom       -
*    dpper       -
*    dsinit      -
*    sgp4        -
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/
"""

def sgp4init(
       whichconst,   opsmode,   satn,   epoch,
       xbstar,   xndot,   xnddot,   xecco,   xargpo,
       xinclo,   xmo,   xno_kozai,
       xnodeo,  satrec,
       ):

     """
     /* ------------------------ initialization --------------------- */
     // sgp4fix divisor for divide by zero check on inclination
     // the old check used 1.0 + cos(pi-1.0e-9), but then compared it to
     // 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
     """
     temp4    =   1.5e-12;

     #  ----------- set all near earth variables to zero ------------
     satrec.isimp   = 0;   satrec.method = 'n'; satrec.aycof    = 0.0;
     satrec.con41   = 0.0; satrec.cc1    = 0.0; satrec.cc4      = 0.0;
     satrec.cc5     = 0.0; satrec.d2     = 0.0; satrec.d3       = 0.0;
     satrec.d4      = 0.0; satrec.delmo  = 0.0; satrec.eta      = 0.0;
     satrec.argpdot = 0.0; satrec.omgcof = 0.0; satrec.sinmao   = 0.0;
     satrec.t       = 0.0; satrec.t2cof  = 0.0; satrec.t3cof    = 0.0;
     satrec.t4cof   = 0.0; satrec.t5cof  = 0.0; satrec.x1mth2   = 0.0;
     satrec.x7thm1  = 0.0; satrec.mdot   = 0.0; satrec.nodedot  = 0.0;
     satrec.xlcof   = 0.0; satrec.xmcof  = 0.0; satrec.nodecf   = 0.0;

     #  ----------- set all deep space variables to zero ------------
     satrec.irez  = 0;   satrec.d2201 = 0.0; satrec.d2211 = 0.0;
     satrec.d3210 = 0.0; satrec.d3222 = 0.0; satrec.d4410 = 0.0;
     satrec.d4422 = 0.0; satrec.d5220 = 0.0; satrec.d5232 = 0.0;
     satrec.d5421 = 0.0; satrec.d5433 = 0.0; satrec.dedt  = 0.0;
     satrec.del1  = 0.0; satrec.del2  = 0.0; satrec.del3  = 0.0;
     satrec.didt  = 0.0; satrec.dmdt  = 0.0; satrec.dnodt = 0.0;
     satrec.domdt = 0.0; satrec.e3    = 0.0; satrec.ee2   = 0.0;
     satrec.peo   = 0.0; satrec.pgho  = 0.0; satrec.pho   = 0.0;
     satrec.pinco = 0.0; satrec.plo   = 0.0; satrec.se2   = 0.0;
     satrec.se3   = 0.0; satrec.sgh2  = 0.0; satrec.sgh3  = 0.0;
     satrec.sgh4  = 0.0; satrec.sh2   = 0.0; satrec.sh3   = 0.0;
     satrec.si2   = 0.0; satrec.si3   = 0.0; satrec.sl2   = 0.0;
     satrec.sl3   = 0.0; satrec.sl4   = 0.0; satrec.gsto  = 0.0;
     satrec.xfact = 0.0; satrec.xgh2  = 0.0; satrec.xgh3  = 0.0;
     satrec.xgh4  = 0.0; satrec.xh2   = 0.0; satrec.xh3   = 0.0;
     satrec.xi2   = 0.0; satrec.xi3   = 0.0; satrec.xl2   = 0.0;
     satrec.xl3   = 0.0; satrec.xl4   = 0.0; satrec.xlamo = 0.0;
     satrec.zmol  = 0.0; satrec.zmos  = 0.0; satrec.atime = 0.0;
     satrec.xli   = 0.0; satrec.xni   = 0.0;

     #  ------------------------ earth constants -----------------------
     #  sgp4fix identify constants and allow alternate values
     #  this is now the only call for the constants
     (satrec.tumin, satrec.mu, satrec.radiusearthkm, satrec.xke,
       satrec.j2, satrec.j3, satrec.j4, satrec.j3oj2) = whichconst;

 	 # -------------------------------------------------------------------------

     satrec.error = 0;
     satrec.operationmode = opsmode;
     satrec.satnum = satn;

     """
     // sgp4fix - note the following variables are also passed directly via satrec.
     // it is possible to streamline the sgp4init call by deleting the "x"
     // variables, but the user would need to set the satrec.* values first. we
     // include the additional assignments in case twoline2rv is not used.
     """
     satrec.bstar   = xbstar;
    # sgp4fix allow additional parameters in the struct
     satrec.ndot    = xndot;
     satrec.nddot   = xnddot;
     satrec.ecco    = xecco;
     satrec.argpo   = xargpo;
     satrec.inclo   = xinclo;
     satrec.mo	    = xmo;
	# sgp4fix rename variables to clarify which mean motion is intended
     satrec.no_kozai= xno_kozai;
     satrec.nodeo   = xnodeo;

    # single averaged mean elements
     satrec.am = 0.0
     satrec.em = 0.0
     satrec.im = 0.0
     satrec.Om = 0.0
     satrec.mm = 0.0
     satrec.nm = 0.0

    # ------------------------ earth constants ----------------------- */
	# sgp4fix identify constants and allow alternate values no longer needed
	# getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );
     ss     = 78.0 / satrec.radiusearthkm + 1.0;
     #  sgp4fix use multiply for speed instead of pow
     qzms2ttemp = (120.0 - 78.0) / satrec.radiusearthkm;
     qzms2t = qzms2ttemp * qzms2ttemp * qzms2ttemp * qzms2ttemp;
     x2o3   =  2.0 / 3.0;

     satrec.init = 'y';
     satrec.t	 = 0.0;

    # sgp4fix remove satn as it is not needed in initl
     (
       satrec.no_unkozai,
       method,
       ainv,  ao,    satrec.con41,  con42, cosio,
       cosio2,eccsq, omeosq, posq,
       rp,    rteosq,sinio , satrec.gsto,
       ) = _initl(
           satrec.xke, satrec.j2, satrec.ecco, epoch, satrec.inclo, satrec.no_kozai, satrec.method,
           satrec.operationmode
         );
     satrec.a    = pow( satrec.no_unkozai*satrec.tumin , (-2.0/3.0) );
     satrec.alta = satrec.a*(1.0 + satrec.ecco) - 1.0;
     satrec.altp = satrec.a*(1.0 - satrec.ecco) - 1.0;

     """
     // sgp4fix remove this check as it is unnecessary
     // the mrt check in sgp4 handles decaying satellite cases even if the starting
     // condition is below the surface of te earth
     //     if (rp < 1.0)
     //       {
     //         printf("# *** satn%d epoch elts sub-orbital ***\n", satn);
     //         satrec.error = 5;
     //       }
     """

     if omeosq >= 0.0 or satrec.no_unkozai >= 0.0:

         satrec.isimp = 0;
         if rp < 220.0 / satrec.radiusearthkm + 1.0:
             satrec.isimp = 1;
         sfour  = ss;
         qzms24 = qzms2t;
         perige = (rp - 1.0) * satrec.radiusearthkm;

         #  - for perigees below 156 km, s and qoms2t are altered -
         if perige < 156.0:

             sfour = perige - 78.0;
             if perige < 98.0:
                 sfour = 20.0;
             #  sgp4fix use multiply for speed instead of pow
             qzms24temp =  (120.0 - sfour) / satrec.radiusearthkm;
             qzms24 = qzms24temp * qzms24temp * qzms24temp * qzms24temp;
             sfour  = sfour / satrec.radiusearthkm + 1.0;

         pinvsq = 1.0 / posq;

         tsi  = 1.0 / (ao - sfour);
         satrec.eta  = ao * satrec.ecco * tsi;
         etasq = satrec.eta * satrec.eta;
         eeta  = satrec.ecco * satrec.eta;
         psisq = fabs(1.0 - etasq);
         coef  = qzms24 * pow(tsi, 4.0);
         coef1 = coef / pow(psisq, 3.5);
         cc2   = coef1 * satrec.no_unkozai * (ao * (1.0 + 1.5 * etasq + eeta *
                        (4.0 + etasq)) + 0.375 * satrec.j2 * tsi / psisq * satrec.con41 *
                        (8.0 + 3.0 * etasq * (8.0 + etasq)));
         satrec.cc1   = satrec.bstar * cc2;
         cc3   = 0.0;
         if satrec.ecco > 1.0e-4:
             cc3 = -2.0 * coef * tsi * satrec.j3oj2 * satrec.no_unkozai * sinio / satrec.ecco;
         satrec.x1mth2 = 1.0 - cosio2;
         satrec.cc4    = 2.0* satrec.no_unkozai * coef1 * ao * omeosq * \
                           (satrec.eta * (2.0 + 0.5 * etasq) + satrec.ecco *
                           (0.5 + 2.0 * etasq) - satrec.j2 * tsi / (ao * psisq) *
                           (-3.0 * satrec.con41 * (1.0 - 2.0 * eeta + etasq *
                           (1.5 - 0.5 * eeta)) + 0.75 * satrec.x1mth2 *
                           (2.0 * etasq - eeta * (1.0 + etasq)) * cos(2.0 * satrec.argpo)));
         satrec.cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 *
                        (etasq + eeta) + eeta * etasq);
         cosio4 = cosio2 * cosio2;
         temp1  = 1.5 * satrec.j2 * pinvsq * satrec.no_unkozai;
         temp2  = 0.5 * temp1 * satrec.j2 * pinvsq;
         temp3  = -0.46875 * satrec.j4 * pinvsq * pinvsq * satrec.no_unkozai;
         satrec.mdot     = satrec.no_unkozai + 0.5 * temp1 * rteosq * satrec.con41 + 0.0625 * \
                            temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4);
         satrec.argpdot  = (-0.5 * temp1 * con42 + 0.0625 * temp2 *
                             (7.0 - 114.0 * cosio2 + 395.0 * cosio4) +
                             temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4));
         xhdot1            = -temp1 * cosio;
         satrec.nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) +
                              2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio;
         xpidot            =  satrec.argpdot+ satrec.nodedot;
         satrec.omgcof   = satrec.bstar * cc3 * cos(satrec.argpo);
         satrec.xmcof    = 0.0;
         if satrec.ecco > 1.0e-4:
             satrec.xmcof = -x2o3 * coef * satrec.bstar / eeta;
         satrec.nodecf = 3.5 * omeosq * xhdot1 * satrec.cc1;
         satrec.t2cof   = 1.5 * satrec.cc1;
         #  sgp4fix for divide by zero with xinco = 180 deg
         if fabs(cosio+1.0) > 1.5e-12:
             satrec.xlcof = -0.25 * satrec.j3oj2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio);
         else:
             satrec.xlcof = -0.25 * satrec.j3oj2 * sinio * (3.0 + 5.0 * cosio) / temp4;
         satrec.aycof   = -0.5 * satrec.j3oj2 * sinio;
         #  sgp4fix use multiply for speed instead of pow
         delmotemp = 1.0 + satrec.eta * cos(satrec.mo);
         satrec.delmo   = delmotemp * delmotemp * delmotemp;
         satrec.sinmao  = sin(satrec.mo);
         satrec.x7thm1  = 7.0 * cosio2 - 1.0;

         #  --------------- deep space initialization -------------
         if 2*pi / satrec.no_unkozai >= 225.0:

             satrec.method = 'd';
             satrec.isimp  = 1;
             tc    =  0.0;
             inclm = satrec.inclo;

             (
                 snodm, cnodm, sinim,  cosim, sinomm,
                 cosomm,day,   satrec.e3,     satrec.ee2,   em,
                 emsq,  gam,   satrec.peo,    satrec.pgho,  satrec.pho,
                 satrec.pinco, satrec.plo,   rtemsq, satrec.se2,   satrec.se3,
                 satrec.sgh2, satrec.sgh3, satrec.sgh4, satrec.sh2, satrec.sh3,
                 satrec.si2,  satrec.si3,  satrec.sl2, satrec.sl3, satrec.sl4,
                 s1,    s2,    s3,     s4,    s5,
                 s6,    s7,    ss1,    ss2,   ss3,
                 ss4,   ss5,   ss6,    ss7,   sz1,
                 sz2,   sz3,   sz11,   sz12,  sz13,
                 sz21,  sz22,  sz23,   sz31,  sz32,
                 sz33,  satrec.xgh2, satrec.xgh3, satrec.xgh4, satrec.xh2,
                 satrec.xh3, satrec.xi2, satrec.xi3, satrec.xl2, satrec.xl3,
                 satrec.xl4,   nm,    z1,     z2,    z3,
                 z11,   z12,   z13,    z21,   z22,
                 z23,   z31,   z32,    z33,   satrec.zmol,
                 satrec.zmos
             ) = _dscom(
                   epoch, satrec.ecco, satrec.argpo, tc, satrec.inclo, satrec.nodeo,
                   satrec.no_unkozai,
                   satrec.e3, satrec.ee2,
                   satrec.peo,  satrec.pgho,   satrec.pho, satrec.pinco,
                   satrec.plo,        satrec.se2, satrec.se3,
                   satrec.sgh2, satrec.sgh3,   satrec.sgh4,
                   satrec.sh2,  satrec.sh3,    satrec.si2, satrec.si3,
                   satrec.sl2,  satrec.sl3,    satrec.sl4,
                   satrec.xgh2, satrec.xgh3,   satrec.xgh4, satrec.xh2,
                   satrec.xh3,  satrec.xi2,    satrec.xi3,  satrec.xl2,
                   satrec.xl3,  satrec.xl4,
                   satrec.zmol, satrec.zmos
                 );
             (satrec.ecco, satrec.inclo, satrec.nodeo, satrec.argpo, satrec.mo
              ) = _dpper(
                   satrec, inclm, satrec.init,
                   satrec.ecco, satrec.inclo, satrec.nodeo, satrec.argpo, satrec.mo,
                   satrec.operationmode
                 );

             argpm  = 0.0;
             nodem  = 0.0;
             mm     = 0.0;

             (
                 em,    argpm,  inclm, mm,
                 nm,    nodem,
                 satrec.irez,  satrec.atime,
                 satrec.d2201, satrec.d2211, satrec.d3210, satrec.d3222,
                 satrec.d4410, satrec.d4422, satrec.d5220, satrec.d5232,
                 satrec.d5421, satrec.d5433, satrec.dedt,  satrec.didt,
                 satrec.dmdt, dndt,  satrec.dnodt, satrec.domdt,
                 satrec.del1,  satrec.del2,  satrec.del3,  satrec.xfact,
                 satrec.xlamo, satrec.xli,   satrec.xni
             ) = _dsinit(
                   satrec.xke,
                   cosim, emsq, satrec.argpo, s1, s2, s3, s4, s5, sinim, ss1, ss2, ss3, ss4,
                   ss5, sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33, satrec.t, tc,
                   satrec.gsto, satrec.mo, satrec.mdot, satrec.no_unkozai, satrec.nodeo,
                   satrec.nodedot, xpidot, z1, z3, z11, z13, z21, z23, z31, z33,
                   satrec.ecco, eccsq, em, argpm, inclm, mm, nm, nodem,
                   satrec.irez,  satrec.atime,
                   satrec.d2201, satrec.d2211, satrec.d3210, satrec.d3222 ,
                   satrec.d4410, satrec.d4422, satrec.d5220, satrec.d5232,
                   satrec.d5421, satrec.d5433, satrec.dedt,  satrec.didt,
                   satrec.dmdt,  satrec.dnodt, satrec.domdt,
                   satrec.del1,  satrec.del2,  satrec.del3,  satrec.xfact,
                   satrec.xlamo, satrec.xli,   satrec.xni
                 );

         #----------- set variables if not deep space -----------
         if satrec.isimp != 1:

           cc1sq          = satrec.cc1 * satrec.cc1;
           satrec.d2    = 4.0 * ao * tsi * cc1sq;
           temp           = satrec.d2 * tsi * satrec.cc1 / 3.0;
           satrec.d3    = (17.0 * ao + sfour) * temp;
           satrec.d4    = 0.5 * temp * ao * tsi * (221.0 * ao + 31.0 * sfour) * \
                            satrec.cc1;
           satrec.t3cof = satrec.d2 + 2.0 * cc1sq;
           satrec.t4cof = 0.25 * (3.0 * satrec.d3 + satrec.cc1 *
                            (12.0 * satrec.d2 + 10.0 * cc1sq));
           satrec.t5cof = 0.2 * (3.0 * satrec.d4 +
                            12.0 * satrec.cc1 * satrec.d3 +
                            6.0 * satrec.d2 * satrec.d2 +
                            15.0 * cc1sq * (2.0 * satrec.d2 + cc1sq));

     """
       /* finally propogate to zero epoch to initialize all others. */
       // sgp4fix take out check to let satellites process until they are actually below earth surface
//       if(satrec.error == 0)
     """
     sgp4(satrec, 0.0, whichconst);

     satrec.init = 'n';

     # sgp4fix return boolean. satrec.error contains any error codes
     return true;


"""
/*-----------------------------------------------------------------------------
*
*                             procedure sgp4
*
*  this procedure is the sgp4 prediction model from space command. this is an
*    updated and combined version of sgp4 and sdp4, which were originally
*    published separately in spacetrack report #3. this version follows the
*    methodology from the aiaa paper (2006) describing the history and
*    development of the code.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    satrec	 - initialised structure from sgp4init() call.
*    tsince	 - time eince epoch (minutes)
*
*  outputs       :
*    r           - position vector                     km
*    v           - velocity                            km/sec
*  return code - non-zero on error.
*                   1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
*                   2 - mean motion less than 0.0
*                   3 - pert elements, ecc < 0.0  or  ecc > 1.0
*                   4 - semi-latus rectum < 0.0
*                   5 - epoch elements are sub-orbital
*                   6 - satellite has decayed
*
*  locals        :
*    am          -
*    axnl, aynl        -
*    betal       -
*    cosim   , sinim   , cosomm  , sinomm  , cnod    , snod    , cos2u   ,
*    sin2u   , coseo1  , sineo1  , cosi    , sini    , cosip   , sinip   ,
*    cosisq  , cossu   , sinsu   , cosu    , sinu
*    delm        -
*    delomg      -
*    dndt        -
*    eccm        -
*    emsq        -
*    ecose       -
*    el2         -
*    eo1         -
*    eccp        -
*    esine       -
*    argpm       -
*    argpp       -
*    omgadf      -
*    pl          -
*    r           -
*    rtemsq      -
*    rdotl       -
*    rl          -
*    rvdot       -
*    rvdotl      -
*    su          -
*    t2  , t3   , t4    , tc
*    tem5, temp , temp1 , temp2  , tempa  , tempe  , templ
*    u   , ux   , uy    , uz     , vx     , vy     , vz
*    inclm       - inclination
*    mm          - mean anomaly
*    nm          - mean motion
*    nodem       - right asc of ascending node
*    xinc        -
*    xincp       -
*    xl          -
*    xlm         -
*    mp          -
*    xmdf        -
*    xmx         -
*    xmy         -
*    nodedf      -
*    xnode       -
*    nodep       -
*    np          -
*
*  coupling      :
*    getgravconst-
*    dpper
*    dpspace
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/
"""

def sgp4(satrec, tsince, whichconst=None):

     mrt = 0.0

     """
     /* ------------------ set mathematical constants --------------- */
     // sgp4fix divisor for divide by zero check on inclination
     // the old check used 1.0 + cos(pi-1.0e-9), but then compared it to
     // 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
     """
     temp4 =   1.5e-12;
     twopi = 2.0 * pi;
     x2o3  = 2.0 / 3.0;
     #  sgp4fix identify constants and allow alternate values
     # tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 = whichconst
     vkmpersec     = satrec.radiusearthkm * satrec.xke/60.0;

     #  --------------------- clear sgp4 error flag -----------------
     satrec.t     = tsince;
     satrec.error = 0;
     satrec.error_message = None

     #  ------- update for secular gravity and atmospheric drag -----
     xmdf    = satrec.mo + satrec.mdot * satrec.t;
     argpdf  = satrec.argpo + satrec.argpdot * satrec.t;
     nodedf  = satrec.nodeo + satrec.nodedot * satrec.t;
     argpm   = argpdf;
     mm      = xmdf;
     t2      = satrec.t * satrec.t;
     nodem   = nodedf + satrec.nodecf * t2;
     tempa   = 1.0 - satrec.cc1 * satrec.t;
     tempe   = satrec.bstar * satrec.cc4 * satrec.t;
     templ   = satrec.t2cof * t2;

     if satrec.isimp != 1:

         delomg = satrec.omgcof * satrec.t;
         #  sgp4fix use mutliply for speed instead of pow
         delmtemp =  1.0 + satrec.eta * cos(xmdf);
         delm   = satrec.xmcof * \
                  (delmtemp * delmtemp * delmtemp -
                  satrec.delmo);
         temp   = delomg + delm;
         mm     = xmdf + temp;
         argpm  = argpdf - temp;
         t3     = t2 * satrec.t;
         t4     = t3 * satrec.t;
         tempa  = tempa - satrec.d2 * t2 - satrec.d3 * t3 - \
                          satrec.d4 * t4;
         tempe  = tempe + satrec.bstar * satrec.cc5 * (sin(mm) -
                          satrec.sinmao);
         templ  = templ + satrec.t3cof * t3 + t4 * (satrec.t4cof +
                          satrec.t * satrec.t5cof);

     nm    = satrec.no_unkozai;
     em    = satrec.ecco;
     inclm = satrec.inclo;
     if satrec.method == 'd':

         tc = satrec.t;
         (
             atime, em,    argpm,  inclm, xli,
             mm,    xni,   nodem,  dndt,  nm,
         ) = _dspace(
               satrec.irez,
               satrec.d2201, satrec.d2211, satrec.d3210,
               satrec.d3222, satrec.d4410, satrec.d4422,
               satrec.d5220, satrec.d5232, satrec.d5421,
               satrec.d5433, satrec.dedt,  satrec.del1,
               satrec.del2,  satrec.del3,  satrec.didt,
               satrec.dmdt,  satrec.dnodt, satrec.domdt,
               satrec.argpo, satrec.argpdot, satrec.t, tc,
               satrec.gsto, satrec.xfact, satrec.xlamo,
               satrec.no_unkozai, satrec.atime,
               em, argpm, inclm, satrec.xli, mm, satrec.xni,
               nodem, nm
             );

     if nm <= 0.0:

         satrec.error_message = ('mean motion {0:f} is less than zero'
                                 .format(nm))
         satrec.error = 2;
         #  sgp4fix add return
         return false, false;

     am = pow((satrec.xke / nm),x2o3) * tempa * tempa;
     nm = satrec.xke / pow(am, 1.5);
     em = em - tempe;

     #  fix tolerance for error recognition
     #  sgp4fix am is fixed from the previous nm check
     if em >= 1.0 or em < -0.001:  # || (am < 0.95)

         satrec.error_message = ('mean eccentricity {0:f} not within'
                                 ' range 0.0 <= e < 1.0'.format(em))
         satrec.error = 1;
         #  sgp4fix to return if there is an error in eccentricity
         return false, false;

     #  sgp4fix fix tolerance to avoid a divide by zero
     if em < 1.0e-6:
         em  = 1.0e-6;
     mm     = mm + satrec.no_unkozai * templ;
     xlm    = mm + argpm + nodem;
     emsq   = em * em;
     temp   = 1.0 - emsq;

     nodem  = nodem % twopi if nodem >= 0.0 else -(-nodem % twopi)
     argpm  = argpm % twopi
     xlm    = xlm % twopi
     mm     = (xlm - argpm - nodem) % twopi

     # sgp4fix recover singly averaged mean elements
     satrec.am = am;
     satrec.em = em;
     satrec.im = inclm;
     satrec.Om = nodem;
     satrec.om = argpm;
     satrec.mm = mm;
     satrec.nm = nm;

     #  ----------------- compute extra mean quantities -------------
     sinim = sin(inclm);
     cosim = cos(inclm);

     #  -------------------- add lunar-solar periodics --------------
     ep     = em;
     xincp  = inclm;
     argpp  = argpm;
     nodep  = nodem;
     mp     = mm;
     sinip  = sinim;
     cosip  = cosim;
     if satrec.method == 'd':

         ep, xincp, nodep, argpp, mp = _dpper(
               satrec, satrec.inclo,
               'n', ep, xincp, nodep, argpp, mp, satrec.operationmode
             );
         if xincp < 0.0:

             xincp  = -xincp;
             nodep = nodep + pi;
             argpp  = argpp - pi;

         if ep < 0.0 or ep > 1.0:

             satrec.error_message = ('perturbed eccentricity {0:f} not within'
                                     ' range 0.0 <= e <= 1.0'.format(ep))
             satrec.error = 3;
             #  sgp4fix add return
             return false, false;

     #  -------------------- long period periodics ------------------
     if satrec.method == 'd':

         sinip =  sin(xincp);
         cosip =  cos(xincp);
         satrec.aycof = -0.5*satrec.j3oj2*sinip;
         #  sgp4fix for divide by zero for xincp = 180 deg
         if fabs(cosip+1.0) > 1.5e-12:
             satrec.xlcof = -0.25 * satrec.j3oj2 * sinip * (3.0 + 5.0 * cosip) / (1.0 + cosip);
         else:
             satrec.xlcof = -0.25 * satrec.j3oj2 * sinip * (3.0 + 5.0 * cosip) / temp4;

     axnl = ep * cos(argpp);
     temp = 1.0 / (am * (1.0 - ep * ep));
     aynl = ep* sin(argpp) + temp * satrec.aycof;
     xl   = mp + argpp + nodep + temp * satrec.xlcof * axnl;

     #  --------------------- solve kepler's equation ---------------
     u    = (xl - nodep) % twopi
     eo1  = u;
     tem5 = 9999.9;
     ktr = 1;
     #    sgp4fix for kepler iteration
     #    the following iteration needs better limits on corrections
     while fabs(tem5) >= 1.0e-12 and ktr <= 10:

         sineo1 = sin(eo1);
         coseo1 = cos(eo1);
         tem5   = 1.0 - coseo1 * axnl - sineo1 * aynl;
         tem5   = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
         if fabs(tem5) >= 0.95:
             tem5 = 0.95 if tem5 > 0.0 else -0.95;
         eo1    = eo1 + tem5;
         ktr = ktr + 1;

     #  ------------- short period preliminary quantities -----------
     ecose = axnl*coseo1 + aynl*sineo1;
     esine = axnl*sineo1 - aynl*coseo1;
     el2   = axnl*axnl + aynl*aynl;
     pl    = am*(1.0-el2);
     if pl < 0.0:

         satrec.error_message = ('semilatus rectum {0:f} is less than zero'
                                 .format(pl))
         satrec.error = 4;
         #  sgp4fix add return
         return false, false;

     else:

         rl     = am * (1.0 - ecose);
         rdotl  = sqrt(am) * esine/rl;
         rvdotl = sqrt(pl) / rl;
         betal  = sqrt(1.0 - el2);
         temp   = esine / (1.0 + betal);
         sinu   = am / rl * (sineo1 - aynl - axnl * temp);
         cosu   = am / rl * (coseo1 - axnl + aynl * temp);
         su     = atan2(sinu, cosu);
         sin2u  = (cosu + cosu) * sinu;
         cos2u  = 1.0 - 2.0 * sinu * sinu;
         temp   = 1.0 / pl;
         temp1  = 0.5 * satrec.j2 * temp;
         temp2  = temp1 * temp;

         #  -------------- update for short period periodics ------------
         if satrec.method == 'd':

             cosisq                 = cosip * cosip;
             satrec.con41  = 3.0*cosisq - 1.0;
             satrec.x1mth2 = 1.0 - cosisq;
             satrec.x7thm1 = 7.0*cosisq - 1.0;

         mrt   = rl * (1.0 - 1.5 * temp2 * betal * satrec.con41) + \
                 0.5 * temp1 * satrec.x1mth2 * cos2u;
         su    = su - 0.25 * temp2 * satrec.x7thm1 * sin2u;
         xnode = nodep + 1.5 * temp2 * cosip * sin2u;
         xinc  = xincp + 1.5 * temp2 * cosip * sinip * cos2u;
         mvt   = rdotl - nm * temp1 * satrec.x1mth2 * sin2u / satrec.xke;
         rvdot = rvdotl + nm * temp1 * (satrec.x1mth2 * cos2u +
                 1.5 * satrec.con41) / satrec.xke;

         #  --------------------- orientation vectors -------------------
         sinsu =  sin(su);
         cossu =  cos(su);
         snod  =  sin(xnode);
         cnod  =  cos(xnode);
         sini  =  sin(xinc);
         cosi  =  cos(xinc);
         xmx   = -snod * cosi;
         xmy   =  cnod * cosi;
         ux    =  xmx * sinsu + cnod * cossu;
         uy    =  xmy * sinsu + snod * cossu;
         uz    =  sini * sinsu;
         vx    =  xmx * cossu - cnod * sinsu;
         vy    =  xmy * cossu - snod * sinsu;
         vz    =  sini * cossu;

         #  --------- position and velocity (in km and km/sec) ----------
         _mr = mrt * satrec.radiusearthkm
         r = (_mr * ux, _mr * uy, _mr * uz)
         v = ((mvt * ux + rvdot * vx) * vkmpersec,
              (mvt * uy + rvdot * vy) * vkmpersec,
              (mvt * uz + rvdot * vz) * vkmpersec)

     #  sgp4fix for decaying satellites
     if mrt < 1.0:

         satrec.error_message = ('mrt {0:f} is less than 1.0 indicating'
                                 ' the satellite has decayed'.format(mrt))
         satrec.error = 6;

     return r, v;

"""
/* -----------------------------------------------------------------------------
*
*                           function gstime
*
*  this function finds the greenwich sidereal time.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    jdut1       - julian date in ut1             days from 4713 bc
*
*  outputs       :
*    gstime      - greenwich sidereal time        0 to 2pi rad
*
*  locals        :
*    temp        - temporary variable for doubles   rad
*    tut1        - julian centuries from the
*                  jan 1, 2000 12 h epoch (ut1)
*
*  coupling      :
*    none
*
*  references    :
*    vallado       2004, 191, eq 3-45
* --------------------------------------------------------------------------- */
"""

def gstime(jdut1):

     tut1 = (jdut1 - 2451545.0) / 36525.0;
     temp = -6.2e-6* tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1 + \
             (876600.0*3600 + 8640184.812866) * tut1 + 67310.54841;  #  sec
     temp = (temp * deg2rad / 240.0) % twopi # 360/86400 = 1/240, to deg, to rad

     #  ------------------------ check quadrants ---------------------
     if temp < 0.0:
         temp += twopi;

     return temp;

# The routine was originally marked private, so make it available under
# the old name for compatibility:
_gstime = gstime

"""
/* -----------------------------------------------------------------------------
*
*                           function getgravconst
*
*  this function gets constants for the propagator. note that mu is identified to
*    facilitiate comparisons with newer models. the common useage is wgs72.
*
*  author        : david vallado                  719-573-2600   21 jul 2006
*
*  inputs        :
*    whichconst  - which set of constants to use  wgs72old, wgs72, wgs84
*
*  outputs       :
*    tumin       - minutes in one time unit
*    mu          - earth gravitational parameter
*    radiusearthkm - radius of the earth in km
*    xke         - reciprocal of tumin
*    j2, j3, j4  - un-normalized zonal harmonic values
*    j3oj2       - j3 divided by j2
*
*  locals        :
*
*  coupling      :
*    none
*
*  references    :
*    norad spacetrack report #3
*    vallado, crawford, hujsak, kelso  2006
  --------------------------------------------------------------------------- */
"""

def getgravconst(whichconst):

       if whichconst == 'wgs72old':
           mu     = 398600.79964;        #  in km3 / s2
           radiusearthkm = 6378.135;     #  km
           xke    = 0.0743669161;
           tumin  = 1.0 / xke;
           j2     =   0.001082616;
           j3     =  -0.00000253881;
           j4     =  -0.00000165597;
           j3oj2  =  j3 / j2;

           #  ------------ wgs-72 constants ------------
       elif whichconst == 'wgs72':
           mu     = 398600.8;            #  in km3 / s2
           radiusearthkm = 6378.135;     #  km
           xke    = 60.0 / sqrt(radiusearthkm*radiusearthkm*radiusearthkm/mu);
           tumin  = 1.0 / xke;
           j2     =   0.001082616;
           j3     =  -0.00000253881;
           j4     =  -0.00000165597;
           j3oj2  =  j3 / j2;

       elif whichconst == 'wgs84':
           #  ------------ wgs-84 constants ------------
           mu     = 398600.5;            #  in km3 / s2
           radiusearthkm = 6378.137;     #  km
           xke    = 60.0 / sqrt(radiusearthkm*radiusearthkm*radiusearthkm/mu);
           tumin  = 1.0 / xke;
           j2     =   0.00108262998905;
           j3     =  -0.00000253215306;
           j4     =  -0.00000161098761;
           j3oj2  =  j3 / j2;

       return tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2
