"""Public API that tries to import C++ module, but falls back to Python."""

from .ext import jday2 as jday

try:
    from .vallado_cpp import Satrec, SatrecArray as _SatrecArray
except ImportError:
    accelerated = False
    from .model import Satrec, SatrecArray
else:
    accelerated = True

    class SatrecArray(_SatrecArray):
        def sgp4(self, jd, fr):
            jd = jd.astype('float64', copy=False)
            fr = fr.astype('float64', copy=False)

            ilength = len(self)
            jlength, = jd.shape
            eshape = ilength, jlength
            fshape = ilength, jlength, 3

            array = type(jd)
            e = array(eshape, 'uint8')
            r = array(fshape, 'float64')
            v = array(fshape, 'float64')
            self._sgp4(jd, fr, e, r, v)
            return e, r, v
