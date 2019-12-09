"""Public API that tries to import C++ module, but falls back to Python."""

# mjd = np.linspace(58805.5, 58806.5, 1000) # TODO
# fr = 0.0

# a = SatrecArray([sat, sat])
# r, v, e = sat.sgp4(jd, fr)

# --------------------

try:
    from .vallado_cpp import Satrec, SatrecArray, jday
except ImportError:
    accelerated = False
    from .ext import jday2 as jday
    from .model import Satrec#, SatrecArray
else:
    accelerated = True
