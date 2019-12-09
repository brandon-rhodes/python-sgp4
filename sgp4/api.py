
# TODO: turn Skyfield time into two JD floats.
# Hmm. That's complicated.

# jd = 58805.5 # TODO
# fr = 0.0

# sat = Satrec.twoline2rv(line1, line2)
# r, v, e = sat.sgp4(jd, fr)

# mjd = np.linspace(58805.5, 58806.5, 1000) # TODO
# fr = 0.0

# a = SatrecArray([sat, sat])
# r, v, e = sat.sgp4(jd, fr)

# --------------------

try:
    from .vallado_cpp import Satrec, SatrecArray
except ImportError:
    accelerated = False
    from .model import Satrec#, SatrecArray
else:
    accelerated = True
