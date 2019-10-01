from sgp4.io import twoline2rv
from sgp4.model import Satellite
from sgp4.earth_gravity import wgs72

class Satellite2(Satellite):
    """
    Subclass of Satellite that uses the twoline2rv function as a constructor.

    This is useful if you're subclassing Satellite and need to modify the constructor (e.g. attach another attribute)
    """
    def __init__(self, longstr1, longstr2, whichconst, afspc_mode=False):
        _Satellite = twoline2rv(longstr1, longstr2, whichconst, afspc_mode)
        for attr in vars(_Satellite):
            setattr(self, attr, vars(_Satellite)[attr])
