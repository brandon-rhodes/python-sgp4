"""Public API that tries to import C++ module, but falls back to Python."""

__all__ = 'Satrec', 'SatrecArray', 'jday'

from .ext import jday2 as jday

try:
    from .vallado_cpp import Satrec, SatrecArray as _SatrecArray
except ImportError:
    accelerated = False
    from .model import Satrec, SatrecArray
else:
    accelerated = True
    from .wrapper import SatrecArray
