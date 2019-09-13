from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules=cythonize(["propagation.pyx", "earth_gravity.pyx", "ext.pyx", "io.pyx", "model.pyx"]),
)
