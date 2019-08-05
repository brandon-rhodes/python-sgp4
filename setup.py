from distutils.core import setup
from textwrap import dedent

cython_installed = False
try:
    from Cython.Build import cythonize
    cython_installed = True
except ImportError:
    pass

import sgp4, sgp4.model

description, long_description = sgp4.__doc__.split('\n', 1)
satdoc = dedent(sgp4.model.Satellite.__doc__.split('\n', 1)[1])
long_description = long_description.replace('entry.', 'entry.' + satdoc)

setup(name = 'sgp4',
      version = '1.4',
      description = description,
      long_description = long_description,
      license = 'MIT',
      author = 'Brandon Rhodes',
      author_email = 'brandon@rhodesmill.org',
      url = 'https://github.com/brandon-rhodes/python-sgp4',
      classifiers = [
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.2',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Topic :: Scientific/Engineering :: Astronomy',
        ],
      packages = ['sgp4'],
      ext_modules = cythonize("sgp4/cpropagation.pyx") if cython_installed else []
      )
