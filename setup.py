import os
import re
import sys
from distutils.core import setup, Extension
from textwrap import dedent

import sgp4
description, long_description = sgp4.__doc__.split('\n', 1)

# Force compilation on Travis CI + Python 3 to make sure it keeps working.
optional = True
if sys.version_info[0] != 2 and os.environ.get('TRAVIS') == 'true':
    optional = False

# It is hard to write C extensions that support both Python 2 and 3, so
# we opt here to support the acceleration only for Python 3.
ext_modules = []
if sys.version_info[0] == 3:
    ext_modules.append(Extension(
        'sgp4.vallado_cpp',
        optional=optional,
        sources=[
            'extension/SGP4.cpp',
            'extension/wrapper.cpp',
        ],

        # TODO: can we safely figure out how to use a pair of options
        # like these, adapted to as many platforms as possible, to use
        # multiple processors when available?
        # extra_compile_args=['-fopenmp'],
        # extra_link_args=['-fopenmp'],
        extra_compile_args=['-ffloat-store'],
    ))

setup(name = 'sgp4',
      version = '2.19',
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
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering :: Astronomy',
        ],
      packages = ['sgp4'],
      package_data = {'sgp4': ['SGP4-VER.TLE', 'sample*', 'tcppver.out']},
      ext_modules = ext_modules,
)
