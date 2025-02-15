import os
import sys
from distutils.core import setup, Extension

import sgp4
description, long_description = sgp4.__doc__.split('\n', 1)

# It is hard to write C extensions that support both Python 2 and 3, so
# we opt here to support the acceleration only for Python 3.

ext_modules = []
if sys.version_info[0] == 3:

    # This lets CI force us to exit with an error if compilation fails,
    # instead of falling back silently to the backup Python code.
    optional = True
    if os.environ.get('SGP4_FORCE_COMPILE') == 'true':
        optional = False

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

# Read the package's "__version__" without importing it.
path = 'sgp4/__init__.py'
with open(path, 'rb') as f:
    text = f.read().decode('utf-8')
text = text.replace('-*- coding: utf-8 -*-', '')  # for Python 2.7
namespace = {}
eval(compile(text, path, 'exec'), namespace)

setup(name = 'sgp4',
      version = namespace['__version__'],
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
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Topic :: Scientific/Engineering :: Astronomy',
      ],
      packages = ['sgp4'],
      package_data = {'sgp4': ['SGP4-VER.TLE', 'sample*', 'tcppver.out']},
      ext_modules = ext_modules,
)
