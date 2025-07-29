# This setup.py supports three modes.
#
# 1. We normally perform a best-effort install: we will try to compile
#    the C++ module, but if it won't compile, then we fall back to
#    installing the pure-Python SGP4 implementation instead.
#
# 2. `PYTHON_SGP4_COMPILE=always` insists on compiling the accelerated
#    module written in C++.  Install fails if it can't be compiled.
#
# 3. `PYTHON_SGP4_COMPILE=never` insists on installing the pure-Python
#    version of SGP4.  Users probably never want this, but it's useful
#    to support local and CI tests that want to make sure the
#    pure-Python version is the one getting tested.

import os
import sys
from distutils.core import setup, Extension

PYTHON_SGP4_COMPILE = os.environ.get('PYTHON_SGP4_COMPILE', '')
ext_modules = []

if sys.version_info[0] == 3 and PYTHON_SGP4_COMPILE != 'never':

    if PYTHON_SGP4_COMPILE == 'always':
        optional = False
    else:
        optional = True

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

# Read the package's docstring and "__version__" without importing it.
path = 'sgp4/__init__.py'
with open(path, 'rb') as f:
    text = f.read().decode('utf-8')
text = text.replace('-*- coding: utf-8 -*-', '')  # for Python 2.7

namespace = {}
eval(compile(text, path, 'exec'), namespace)

version = namespace['__version__']
description, long_description = namespace['__doc__'].split('\n', 1)

# Avoid "`long_description` has syntax errors in markup" from extra
# whitespace that somehow creeps in.
description = description.strip()
long_description = long_description.strip() + '\n'

setup(
    name = 'sgp4',
    version = version,
    description = description,
    long_description = long_description,
    long_description_content_type = 'text/x-rst',
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
        'Programming Language :: Python :: 3.13',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
    packages = ['sgp4'],
    package_data = {'sgp4': ['SGP4-VER.TLE', 'sample*', 'tcppver.out']},
    provides = ['sgp4'],
    ext_modules = ext_modules,
)
