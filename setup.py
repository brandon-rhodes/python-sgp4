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
import sysconfig
from distutils.ccompiler import get_default_compiler
from distutils.core import setup, Extension

if sys.version_info[0] == 2:
    message = (
        'alas, this newer version of "sgp4" no longer supports Python 2,'
        ' but you can install "sgp4==2.27" instead'
    )
    raise RuntimeError(message)

PYTHON_SGP4_COMPILE = os.environ.get('PYTHON_SGP4_COMPILE', '')

USE_PY_LIMITED_API = (
    # Py_buffer is required
    sys.version_info >= (3, 11)
    # LIMITED_API is not compatible with free-threading (as of CPython 3.14)
    and not sysconfig.get_config_var("Py_GIL_DISABLED")
)
ABI3_TARGET_VERSION = "".join(str(_) for _ in sys.version_info[:2])
ABI3_TARGET_HEX = hex(sys.hexversion & 0xFFFF00F0)

ext_modules = []
define_macros = []

if USE_PY_LIMITED_API:
    define_macros.append(("Py_LIMITED_API", ABI3_TARGET_HEX))

if PYTHON_SGP4_COMPILE != 'never':

    if PYTHON_SGP4_COMPILE == 'always':
        optional = False
    else:
        optional = True

    # TODO: can we safely figure out how to use a pair of options
    # like these, adapted to as many platforms as possible, to use
    # multiple processors when available?
    # extra_compile_args=['-fopenmp'],
    # extra_link_args=['-fopenmp'],
    extra_compile_args = ['-ffloat-store']
    if get_default_compiler() == "msvc":
        extra_compile_args.append("/std:c++20")

    ext_modules.append(Extension(
        'sgp4.vallado_cpp',
        optional=optional,
        sources=[
            'extension/SGP4.cpp',
            'extension/wrapper.cpp',
        ],
        extra_compile_args=extra_compile_args,
        define_macros=define_macros,
        py_limited_api=USE_PY_LIMITED_API,
    ))

# Read the package's docstring and "__version__" without importing it.

path = 'sgp4/__init__.py'
with open(path, 'rb') as f:
    text = f.read().decode('utf-8')

namespace = {}
eval(compile(text, path, 'exec'), namespace)

version = namespace['__version__']
description = namespace['__doc__'].split('\n', 1)[0].strip()

setup(
    version = version,
    description = description,
    author = 'Brandon Rhodes',
    author_email = 'brandon@rhodesmill.org',
    url = 'https://github.com/brandon-rhodes/python-sgp4',
    classifiers = [
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
        'Programming Language :: Python :: 3.14',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
    packages = ['sgp4'],
    package_data = {'sgp4': ['SGP4-VER.TLE', 'sample*', 'tcppver.out']},
    provides = ['sgp4'],
    ext_modules = ext_modules,
    options={
        "bdist_wheel": {"py_limited_api": "cp%s" % ABI3_TARGET_VERSION}
    } if USE_PY_LIMITED_API else {},
)
