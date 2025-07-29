#!/bin/bash
#
# We 'cd' to 'ci/' first so that Python doesn't import the 'sgp4'
# package source code from the current directory.

set -e
export PYTHONDONTWRITEBYTECODE=dont

# Python 2.

(
    set -e
    python2.7 setup.py install --prefix ci/legacy
    cd ci
    PYTHONPATH=legacy/lib/python2.7/site-packages python2.7 -m sgp4.tests
)

# Python 3, using pure-Python SGP4 implementation.

PYTHON_SGP4_COMPILE=never uv run \
   --directory ci \
   --with ..,numpy \
   --reinstall-package sgp4 \
   -m unittest sgp4.tests

# Python 3, using the compiled C++ SGP4 implementation.

PYTHON_SGP4_COMPILE=always uv run \
   --directory ci \
   --with ..,numpy \
   --reinstall-package sgp4 \
   -m unittest sgp4.tests
