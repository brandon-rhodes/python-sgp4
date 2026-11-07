#!/bin/bash

set -e
cd "$(readlink -f $(dirname "${BASH_SOURCE[0]}"))"

export PYTHONDONTWRITEBYTECODE=dont

echo
echo Building source distribution
echo

python3 setup.py sdist
SDIST_TAR_GZ="$(ls -t dist | head -1)"
echo "$SDIST_TAR_GZ"

# Change directory to the mostly-empty `ci/` directory, so Python can't
# cheat and import 'sgp4' directly from the source tree.

cd ci/

echo
echo Testing pure-Python version
echo

PYTHON_SGP4_COMPILE=never uv run \
   --isolated \
   --no-project \
   --reinstall-package sgp4 \
   --with ../dist/$SDIST_TAR_GZ,numpy \
   -m unittest sgp4.tests

echo
echo Testing compiled extension
echo

# Okay, this is crazy, and just burned about an hour of my time: despite
# the pile of options I'm using here, plus "PYTHON_SGP4_COMPILE" being
# listed in "pyproject.toml" as a cache key, the ONLY way I've managed
# to get "uv" to rebuild "sgp4" for this second run is to touch the
# tar.gz file to update its mtime.  Gads.

touch ../dist/$SDIST_TAR_GZ

PYTHON_SGP4_COMPILE=always uv run \
   --isolated \
   --no-project \
   --reinstall-package sgp4 \
   --with ../dist/$SDIST_TAR_GZ,numpy \
   -m unittest sgp4.tests

echo
echo Running doctests
echo

PYTHON_SGP4_COMPILE=always uv run \
   --isolated \
   --no-project \
   --with ../dist/$SDIST_TAR_GZ,numpy \
   -m doctest -o ELLIPSIS ../README.rst

# Since "doctest" itself is silent upon success:

echo Doctests successful!
echo
