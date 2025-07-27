#!/bin/bash
#
# We 'cd' to 'ci/' first so that Python doesn't import the 'sgp4'
# package source code from the current directory.

set -e
export PYTHONDONTWRITEBYTECODE=dont

(
    set -e
    python2.7 setup.py install --prefix ci/legacy
    cd ci
    PYTHONPATH=legacy/lib/python2.7/site-packages python2.7 -m sgp4.tests
)

uv run \
   --directory ci \
   --with ..,numpy \
   -m unittest sgp4.tests
