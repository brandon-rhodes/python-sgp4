#!/usr/bin/env sh

cd $SRC_DIR

if [[ "$PY_VER" = '2.6' ]]; then
  unit2 discover sgp4
else
  python -m unittest discover sgp4
fi
