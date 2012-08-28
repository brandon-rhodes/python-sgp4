#!/bin/bash

set -e
cd $(dirname "$0")

# Determine how much of the original C++ we have retained.

CPP=../AIAA-2006-6753/cpp/sgp4unit.cpp

if ! [ -f "$CPP" ]
then
    echo
    echo "Error: cannot find path $CPP"
    echo
    echo "Please unzip AIAA-2006-6753.zip in the parent directory to this one"
    echo
    exit 2
fi

echo
echo "The number of lines of C++ that have disappeared is:"
diff -b -U 9999 "$CPP" sgp4/propagation.py | grep '^-' | wc -l
echo
echo "The number of Python lines that have been added is:"
diff -b -U 9999 "$CPP" sgp4/propagation.py | grep '^\+' | wc -l
echo
echo "The number of C++ lines that are still preserved is:"
diff -b -U 9999 "$CPP" sgp4/propagation.py | grep '^ ' | wc -l
echo
