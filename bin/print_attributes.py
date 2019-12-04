#!/usr/bin/env python

from __future__ import print_function

from fileinput import input
from sgp4.vallado_cpp import Satrec

def main():
    lines = iter(input())
    for line in lines:
        name = line
        line1 = next(lines)
        line2 = next(lines)
        sat = Satrec.twoline2rv(line1, line2)
        for name in dir(sat):
            if name.startswith('_') or name in ('sgp4', 'twoline2rv'):
                continue
            value = getattr(sat, name)
            print(name, value)
        print()

if __name__ == '__main__':
    try:
        main()
    except BrokenPipeError:
        pass


