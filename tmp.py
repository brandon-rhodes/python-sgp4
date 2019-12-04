# ----------------------------------------------------
# Based on the README at https://github.com/bwinkel/cysgp4

import numpy as np
from cysgp4 import PyTle, PyObserver, propagate_many

with open('science.txt') as f:
    text = f.read()

all_lines = text.splitlines()

# Need to convert them to a list of tuples (each tuple consisting
# of the three TLE strings)
tle_list = list(zip(*tuple(
    all_lines[idx::3] for idx in range(3)
    )))
# Create an array of PyTle and PyObserver objects, and MJDs
tles = np.array([
    PyTle(*tle) for tle in tle_list
    ])[np.newaxis, np.newaxis, :20]  # use first 20 TLEs
mjds = np.linspace(
    58805.5, 58806.5, 1000  # 1000 time steps
    )[:, np.newaxis, np.newaxis]

t0 = __import__('time').time()
result = propagate_many(
    mjds, tles, None,
    do_eci_pos=True, do_eci_vel=True, do_geo=False, do_topo=False
    )
print('TIME:', __import__('time').time() - t0)

print(result.keys())
print(result['eci_pos'].shape)
print(result['eci_pos'][0])

# -------------------------------------

import sys
del sys.path[0]  # prevent importing from "sgp4/" in this directory

from sgp4.vallado_cpp import Satrec, SatrecArray

sats = []
lines = iter(text.splitlines())
for name in lines:
    line1 = next(lines)
    line2 = next(lines)
    sat = Satrec.twoline2rv(line1, line2)
    # print(sat)
    # print(sat.satnum)
    # print(sat.method)
    sats.append(sat)

a = SatrecArray(sats[:20])
print("Number of satellites:", len(a))

mjd = np.linspace(58805.5, 58806.5, 1000)  # 1000 time steps

# t = array([0.0, 1, 2, 3, 49999999])
t = mjd
r = np.ndarray((len(a), len(t), 3))
v = np.ndarray((len(a), len(t), 3))
e = np.ndarray((len(a), len(t)), 'uint8')

t0 = __import__('time').time()
print(a.sgp4(t, r, v, e))
print(__import__('time').time() - t0)

print(r.shape)
print(r[0])
# print(v)
# print(e)
