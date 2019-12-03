import sys
del sys.path[0]
import sgp4
import sgp4.vallado_cpp
print(dir(sgp4.vallado_cpp))

# s = sgp4.vallado_cpp.Satrec()
# print(s)
# print("len", len(s))

line1 = '1 25544U 98067A   14020.93268519  .00009878  00000-0  18200-3 0  5082'
line2 = '2 25544  51.6498 109.4756 0003572  55.9686 274.8005 15.49815350868473'

s = sgp4.vallado_cpp.Satrec.twoline2rv(line1, line2)

print(s)
print(s.satnum)
print(s.method)

print(s.sgp4(0))
print(s.sgp4(1))

a = sgp4.vallado_cpp.SatrecArray([s, s])
print("len", len(a))

from numpy import array, ndarray
t = array([0.0, 1, 2, 3, 49999999])
r = ndarray((len(a), len(t), 3))
v = ndarray((len(a), len(t), 3))
e = ndarray((len(a), len(t)), 'uint8')
print(v)
print(a.sgp4(t, r, v, e))
print(r)
print(v)
print(e)
