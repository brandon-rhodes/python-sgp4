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
print("len", len(s))
print(s.satnum)
print(s.method)

t0 = __import__('time').time()
print(s.sgp4(0))
print(__import__('time').time() - t0)

t0 = __import__('time').time()
print(s.sgp4(1))
print(__import__('time').time() - t0)

t0 = __import__('time').time()
print(s.sgp4(2))
print(__import__('time').time() - t0)

t0 = __import__('time').time()
print(s.sgp4(3))
print(__import__('time').time() - t0)

