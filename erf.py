import math

for x in range(1, 20):
  y = x / 10.0
  print 0.5 * math.sqrt(math.pi / y) * math.erf(math.sqrt(y))
