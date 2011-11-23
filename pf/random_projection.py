#!/usr/bin/python

from numpy import *
from pylab import *
from numpy.random import *

n = 100
gamma = 0.1
d = 300
k = int(log(n)/gamma**2)
print "k: ", k

sq3 = sqrt(3/float(k))
R = zeros( (k, d), float)
for i in range(k):
    for j in range(d):
        t = randint(1, 7)       # 7 exclusive
        if (t == 1):
            R[i, j] = sq3
        elif (t == 6):
            R[i, j] = -sq3

for c in range(n):
    x  = [0 for i in range(d)]
    dist = 0 
    for i in range(d):
        x[i] = random()
        dist += x[i]**2
    dist = sqrt(dist)
    
    xn = dot(R, x)
    distn = 0
    for i in range(k):
        distn += xn[i]**2
    distn = sqrt(distn)

    print "dist: ", dist, " distn: ", distn

