#!/usr/bin/python

from numpy import *
from pylab import *
from numpy.random import *

gamma = 1
d = 2000
k = int(log(d)/gamma**2) + 10
print "k: ", k

sq3 = sqrt(3)
R = zeros( (k, d), float)
for i in range(k):
    for j in range(d):
        t = randint(1, 7)       # 7 exclusive
        if (t == 1):
            R[i, j] = sq3
        elif (t == 6):
            R[i, j] = -sq3

for c in range(1):
    x  = [0 for i in range(d)]
    y  = [0 for i in range(d)]
    dist = 0 
    for i in range(d):
        x[i] = random()
        y[i] = random()
        dist += (x[i] - y[i])**2
    dist = sqrt(dist)
    
    xn = dot(R, x)
    yn = dot(R, y)
    distn = 0
    for i in range(k):
        distn += (xn[i] - yn[i])**2
    distn = sqrt(float(d)/k)*sqrt(distn)

    print "dist: ", dist, " distn: ", distn


