#!/usr/bin/python

from numpy import *
from pylab import *
from numpy.random import *

n = 1000
e = 0.8
d = n
beta = 2
k = int(log(n)/(e*e/2 - e*e*e/3)*(4 + 2*beta))
print "k: ", k

sq3 = sqrt(3/float(k))
R = zeros( (d, k), float)
for i in range(d):
    for j in range(k):
        t = randint(1, 7)       # 7 exclusive
        if (t == 1):
            R[i, j] = sq3
        elif (t == 6):
            R[i, j] = -sq3

print R.shape

old_d = []
new_d = []
for c in range(n):
    x  = [0 for i in range(d)]
    dist = 0 
    for i in range(d):
        x[i] = random()
        dist = dist + x[i]*x[i]
    dist = sqrt(dist)
    old_d.append(dist)

    xn = dot(x, R)
    distn = 0
    for i in range(k):
        distn = distn + xn[i]*xn[i]
    distn = sqrt(distn)
    new_d.append(distn)

    #print "dist: ", dist, " distn: ", distn

old_d = array(old_d)
new_d = array(new_d)
err = (old_d - new_d)/old_d
mean = average(err)
stddev = std(err)
print "mean_err: ", mean, " stddev: ", stddev

#plot(err, 'b-')
#show()
