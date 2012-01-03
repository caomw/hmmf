#!/usr/bin/python

from numpy import *
from pylab import *

dT = 0.0001

x = linspace(0,1,1/dT)
y = 0*x
for i in range(len(x)-1):
    y[i+1] = y[i] + normal(0, dT)

fig = figure(1)
plot(x,y)
grid()
show()
fig.savefig('brownian.pdf')
