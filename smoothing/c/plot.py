#!/usr/bin/python

from numpy import *
from pylab import *

truth = loadtxt('truth.dat')
obs = loadtxt('observations.dat')
est = loadtxt('filter.dat')

figure(1)
plot(truth[:,0], 'r-')
plot(obs[:,0], 'b-')
plot(est[:,0], 'g-')
grid()

show()
