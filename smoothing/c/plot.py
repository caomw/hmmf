#!/usr/bin/python

import sys
from numpy import *
from pylab import *

to_save = False
if len(sys.argv) >= 3:
    to_save = bool(sys.argv[1])

truth = loadtxt('truth.dat')
obs = loadtxt('observations.dat')
fest = loadtxt('filter.dat')
sest = loadtxt('smoothing.dat')

d = 1
if len(sys.argv) >= 2:
    d = float(sys.argv[1])
dt = d/float(len(truth[:,0]))
times = linspace(dt, d+dt, len(truth[:,0]))

#print "dt: ", dt
#print times

figure(1)
plot(times, truth[:,0], 'r-', label='sys')
#plot(times, obs[:,0], 'b-', label='obs')
plot(times, fest[:,0], 'g-', label='filter')
plot(times, sest[:,0], 'c-', label='smoothing')
axis('tight')
grid()
xlabel('t [s]')
legend(loc=4)
title('vanderpol_x')
if to_save:
    savefig('smoothing_vanderpol_x.pdf', bbox_inches='tight')

figure(2)
plot(times, truth[:,1], 'r-', label='sys')
#plot(times, obs[:,1], 'b-', label='obs')
plot(times, fest[:,1], 'g-', label='filter')
plot(times, sest[:,1], 'c-', label='smoothing')
axis('tight')
grid()
xlabel('t [s]')
legend(loc=4)
title('vanderpol_x_dot')
if to_save:
    savefig('smoothing_vanderpol_x_dot.pdf', bbox_inches='tight')

"""
figure(3)
plot(times, truth[:,2], 'r-', label='sys')
plot(times, fest[:,2], 'g-', label='filter')
plot(times, sest[:,2], 'c-', label='smoothing')
axis('tight')
grid()
xlabel('t [s]')
legend(loc=4)
title('vanderpol_mu')
if to_save:
    savefig('smoothing_vanderpol_mu.pdf', bbox_inches='tight')
"""

ferr = truth - fest
serr = truth - sest
n1 = [norm(ferr[i,:])*norm(ferr[i,:]) for i in range(len(ferr[:,0]))]
n2 = [norm(serr[i,:])*norm(serr[i,:]) for i in range(len(serr[:,0]))]

print "ferr: ", mean(n1)*d, " serr: ", mean(n2)*d

show()
