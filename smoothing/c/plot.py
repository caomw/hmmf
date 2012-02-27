#!/usr/bin/python

import sys
from numpy import *
from pylab import *

rc('font', family='serif')
rc('text', usetex='True')

to_save = False
if len(sys.argv) >= 3:
    to_save = bool(sys.argv[1])

truth = loadtxt('truth.dat')
obs = loadtxt('observations.dat')
fest = loadtxt('filter.dat')
sest = loadtxt('smoothing.dat')
pfest = loadtxt('pfilter.dat')

d = 1
if len(sys.argv) >= 2:
    d = float(sys.argv[1])
dt = d/float(len(truth[:,0]))
times = linspace(dt, d+dt, len(truth[:,0]))

#print "dt: ", dt
#print times

if len(truth[0,:]) > 0:
    figure(1)
    plot(times, truth[:,0], 'r-', label='sys')
    #plot(times, obs[:,0], 'b-', label='obs')
    plot(times, fest[:,0], 'g-', label='filter')
    #plot(times, sest[:,0], 'c-', label='smoothing')
    plot(times, pfest[:,0], 'y-', label='pf')
    axis('tight')
    grid()
    xlabel(r'$t(s)$')
    legend(loc=4)
    title(r'$x_1(t)$')
    if to_save:
        savefig('x1.pdf', bbox_inches='tight')

if len(truth[0,:]) > 1:
    figure(2)
    plot(times, truth[:,1], 'r-', label='sys')
    #plot(times, obs[:,1], 'b-', label='obs')
    plot(times, fest[:,1], 'g-', label='filter')
    #plot(times, sest[:,1], 'c-', label='smoothing')
    plot(times, pfest[:,1], 'y-', label='pf')
    axis('tight')
    grid()
    xlabel(r'$t(s)$')
    legend(loc=1)
    title(r'$x_2(t)$')
    if to_save:
        savefig('x2.pdf', bbox_inches='tight')
    
"""
if len(truth[0,:]) > 2:
    figure(3)
    plot(times, truth[:,2], 'r-', label='sys')
    plot(times, fest[:,2], 'g-', label='filter')
    #plot(times, sest[:,2], 'c-', label='smoothing')
    plot(times, pfest[:,2], 'y-', label='pf')
    axis('tight')
    grid()
    xlabel(r'$t(s)$')
    legend(loc=4)
    title(r'$x_3(t)$')
    if to_save:
        savefig('x3.pdf', bbox_inches='tight')

if len(truth[0,:]) > 3:
    figure(4)
    plot(times, truth[:,3], 'r-', label='sys')
    plot(times, fest[:,3], 'g-', label='filter')
    #plot(times, sest[:,3], 'c-', label='smoothing')
    plot(times, pfest[:,3], 'y-', label='pf')
    axis('tight')
    grid()
    xlabel(r'$t(s)$')
    legend(loc=4)
    title(r'$x_4(t)$')
    if to_save:
        savefig('x4.pdf', bbox_inches='tight')
"""

figure(5)
plot(truth[:,0], truth[:,1], 'r-', label='sys')
#plot(fest[:,0], fest[:,1], 'g-', label='filter')
plot(pfest[:,0], pfest[:,1], 'y-', label='pf')
axis('tight')
grid()
xlabel(r'$x(t)$')
ylabel(r'$y(t)$')
legend(loc=1)


ferr = truth - fest
#serr = truth - sest
pferr = truth - pfest
n1 = [norm(ferr[i,:])*norm(ferr[i,:]) for i in range(len(ferr[:,0]))]
#n2 = [norm(serr[i,:])*norm(serr[i,:]) for i in range(len(serr[:,0]))]
n3 = [norm(pferr[i,:])*norm(pferr[i,:]) for i in range(len(pferr[:,0]))]

n2 = 0
print "ferr: ", mean(n1)*d, " serr: ", mean(n2)*d, " pferr: ", mean(n3)*d

show()
