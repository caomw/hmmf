#!/usr/bin/python

import sys
from numpy import *
from pylab import *

rc('font', family='serif')
rc('text', usetex='True')

to_save = False

truth = loadtxt('data/truth.dat')
obs = loadtxt('data/observations.dat')
try:
    fest = loadtxt('data/filter.dat')
except IOError:
    print 'no filter data'
try:
    sest = loadtxt('data/smoothing.dat')
except IOError:
    print 'no smoothing data'
try:
    pfest = loadtxt('data/pfilter.dat')
except IOError:
    print 'no pf data'
try:
    kfest = loadtxt('data/kfilter.dat')
except IOError:
    print 'no kf data'


d = 1
if len(sys.argv) > 1:
    if int(sys.argv[1]) == 1:
        to_save = True
    else:
        to_save = False
    d = float(sys.argv[2])
if truth.ndim == 1:
    dt = d/float(len(truth))
    times = linspace(dt, d+dt, len(truth))
else:
    dt = d/float(len(truth[:,0]))
    times = linspace(dt, d+dt, len(truth[:,0]))
#print "dt: ", dt
#print times

if truth.ndim == 1:
    fig = figure(1)
    ax = fig.add_subplot(111, aspect='equal')
    setp(ax.get_xticklabels(), fontsize=20)
    setp(ax.get_yticklabels(), fontsize=20)
    plot(times, truth, 'r-', label='sys')
    plot(times, obs, 'b-', label='obs')
    plot(times, fest, 'g-', label='filter')
    #plot(times, kfest, 'c-', label='kf')
    plot(times, sest, 'c-', label='smoothing')
    plot(times, pfest, 'y-', label='pf')
    axis('tight')
    grid()
    xlabel(r'$t(s)$')
    legend(loc=1)
    title(r'$x(t)$')
    if to_save:
        savefig('x1.pdf', bbox_inches='tight')

else:
    fig = figure(1)
    ax = fig.add_subplot(111, aspect='equal')
    setp(ax.get_xticklabels(), fontsize=20)
    setp(ax.get_yticklabels(), fontsize=20)
    plot(times, truth[:,0], 'r-', label='sys')
    #plot(times, obs[:,0], 'b-', label='obs')
    plot(times, fest[:,0], 'g-', label='filter')
    #plot(times, kfest[:,0], 'c-', label='kf')
    plot(times, sest[:,0], 'c-', label='smoothing')
    plot(times, pfest[:,0], 'y-', label='pf')
    axis('tight')
    grid()
    xlabel(r'$t(s)$')
    legend(loc=4)
    title(r'$x_1(t)$')
    if to_save:
        savefig('x1.pdf', bbox_inches='tight')

    if len(truth[0,:]) > 1:
        fig = figure(2)
        ax = fig.add_subplot(111, aspect='equal')
        setp(ax.get_xticklabels(), fontsize=20)
        setp(ax.get_yticklabels(), fontsize=20)
        plot(times, truth[:,1], 'r-', label='sys')
        #plot(times, obs[:,1], 'b-', label='obs')
        plot(times, fest[:,1], 'g-', label='filter')
        #plot(times, kfest[:,1], 'c-', label='kf')
        plot(times, sest[:,1], 'c-', label='smoothing')
        plot(times, pfest[:,1], 'y-', label='pf')
        axis('tight')
        grid()
        xlabel(r'$t(s)$')
        legend(loc=3)
        title(r'$x_2(t)$')
        if to_save:
            savefig('x2.pdf', bbox_inches='tight')

    """
    if len(truth[0,:]) > 2:
        fig = figure(3)
        ax = fig.add_subplot(111, aspect='equal')
        setp(ax.get_xticklabels(), fontsize=20)
        setp(ax.get_yticklabels(), fontsize=20)
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
        fig = figure(4)
        ax = fig.add_subplot(111, aspect='equal')
        setp(ax.get_xticklabels(), fontsize=20)
        setp(ax.get_yticklabels(), fontsize=20)
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
    """
    fig = figure(5)
    ax = fig.add_subplot(111, aspect='equal')
    setp(ax.get_xticklabels(), fontsize=20)
    setp(ax.get_yticklabels(), fontsize=20)
    plot(truth[:,0], truth[:,1], 'r-', lw=1.5, label='sys')
    plot(fest[:,0], fest[:,1], 'g-', lw=1.5, label='hmmf')
    plot(kfest[:,0], kfest[:,1], 'c-', lw=1.5, label='kf')
    plot(pfest[:,0], pfest[:,1], 'y-', lw=1.5, label='pf')
    axis('tight')
    grid()
    xlabel(r'$x(t)$')
    ylabel(r'$y(t)$')
    legend(loc=4)
    if to_save:
        savefig('ship.pdf', bbox_inches='tight')
    """

    ferr = truth - fest
    serr = truth - sest
    pferr = truth - pfest
    n1 = [norm(ferr[i,:])*norm(ferr[i,:]) for i in range(len(ferr[:,0]))]
    n2 = [norm(serr[i,:])*norm(serr[i,:]) for i in range(len(serr[:,0]))]
    n3 = [norm(pferr[i,:])*norm(pferr[i,:]) for i in range(len(pferr[:,0]))]

    #n2 = 0
    print "ferr: ", mean(n1)*d, " serr: ", mean(n2)*d, " pferr: ", mean(n3)*d

show()
