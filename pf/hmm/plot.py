#!/usr/bin/python

from pylab import *

if __name__ == "__main__":

    rrgp = open("rrgp.dat", 'r')
    traj = open("traj.dat", 'r')

    vx = []
    vy = []

    """
    rrg = open("rrg.dat", 'r')
    if rrg:
        lines = rrg.readlines(100000)
        for l in lines:
            s = l.split('\t')
            tx = [float(s[0]),float(s[2])]
            ty = [float(s[1]), float(s[3])]
            plot(tx, ty, 'yo', lw=0.5)

    rrg.close()
    """

    rrgpx = []
    rrgpy = []
    rrgp = open("rrgp.dat", 'r')
    if rrgp:
        lines = rrgp.readlines(100000)
        for l in lines:
            s = l.split('\t')
            rrgpx.append(float(s[0]))
            rrgpy.append(float(s[1]))

    rrgp.close()

    if traj:
        lines = traj.readlines(1000000)
        which = 0
        sysx = []
        sysy = []
        obsx = []
        obsy = []
        bpx = []
        bpy = []
        ax = []
        ay = []
        aa = []

    for l in lines:
        s= l.split('\t')
        if len(s) == 1:
            if s[0] == "system\n":
                which = 0
            elif s[0] == "observation\n":
                which = 1
            elif s[0] == "best_path\n":
                which = 2
            elif s[0] == "alpha\n":
                which = 3

        if len(s) > 1:
            if which == 0:
                sysx.append(float(s[0]))
                sysy.append(float(s[1]))
            elif which == 1:
                obsx.append(float(s[0]))
                obsy.append(float(s[1]))
            elif which == 2:
                bpx.append(float(s[0]))
                bpy.append(float(s[1]))
            elif which == 3:
                ax.append( float(s[0]))
                ay.append( float(s[1]))
                aa.append( float(s[2]))

    traj.close()

    fig = figure()
    axplot = fig.add_subplot(111)

    """
    avgx, avgy = 0.0, 0.0
    for i in range(len(aa)):
        avgx = avgx + ax[i]*aa[i]
        avgy = avgy + ay[i]*aa[i]
        circle = Circle( (ax[i], ay[i]), 0.01, fc='blue', alpha = 10*aa[i])
        axplot.add_patch(circle)
    
    circle = Circle( (avgx, avgy), 0.01, fc='green', alpha = 0.5)
    axplot.add_patch(circle)
    """
    
    plot(rrgpx, rrgpy, 'yo', lw=0.5)
    plot(sysx, sysy, 'ro-')
    plot(obsx, obsy, 'bo-')
    plot(bpx, bpy, 'go-')

    grid()
    show()

