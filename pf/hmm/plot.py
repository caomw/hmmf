#!/usr/bin/python

from pylab import *

if __name__ == "__main__":

    fig = figure()
    rrgp = open("rrgp.dat", 'r')
    traj = open("traj.dat", 'r')
    mini = open("miniout.dat", 'r')

    minix = []
    miniy = []
    if mini:
        lines = mini.readlines()
        for l in lines:
            s = l.split('\t')
            minix.append( float(s[0]))
            miniy.append( float(s[1]))

        # print len(minix)
        mini.close()
    """   
    #t1, x1, t2, x2, prob, delt
    rrg = open("rrg.dat", 'r')
    if rrg:
        lines = rrg.readlines()
        for l in lines:
            s = l.split('\t')
            tx = [float(s[0]),float(s[2])]
            ty = [float(s[1]), float(s[3])]
            plot(tx, ty, 'k-', lw=0.5, alpha=float(s[4]))

    rrg.close()

    rrgpx = []
    rrgpy = []
    rrgp = open("rrgp.dat", 'r')
    if rrgp:
        lines = rrgp.readlines()
        for l in lines:
            s = l.split('\t')
            rrgpx.append(float(s[0]))
            rrgpy.append(float(s[1]))

    rrgp.close()
    """

    if traj:
        lines = traj.readlines()
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
        simx = []
        simy = []
        kfx = []
        kfy = []

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
            elif s[0] == "sim\n":
                which = 4
            elif s[0] == "kf_path\n":
                which = 5

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
            elif which == 4:
                simx.append(float(s[0]))
                simy.append(float(s[1]))
            elif which == 5:
                kfx.append(float(s[0]))
                kfy.append(float(s[1]))

    traj.close()

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
    
    #plot(rrgpx, rrgpy, 'yo', ms=3.0)
    plot(sysx, sysy, 'r-', label='sys')
    plot(obsx, obsy, 'bo-', label='sys')
    plot(bpx, bpy, 'g-', label='sys')
    plot(kfx, kfy, 'c-', label='sys')
    
    xlabel('t [s]')
    ylabel('x (t)')
    title('Filter output')
    legend() 
    grid()
    fig.savefig("run.png")
    show()

