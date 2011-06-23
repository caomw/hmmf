#!/usr/bin/python

from sys import *
from pylab import *

if __name__ == "__main__":

    rrgp = open("rrgp.dat", 'r')
    traj = open("traj.dat", 'r')

    if len(argv) > 1:
        NUM_DIM = int(argv[1])
    else:
        NUM_DIM = 2

    draw_edges = 0
    if draw_edges:
        #t1, x1, t2, x2, prob, delt
        rrg = open("rrg.dat", 'r')
        if rrg:
            figure(1)
            lines = rrg.readlines()
            for l in lines:
                s = l.split('\t')
                tx = [float(s[0]),float(s[2])]
                ty = [float(s[1]), float(s[3])]
                if float(s[4]) > 2.0:
                    tmp_alpha = 1.0
                else:
                    tmp_alpha = float(s[4])/2.0
                plot(tx, ty, 'k-', alpha = tmp_alpha, lw=0.5 )

        rrg.close()

    rrgp = []
    rrgpf = open("rrgp.dat", 'r')
    if rrgpf:
        lines = rrgpf.readlines()
        for l in lines:
            s = l.split('\t')
            to_put = [ float(s[i]) for i in range(NUM_DIM) ]
            rrgp.append( to_put )

    rrgpf.close()

    if traj:
        lines = traj.readlines()
        which = 0

        sys = []
        obs = []
        bp = []
        kf = []

        for l in lines:
            s= l.split('\t')
            if len(s) == 1:
                if s[0] == "system\n":
                    which = 0
                elif s[0] == "observation\n":
                    which = 1
                elif s[0] == "best_path\n":
                    which = 2
                elif s[0] == "kf_path\n":
                    which = 3

            if len(s) > 1:
                to_put = [ float(s[i]) for i in range(NUM_DIM) ]
                if which == 0:
                    sys.append(to_put)
                elif which == 1:
                    obs.append(to_put)
                elif which == 2:
                    bp.append(to_put)
                elif which == 3:
                    kf.append(to_put)

    traj.close()



    """
    axplot = fig.add_subplot(111)
    avgx, avgy = 0.0, 0.0
    for i in range(len(aa)):
        avgx = avgx + ax[i]*aa[i]
        avgy = avgy + ay[i]*aa[i]
        circle = Circle( (ax[i], ay[i]), 0.01, fc='blue', alpha = 10*aa[i])
        axplot.add_patch(circle)

    circle = Circle( (avgx, avgy), 0.01, fc='green', alpha = 0.5)
    axplot.add_patch(circle)
    """
    
    rrgp = array (rrgp)
    sys = array(sys)
    obs = array(obs)
    bp = array(bp)
    kf = array(kf)
    
    PLOT_RRG = 1
    SAVE = 0

    for i in range(NUM_DIM-1):
        
        if i == 0:
            subplot(311, aspect='auto')
            if PLOT_RRG:
                plot(rrgp[:,0], rrgp[:,1], 'yo', ms=3.0)
            ylabel('x')
        elif i == 1:
            subplot(312, aspect='auto')
            if PLOT_RRG:
                plot(rrgp[:,0], rrgp[:,2], 'yo', ms=3.0)
            ylabel('y')
        elif i == 2:
            subplot(313, aspect='auto')
            if PLOT_RRG:
                plot(rrgp[:,0], rrgp[:,3], 'yo', ms=3.0)
            ylabel('th')

        plot( sys[:,0], sys[:,i+1], 'r-', label='sys')
        plot( obs[:,0], obs[:,i+1], 'b-', label='obs')
        
        if len(bp) != 0:
            plot( bp[:,0], bp[:,i+1], 'g-', label='hmm')
        if len(kf) != 0:
            plot( kf[:,0], kf[:,i+1], 'c-', label='kf')
        

        #legend() 
        grid()

    if SAVE:
        savefig("run.png")
    show()
    
    """
    close_flag = raw_input("Press q to close")
    if(close_flag == 'q'):
        close('all')
    exit(0)
    """
