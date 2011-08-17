#!/usr/bin/python

from sys import *
from pylab import *

SAVE = 1

if len(argv) > 1:
    NUM_DIM = int(argv[1])
else:
    NUM_DIM = 2


fig = figure(1)

def draw_obstacles():
    axplot = fig.add_subplot(111)
    rect = Rectangle( (.127, 0), .26-.127, .217, fc='blue', alpha = 0.2)
    axplot.add_patch(rect)
    rect = Rectangle( (.1, .32), .2 - .1, .5 - .32, fc='blue', alpha = 0.2)
    axplot.add_patch(rect)

def draw_edges():
    draw_edges = 0
    if draw_edges:
        #t1, x1, t2, x2, prob, delt
        rrg = open("rrg.dat", 'r')
        if rrg:
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

def plot_graph():

    rrgp = []
    rrgpf = open("rrgp.dat", 'r')
    if rrgpf:
        lines = rrgpf.readlines()
        for l in lines:
            s = l.split('\t')
            to_put = [ float(s[i]) for i in range(NUM_DIM) ]
            rrgp.append( to_put )
        
        rrgp_is_exists = 1;
    rrgpf.close()
    
    rrgp = array (rrgp)

    for i in range(NUM_DIM-1):
        
        subplot(NUM_DIM-1,1,i+1, aspect='auto')
        plot(rrgp[:,0], rrgp[:,i+1], 'yo', ms=3.0, alpha = 0.6)
        grid()

def plot_trajs():
    
    traj = open("traj.dat", 'r')

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

    sys = array(sys)
    obs = array(obs)
    bp = array(bp)
    kf = array(kf)
    
    for i in range(NUM_DIM-1):
        
        subplot(NUM_DIM-1,1,i+1, aspect='auto')
        grid()

        if len(sys) != 0:
            plot( sys[:,0], sys[:,i+1], 'r-', label='sys', lw=2.0)
        if len(obs) != 0:
            plot( obs[:,0], obs[:,i+1], 'b-', label='obs')
        
        if len(bp) != 0:
            plot( bp[:,0], bp[:,i+1], 'go-', label='hmm')
        if len(kf) != 0:
            plot( kf[:,0], kf[:,i+1], 'c-', label='kf')
        

def plot_sim_trajs():

    mc = open("monte_carlo.dat", 'r')
    
    probs = []
    curr_traj = []

    if mc:
        lines = mc.readlines()
        last_prob = 1
        curr_prob = 1
        for l in lines:
            s = l.split('\t')

            if(len(s) ==3):
                last_prob = curr_prob
                curr_prob = 5*float(s[1])
                to_plot = array(curr_traj)
                
                print curr_prob
                if( len(to_plot) > 0):
                    
                    for i in range(NUM_DIM-1):
                        subplot(NUM_DIM-1,1,i+1, aspect='auto')
                        plot(to_plot[:,0], to_plot[:,i+1], 'm-', alpha=last_prob)
                        grid()

                curr_traj = []
            elif(len(s) == 5):
                to_put = [float(s[i]) for i in range(NUM_DIM)]
                curr_traj.append( to_put )

        mc.close()

if __name__ == "__main__":

    #plot_graph()
    plot_trajs()
    plot_sim_trajs()
    

    xlabel('t [s] ')
    ylabel('x')
    #legend()
    
    if SAVE:
        fig.savefig("run.pdf")

    show()

