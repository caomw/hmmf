#!/usr/bin/python

from sys import *
from pylab import *

NUM_DIM = 2

if len(argv) > 1:
    save_name = argv[1]
else:
    save_name = "none"


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
    prob = []
    rrgpf = open("rrgp.dat", 'r')
    if rrgpf:
        lines = rrgpf.readlines()
        for l in lines:
            s = l.split('\t')
            to_put = [ float(s[i]) for i in range(NUM_DIM) ]
            rrgp.append( to_put )
            prob.append( float(s[NUM_DIM]) )
        
    rrgpf.close()
    
    rrgp = array (rrgp)
    prob = array(prob)
    
    if( len(rrgp) > 0):
        for i in range(NUM_DIM-1):
        
            subplot(NUM_DIM-1,1,i+1, aspect='auto')
            plot(rrgp[:,0], rrgp[:,i+1], 'yo', ms=5.0, alpha = 0.1 )
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
        tsys = []
        tobs = []
        tbp = []
        tkf = []

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
                time = float(s[0])
                to_put = [ float(s[i+1]) for i in range(NUM_DIM) ]
                if which == 0:
                    tsys.append(time)
                    sys.append(to_put)
                elif which == 1:
                    tobs.append(time)
                    obs.append(to_put)
                elif which == 2:
                    tbp.append(time)
                    bp.append(to_put)
                elif which == 3:
                    tkf.append(time)
                    kf.append(to_put)

    traj.close()

    sys = array(sys)
    obs = array(obs)
    bp = array(bp)
    kf = array(kf)
    
    tsys = array(tsys)
    tobs = array(tobs)
    tbp = array(tbp)
    tkf = array(tkf)
    

    for i in range(NUM_DIM):
        
        subplot(NUM_DIM,1,i+1, aspect='auto')
        grid()

        if len(sys) != 0:
            plot( tsys[:], sys[:,i], 'r-', label='sys', lw=2.0)
        #if len(obs) != 0:
            #plot( tobs[:], obs[:,i], 'b-', label='obs')
        
        if len(bp) != 0:
            plot( tbp[:], bp[:,i], 'g-', label='hmm', lw=2.0)
        if len(kf) != 0:
            plot( tkf[:], kf[:,i], 'c-', label='kf', lw=2.0)
        

def plot_sim_trajs():

    mc = open("monte_carlo.dat", 'r')
    
    probs = []
    curr_traj = []
    curr_times = []

    if mc:
        lines = mc.readlines()
        last_prob = 1
        curr_prob = 1
        for l in lines:
            s = l.split('\t')

            if(len(s) ==3):
                last_prob = curr_prob
                curr_prob = 50*float(s[1])
                to_plot = array(curr_traj)
                to_plot_time = array(curr_times)

                if( len(to_plot) > 0):
                    
                    for i in range(NUM_DIM):
                        subplot(NUM_DIM,1,i+1, aspect='auto')
                        grid()
                        plot(to_plot_time[:], to_plot[:,i], 'm-', alpha=0.01)

                curr_traj = []
                curr_times = []
            elif(len(s) == 6):
                curr_times.append( float(s[0]) )
                to_put = [float(s[i+1]) for i in range(NUM_DIM)]
                curr_traj.append( to_put )

        mc.close()

if __name__ == "__main__":

    plot_trajs()
    plot_sim_trajs()
    
    #plot_graph()

    #legend()
    
    if save_name != "none":
        fig.savefig(save_name)

    show()

