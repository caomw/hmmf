#!/usr/bin/python

from sys import *
from pylab import *

times = []
NUM_DIM = 2

if len(argv) > 1:
    save_name = argv[1]
else:
    save_name = "none"


fig = figure(1)

def draw_obstacles():
    fig = figure(2)
    axplot = fig.add_subplot(111)
    rect = Rectangle( (0.56, 0.4), 0.14, 0.2, fc='blue', alpha = 0.2)
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
    
    global times
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
    
    """ 
    figure(1)    
    if( len(rrgp) > 0):
        for i in range(len(times)):
            for j in range(NUM_DIM):
                
                tmp = [times[i] for x in range(len(rrgp))]
                tmp = array(tmp)
                subplot(NUM_DIM,1,j+1, aspect='auto')
                plot(tmp, rrgp[:,j], 'yo', ms=5.0, alpha = 0.1 )
                grid()
    """
    """
    figure(2)
    plot(rrgp[:,0], rrgp[:,1], 'yo', ms=5.0, alpha=0.1)
    grid()
    """

def plot_trajs():
    
    global times
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
   
    times = tobs

    figure(1)    
    for i in range(NUM_DIM):
        
        subplot(NUM_DIM,1,i+1, aspect='auto')
        grid()

        if len(sys) != 0:
            plot( tsys[:], sys[:,i], 'r-', label='sys', lw=1.5)
        #if len(obs) != 0:
            #plot( tobs[:], obs[:,i], 'bx', label='obs')
        
        if len(bp) != 0:
            plot( tbp[:], bp[:,i], 'g-', label='hmm', lw=1.5)
        if len(kf) != 0:
            plot( tkf[:], kf[:,i], 'c-', label='kf', lw=1.5)
    
    """
    figure(2)
    if len(sys) != 0:
        plot( sys[:,0], sys[:,1], 'r-', label='sys', lw=1.5)
    if len(bp) != 0:
        plot( bp[:,0], bp[:,1], 'g-', label='sys', lw=1.5)
    if len(kf) != 0:
        plot( kf[:,0], kf[:,1], 'c-', label='kf', lw=1.5)
    #if len(obs) != 0:
        #plot( obs[:,0], obs[:,1], 'bx', label='obs')
    grid()
    """

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
                    
                    figure(1)    
                    for i in range(NUM_DIM):
                        subplot(NUM_DIM,1,i+1, aspect='auto')
                        grid()
                        plot(to_plot_time[:], to_plot[:,i], 'm-', alpha=last_prob)
                   
                    """
                    figure(2)
                    plot(to_plot[:,0], to_plot[:,1], 'mo-', alpha=0.1)
                    grid()
                    """

                curr_traj = []
                curr_times = []
            elif(len(s) == 6):
                curr_times.append( float(s[0]) )
                to_put = [float(s[i+1]) for i in range(NUM_DIM)]
                curr_traj.append( to_put )

        mc.close()

def do_timing_plot():

    tp = open("timing_vanderpol.txt", "r")

    times=[]
    vert=[]
    bpe = []
    kfe=[]

    if tp:
        lines = tp.readlines()
        for l in lines:
            s = l.split('\t')

            vert.append( float(s[0]))
            times.append( float(s[1]))
            bpe.append( float(s[2]))
            kfe.append( float(s[3]))

        vert = array(vert)
        times = array(times)
        bpe = array(bpe)
        kfe = array(kfe)

        figure(3)
        grid()
        plot(vert[:], times[:], 'b-', lw=1.5)
        figure(4)
        grid()
        plot(vert[:], bpe[:], 'g-', lw=1.5)
        plot(vert[:], kfe[:], 'r-', lw=1.5)


if __name__ == "__main__":

    plot_trajs()
    plot_sim_trajs()
    #draw_obstacles()    
    
    #do_timing_plot()

    #plot_graph()
    
    #figure(2)
    #axis([-10, 10, -10, 10])

    #legend()
    
    if save_name != "none":
        fig.savefig(save_name)

    show()

