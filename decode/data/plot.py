#!/usr/bin/python

from sys import *
from pylab import *
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as mcolor

rc('font', family='serif')
rc('text', usetex='True')

times = []
NUM_DIM = 1

if len(argv) > 1:
    save_name = argv[1]
else:
    save_name = "none"


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

def plot_graph(fname):
    
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
    
    rrgp = np.array (rrgp)
    prob = np.array(prob)
    
    figure(1)    
    if( len(rrgp) > 0):
        for i in range(NUM_DIM-1):

            subplot(NUM_DIM-1,1,i+1, aspect='auto')
            plot(rrgp[:,0], rrgp[:,i+1], 'yo', ms=5.0, alpha = 0.1 )
            grid()

    """
    fig = figure(2)
    ax = fig.add_subplot(111, aspect=1.0)
    
    scatter( rrgp[:,0], rrgp[:,1], marker='o', c= prob[:], s=30, alpha=0.8)
    grid()
    axis([-2, 2, -2, 2])
    
    fig.savefig(fname)
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
        pf = []
        kf_covar = []

        tsys = []
        tobs = []
        tbp = []
        tkf = []
        tpf = []
        tkf_covar = []

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
                elif s[0] == "pf_path\n":
                    which = 4
                elif s[0] == "kf_covar\n":
                    which = 5

            if len(s) > 1:
                time = float(s[0])
                to_put = [ float(s[i+2]) for i in range(NUM_DIM) ]
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
                elif which == 4:
                    tpf.append(time)
                    pf.append(to_put)
                elif which == 5:
                    tkf_covar.append(time)
                    kf_covar.append(to_put)

    traj.close()

    sys = array(sys)
    obs = array(obs)
    bp = array(bp)
    kf = array(kf)
    pf = array(pf)
    kf_covar = array(kf_covar)
    
    tsys = array(tsys)
    tobs = array(tobs)
    tbp = array(tbp)
    tkf = array(tkf)
    tpf = array(tpf)
    tkf_covar = array(tkf_covar)
   
    times = tobs

    fig = figure(1)
    ax = fig.add_subplot(111, aspect='equal')
    setp(ax.get_xticklabels(), fontsize=20)
    setp(ax.get_yticklabels(), fontsize=20)
    if len(sys) != 0:
        plot( tsys, sys[:,0], 'r-', label='sys', lw=1.5)
    if len(bp) != 0:
        plot( tbp, bp[:,0], 'go-', label='decode', lw=1.5)
    #if len(obs) != 0:
        #plot( tobs, obs[:,0], 'bx', label='obs')
    grid()
    legend()
    axis('tight')
    xlabel(r'$t(s)$')
    ylabel(r'$x(t)$')
    title(r'$20,000$ samples')
    savefig('x1.pdf', bbox_inches='tight')

    if NUM_DIM == 2:
        fig = figure(2)
        ax = fig.add_subplot(111, aspect='equal')
        setp(ax.get_xticklabels(), fontsize=20)
        setp(ax.get_yticklabels(), fontsize=20)
        if len(sys) != 0:
            plot( tsys, sys[:,1], 'r-', label='sys', lw=1.5)
        if len(bp) != 0:
            plot( tbp, bp[:,1], 'go-', label='decode', lw=1.5)
        #if len(obs) != 0:
            #plot( tobs, obs[:,1], 'bx', label='obs')
        grid()
        legend()
        axis('tight')
        xlabel(r'$t$(s)')
        ylabel(r'$x(t)$')
        title(r'$20,000$ samples')

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
                    for i in range(NUM_DIM-1):
                        subplot(NUM_DIM-1,1,i+1, aspect='auto')
                        grid()
                        plot(to_plot_time[:], to_plot[:,i+1], 'm-', alpha=last_prob)
                   
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

    ep = open("err_vanderpol.txt", "r")
    tp = open("timing_inc.dat", "r")

    times=[]
    vert=[]
    bpe = []
    kfe=[]

    if ep:
        lines = ep.readlines()
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

        figure(4)
        grid()
        plot(vert[:], bpe[:], 'g-', lw=1.5)
        plot(vert[:], kfe[:], 'r-', lw=1.5)
    
    times = []
    vert = []
    
    if tp:
        lines = tp.readlines()
        for l in lines:
            s = l.split('\t')
            
            n = float(s[0])
            vert.append(n)
            times.append( float(s[1])/n/pow(log(n),4))

        vert = array(vert)
        times = array(times)
        
        figure(3)
        grid()
        plot(vert[:], times[:], 'b-', lw=1.5)

def plot_density(fname):
    
    denp = []
    denf = open("data/density.dat", 'r')
    if denf:
        lines = denf.readlines()
        for l in lines:
            s = l.split('\t')
            to_put = [ float(s[i]) for i in range(NUM_DIM+1) ]
            denp.append( to_put )
        
    denf.close()
    
    denp = np.array (denp)
    minp = min(denp[:,0])
    maxp = max(denp[:,0])
    size = [ (3+(x-minp)/(maxp-minp)*1000) for x in denp[:,0] ]
    
    """
    fig = figure(2)
    ax = fig.add_subplot(111, aspect=1.0)
    scatter( denp[:,1], denp[:,2], marker='o', c='b', s= 25, alpha=0.7)
    axis([0, 1, 0, 1])
    grid()

    fig.savefig(fname)
    """
    
    fig = figure(2)
    hexbin(denp[:,1], denp[:,2], bins='log', cmap=cm.get_cmap('Jet'), alpha=0.9, mincnt=1)
    axis([-0.5, 1, -0.5, 1])
    #colorbar()
    fig.savefig(fname)

if __name__ == "__main__":

    plot_trajs()
    #plot_graph(save_name)
    #plot_sim_trajs()
    #draw_obstacles()    
    
    #do_timing_plot()

    #plot_density(save_name)

    #legend()
    
    """
    if save_name != "none":
        fig.savefig(save_name)
    """

    show()

