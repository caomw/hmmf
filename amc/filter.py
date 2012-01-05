#!/usr/bin/python

import sys
from time import *
from numpy import *
from random import *
from math import *
from scipy.spatial import kdtree
from pylab import *

try:
    import cPickle as pickle
except:
    import pickle

NUM_DIM=1
zmin = [0.0]
zmax = [1.0]
init_state = [0.5]
init_var = array(diag([0.001]))
consant_holding_time = 0

def wstd(arrin, weights_in):
    arr = array(arrin, ndmin=1, copy=False)
    weights = array(weights_in, ndmin=1, dtype='f8', copy=False)
  
    wtot = weights.sum()
    wmean = ( weights*arr ).sum()/wtot

    wvar = ( weights*(arr-wmean)**2 ).sum()/wtot
    wsdev = sqrt(wvar)
    return wmean,wsdev

def fdt(z, u, dt):
    return array((-z +u)*dt)
def get_process_var(dt):
    return diag(array([0.001*dt]))
def get_observation_var():
    return [0.001]
def get_observation(z):
    var = get_observation_var()
    return [z[i] + normal(0, var[i]) for i in range(NUM_DIM)]
def get_holding_time(z, u, bowlr):
    bowlr = bowlr*(zmax[0]-zmin[0])
    dt = 1
    if NUM_DIM >=  2:
        return bowlr*bowlr/(bowlr*linalg.norm(fdt(z,u,dt)) + trace(get_process_var(dt)))
    else:
        return bowlr*bowlr/(bowlr*linalg.norm(fdt(1,0,dt)) + get_process_var(dt)[0,0])
        #return bowlr*bowlr/(bowlr*linalg.norm(fdt(z,u,dt)) + get_process_var(dt)[0,0])

def normal_val(x, mu, var):
    dim = len(mu)
    delx = x -mu;
    if len(delx) >= 2:
        det = sqrt(linalg.det(var))
        toret = 1/pow(2*pi,dim/2.0)/det*exp(-0.5*dot(delx,dot(linalg.inv(var),delx)))
        return toret
    else:
        det = sqrt(var[0])
        toret = 1/pow(2*pi,dim/2.0)/det*exp(-0.5*dot(delx,delx)/det/det)
        return toret[0]


class Edge:
    transition_probability = 1
    transition_time = 0
    transition_input = 0;

    def __init__(self, prob, time):
        transition_probability = prob
        transition_time = time
        transition_input = []

class Node:
    

    def __init__(self, is_init=False):
        self.x = []
        self.density = 0
        self.edge_probs = []
        self.edge_times = []
        self.edge_controls = []
        self.pself = 1

        if is_init == False:
            self.x = array([zmin[i] + random()*(zmax[i]-zmin[i]) for i in range(NUM_DIM)])
        else:
            self.x = array( [init_state[i] + normal(0, sqrt(init_var[i,i])) for i in range(NUM_DIM)] )
        # self.density = normal_val(self.x, array(init_state), init_var)
        # print self.x
    
class Graph:

    def key(self, z):
        return [(z[i] - zmin[i])/(zmax[i]-zmin[i]) for i in range(NUM_DIM)]

    def __init__(self, vert):
        self.num_vert = vert
        self.bowlr = 2.1*pow(log(self.num_vert)/float(self.num_vert),1/float(NUM_DIM))
        self.delta = 0.1

        self.nodes = []
        self.point = []

        self.tree = []
        self.mydict = []

        self.truth = []
        self.observations = []
        self.estimates = []

        for i in range(self.num_vert):
            self.nodes.append(Node(False))
        
        # initial variance
        for i in range(self.num_vert):
            n1 = self.nodes[i]
            n1.density = normal_val(n1.x, array(init_state), init_var)
        
        self.normalize_density()

        self.points = [self.key(mynode.x) for mynode in self.nodes]
        self.tree = kdtree.KDTree(self.points)
        
    def draw_edges(self, n1):
        # clear old edges
        n1.edge_probs = []
        n1.edge_times = []

        # n1 is the key of mydict
        probs = []
        holding_time = get_holding_time(n1.x, 0, self.bowlr)
        pself = holding_time/(holding_time + self.delta)
        n1.pself = pself

        for n2 in self.mydict[n1]:
                mu = n1.x + fdt(n1.x, 0, holding_time)
                curr_prob = normal_val(n2.x, mu, 
                                    get_process_var(holding_time))
                probs.append(curr_prob)
 
        tot_prob = sum(probs)
        probs = probs/tot_prob*(1 - pself)
        probs = list(probs)
        probs.append(pself)

        all_neighbors = self.mydict[n1]
        all_neighbors.append(n1)
        self.mydict[n1] = all_neighbors
        n1.edge_probs.append(probs)

        transition_times = [holding_time*self.delta/(holding_time + self.delta) for i in range(len(probs))]
        n1.edge_times.append(transition_times)

    def connect(self):
        for n1 in self.nodes:
            neighbors = []
            neighbors_index = self.tree.query_ball_point(self.key(n1.x), self.bowlr)
            for n2_index in neighbors_index:
                    if n1 != self.nodes[n2_index]:
                        neighbors.append( self.nodes[n2_index])
            self.mydict.append((n1, neighbors))
        
        self.mydict = dict(self.mydict)
        
        c = 0
        for n1 in self.mydict.keys():
            self.draw_edges(n1)
            c = c+1
            """
            if c%100 == 0:
                print c
            """
        #print "finished connecting"
        count = 0
        count = count+1
        
    def print_graph(self):

        for n1 in self.mydict.keys():
            print n1.x
            count = 0
            for n2 in self.mydict[n1]:
                print "\t", n2.x,'\t', n1.edge_times[0][self.mydict[n1].index(n2)], " ", n1.edge_probs[0][self.mydict[n1].index(n2)]
                count = count + 1
            raw_input()

    def find_closest_index(self, mylist, myvar):
        tmp = [ abs(mylist[i] - myvar) for i in range(len(mylist))]
        return tmp.index(min(tmp))

    def find_nearest_node(self, state):
        key = self.key(state)
        dist, index = self.tree.query(key, 1)
        return self.nodes[index]

    def simulate_trajectories(self):
        # utraj is a array of (t, u)
        trajs = []
        traj_probs = []
        traj_times = []
        max_time = 0.5
        for traj_index in range(1000):
            
            #if traj_index%100 == 0:
            #    print traj_index

            traj_curr = []
            curr_time = 0
            traj_time = []
            start_state = normal(init_state, sqrt(init_var[0]))
            node_curr = self.find_nearest_node(start_state)
            traj_prob = normal_val(node_curr.x, array(init_state), init_var)
                
            traj_curr.append(list(node_curr.x))
            traj_time.append(curr_time)
            # print "changed node to: ", node_curr.x, " time to: ", curr_time

            while curr_time < max_time:
                
                tmp_probs = node_curr.edge_probs[0]

                cum_probs = cumsum(tmp_probs)
                coin_toss = random()
                next_index = 0
                for i in range(len(cum_probs)):
                    if coin_toss >= cum_probs[i]:
                        next_index = next_index+1
                    else:
                        break
                
                #print len(self.mydict[node_curr]), next_index
                traj_prob = traj_prob * node_curr.edge_probs[0][next_index]
                curr_time = curr_time + node_curr.edge_times[0][next_index]
                node_curr = self.mydict[node_curr][next_index]
                # print "changed node to: ", node_curr.x, " time to: ", curr_time

                traj_curr.append(list(node_curr.x))
                traj_time.append(curr_time)

            to_put = [item for sublist in traj_curr for item in sublist]
            while len(to_put) < 300:
                to_put.append(to_put[-1])
                traj_time.append(traj_time[-1])

            trajs.append(to_put)
            traj_probs.append(traj_prob)
            traj_times.append(traj_time)
            
            #plot(traj_time, to_put,'b--')
        
        trajs = array(trajs)
        traj_probs = array(traj_probs)
        #traj_avg = average(trajs, axis=0, weights=traj_probs)
        traj_avg = array([wstd(trajs[:,i], traj_probs)[0] for i in range(len(trajs[0,:]))])
        traj_std = array([wstd(trajs[:,i], traj_probs)[1] for i in range(len(trajs[0,:]))])
        
        return traj_avg[-1], traj_std[-1]

        """
        fig = figure(1)
        grid()
        plot(traj_times[0], traj_avg, 'b-', label='mean')
        plot(traj_times[0], traj_avg-traj_std, 'b--', label='+/- std')
        plot(traj_times[0], traj_avg+traj_std, 'b--')

        cont_time = linspace(0,max_time,1000)
        cont_mean = init_state[0]*exp(-linspace(0,max_time,1000))
        cont_std = sqrt(array([0.001/2 + 0.001/2*exp(-2*curr_time) for curr_time in cont_time]))
        plot(cont_time, cont_mean, 'r-', label='cont. mean')
        plot(cont_time, cont_mean+cont_std, 'r--', label='cont. +/- std')
        plot(cont_time, cont_mean-cont_std, 'r--')
        axis('tight')
        xlabel('t [s]')
        ylabel( 'x(t)')
        legend()
        show()
        #fig.savefig('singleint_mc_convergence.pdf', bbox_inches='tight') 
        """

    def update_conditional_density(self, curr_observation):
        for n1 in self.nodes:
            transition_prob = 0
            for n2 in self.mydict[n1]:
                if n2 != n1:
                    transition_prob = transition_prob + n2.density*n2.edge_probs[0][self.mydict[n2].index(n1)]
        
            transition_prob = transition_prob + n1.density
            n1.density = transition_prob*n1.pself*normal_val(n1.x, array(curr_observation), array(diag(get_observation_var())))
        
    def normalize_density(self):
        tot_prob = 0
        for i in range(self.num_vert):
            tot_prob = tot_prob + self.nodes[i].density
        for i in range(self.num_vert):
            self.nodes[i].density = self.nodes[i].density/tot_prob
   
    def get_mean_std(self):
        mean_std = []
        weights = array([n1.density for n1 in self.nodes])
        for i in range(NUM_DIM):
            xdim = array([n1.x[i] for n1 in self.nodes])
            mean_std.append(wstd(xdim, weights))
        return mean_std

    def run_filter(self):
        max_time = 1
        curr_time = 0
        integration_delta = 1e-2

        self.truth.append(init_state)
        mean_std = self.get_mean_std()
        self.estimates.append(mean_std[0][0])
        self.observations.append( get_observation(init_state))
        while curr_time < max_time:
            curr_state = self.truth[-1]
            next_state = curr_state

            runner_time = 0
            while runner_time < self.delta:
                next_state = [next_state[i] + (-next_state[i]*integration_delta + normal(0, sqrt(get_process_var(integration_delta))[0]) ) for i in range(NUM_DIM)]
                runner_time = runner_time + integration_delta

            self.truth.append(next_state)
            self.observations.append( get_observation(next_state))

            self.update_conditional_density(self.observations[-1])
            self.normalize_density()

            mean_std = self.get_mean_std()
            self.estimates.append(mean_std[0][0])

            curr_time = curr_time + self.delta


        plot(self.truth, 'ro-')
        plot(self.observations, 'bo-')
        plot(self.estimates, 'go-')
        for n1 in self.nodes:
            plot(0.5, n1.x[0],'yo')
        grid()
        show()

def do_convergence_error_plot():
    num_nodes = linspace(100, 10000, 10)
    state_mean = 0*num_nodes
    state_std = 0*num_nodes
    for i in range(len(num_nodes)):
        graph = Graph(int(num_nodes[i]))
        graph.connect()
        state_mean[i], state_std[i] = graph.simulate_trajectories()
        print i

    plot(num_nodes, state_mean, 'r-', label='mean')
    plot(num_nodes, state_mean+state_std, 'r--', label='+/- std')
    plot(num_nodes, state_mean-state_std, 'r--')
    grid()
    show()

if __name__ == "__main__":
        
    # patch kdtree
    kdtree.node = kdtree.KDTree.node
    kdtree.leafnode = kdtree.KDTree.leafnode
    kdtree.innernode = kdtree.KDTree.innernode

    if len(sys.argv) >=3:
        if sys.argv[1] == 'w':
            tic = clock()
            graph = Graph(int(sys.argv[2]))
            graph.connect()
            
            # graph.print_graph()
            pickle.dump(graph, open('graph.pkl','wb'))
            print clock() - tic, '[s]'

            # graph.print_graph()

    else:
        # graph = pickle.load(open('graph.pkl','rb'))

        # graph.simulate_trajectories()
        # graph.run_filter()
        do_convergence_error_plot()

