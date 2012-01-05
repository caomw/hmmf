#!/usr/bin/python

"""
Implements smoothing using approximating Markov chains
1.  sample n states, create basic MC
2.  copy this to m different times partitioning [0, t_max]
3.  run alphas from left, betas from right to get smoothened estimate
4.  add new state + refine_state to all m times, approximate alpha + beta for them all
5.  add new time, copy S_n for this time, update all alphas forward, betas backward from this time
"""

import sys
from time import *
from numpy import *
from random import *
from math import *
from scipy.spatial import kdtree
from pylab import *
import cPickle as pickle

# without the time
NUM_DIM=1
zmin = [-2.0]
zmax = [2.0]
init_state = [0.5]
init_var = array(diag([0.01]))

# generic
def find_closest_index(self, mylist, myvar):
    tmp = [ abs(mylist[i] - myvar) for i in range(len(mylist))]
    return tmp.index(min(tmp))
def wstd(arrin, weights_in):
    arr = array(arrin, ndmin=1, copy=False)
    weights = array(weights_in, ndmin=1, dtype='f8', copy=False)
    wtot = weights.sum()
    wmean = ( weights*arr ).sum()/wtot
    wvar = ( weights*(arr-wmean)**2 ).sum()/wtot
    wsdev = sqrt(wvar)
    return wmean,wsdev
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

# system
def fdt(z, u, dt):
    return array((-z +u)*dt)
def get_process_var(dt):
    return diag(array([0.01*dt]))
def get_observation_var():
    return [0.01]
def get_observation(z, is_clean):
    if not is_clean:
        var = get_observation_var()
        return [z[i] + normal(0, var[i]) for i in range(NUM_DIM)]
    else:
        return [z[i] for i in range(NUM_DIM)]
def get_holding_time(z, u, bowlr):
    bowlr = bowlr*(zmax[0]-zmin[0])
    dt = 1
    if NUM_DIM >=  2:
        return bowlr*bowlr/(bowlr*linalg.norm(fdt(z,u,dt)) + trace(get_process_var(dt)))
    else:
        # return bowlr*bowlr/(bowlr*linalg.norm(fdt(0,1,dt)) + get_process_var(dt)[0,0])
        return bowlr*bowlr/(bowlr*linalg.norm(fdt(z,u,dt)) + get_process_var(dt)[0,0])

class Edge:
    def __init__(self, n_from, n_to, e_prob, e_time):
        self.nfrom = n_from
        self.nto = n_to
        self.prob = e_prob
        self.time = e_time

class Node:
    def __init__(self, is_init=False):
        self.index = -1
        self.x = []
        self.density = 0
        self.node_times = []
        self.node_alphas = []
        self.node_betas = []
        self.pself = 1
        self.holding_time = -1

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
        self.delta = 0.01

        self.nodes = []
        self.edges = []
        self.tree = []
        
        #mydict = (key = node, data = [ [in_edge_1, ... ], [ out_edge_1, ...]])
        self.mydict = dict()

        for i in range(self.num_vert):
            n1 = Node()
            self.nodes.append(n1)
        
        # initial variance
        for i in range(self.num_vert):
            n1 = self.nodes[i]
            n1.density = normal_val(n1.x, array(init_state), init_var)
        self.normalize_density()

        self.points = [self.key(mynode.x) for mynode in self.nodes]
        self.tree = kdtree.KDTree(self.points)
        
    def draw_edges(self, n1):
        is_new_node = True
        if n1 in self.mydict:
            is_new_node = False
            self.mydict[n1][0] = []
        else:
            self.mydict[n1] = [[],[]]

        probs = []
        holding_time = get_holding_time(n1.x, 0, self.bowlr)
        n1.holding_time = holding_time
        pself = holding_time/(holding_time + self.delta)
        n1.pself = pself
        mu = n1.x + fdt(n1.x, 0, holding_time)
        var = get_process_var(holding_time)
        
        trans_probs = []
        n1_out_edges = []
        neighbors_index = self.tree.query_ball_point(self.key(n1.x), self.bowlr)
        for n2_index in neighbors_index:
            n2 = self.nodes[n2_index]
            trans_probs.append(normal_val(n2.x, mu, var))
 
        tot_prob = sum(trans_probs)
        trans_probs = trans_probs/tot_prob
        
        runner = 0
        for n2_index in neighbors_index:
            n2 = self.nodes[n2_index]
            etmp = Edge(n1, n2, trans_probs[runner], holding_time)
            n1_out_edges.append(etmp)
            runner = runner + 1

        self.mydict[n1][0] = n1_out_edges

    def draw_reverse_edges(self, n1):
        n1_in_edges = []
        for etmp in self.mydict[n1][0]:
            trans_prob = 1
            for to_edge in self.mydict[etmp.nto][0]:
                if to_edge.nto == n1:
                    trans_prob = to_edge.prob
            enew = Edge(etmp.nto, n1, trans_prob, etmp.nto.holding_time)
            n1_in_edges.append(enew)
        self.mydict[n1][1] = n1_in_edges

    def connect(self):
        for n1 in self.nodes:
            self.draw_edges(n1)
            for etmp in self.mydict[n1][0]:
                if not (etmp.nto in self.mydict):
                    self.draw_edges(etmp.nto)
            self.draw_reverse_edges(n1)
        print "finished connecting"

    def print_graph(self):
        for n1 in self.mydict.keys():
            print n1.x
            for etmp in self.mydict[n1][0]:
                print "\t", etmp.nto.x,'\t', etmp.prob, " ", etmp.time
            print
            for etmp in self.mydict[n1][1]:
                print "\t", etmp.nfrom.x,'\t', etmp.prob, " ", etmp.time
            raw_input()
    
    
    
    """
    def simulate_trajectories(self):
        # utraj is a array of (t, u)
        fig = figure(1)
        trajs = []
        traj_probs = []
        traj_times = []
        max_time = 2
        for traj_index in range(5000):
            traj_curr = []
            curr_time = 0
            traj_time = []
            node_curr = self.nodes[0]
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
                curr_time = curr_time + node_curr.edge_transition_times[0][next_index]
                node_curr = self.mydict[node_curr][next_index]
                # print "changed node to: ", node_curr.x, " time to: ", curr_time

                traj_curr.append(list(node_curr.x))
                traj_time.append(curr_time)

            to_put = [item for sublist in traj_curr for item in sublist]
            while len(to_put) < 100:
                to_put.append(to_put[-1])
                traj_time.append(traj_time[-1])

            trajs.append(to_put)
            traj_probs.append(traj_prob)
            traj_times.append(traj_time)
            
            #plot(traj_time, to_put,'bo-')
            #plot(traj_times,traj_controls,'r--')
        
        #print traj_probs
        #print traj_times
        trajs = array(trajs)
        traj_probs = array(traj_probs)
        traj_avg = average(trajs, axis=0, weights=traj_probs)
        traj_std = array([wstd(trajs[:,i], traj_probs) for i in range(len(trajs[0,:]))])
        
        grid()
        plot(traj_times[0], traj_avg, 'b-', label='mean')
        plot(traj_times[0], traj_avg-traj_std, 'b--', label='+/- std')
        plot(traj_times[0], traj_avg+traj_std, 'b--')
        plot(linspace(0,max_time,1000), exp(-linspace(0,max_time,1000)), 'r-', label='cont. mean')
        xlabel('t [s]')
        ylabel( 'x(t), xd(t)')
        #plot(u_traj[:,0], u_traj[:,1], 'r-')
        show()
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

class Smoother:
    def __init__(self, num_vert):
        self.max_time = 1
        self.times = []
        self.truth = []
        self.observations = []
        self.estimates = []
            
        self.graph = Graph(num_vert)
        self.graph.connect()
        self.propagate_system()

    def propagate_system(self):
        curr_time = 0
        integration_delta = 1e-3
        self.truth.append(init_state)
        self.observations.append( get_observation(init_state, False))
        self.times.append(curr_time)

        while curr_time < self.max_time:
            curr_state = self.truth[-1]
            next_state = curr_state
            next_state = [next_state[i] + (-next_state[i]*integration_delta + 
                                           normal(0, sqrt(get_process_var(integration_delta))[0]) ) for i in range(NUM_DIM)]
            self.truth.append(next_state)
            self.observations.append( get_observation(next_state, False))
            curr_time = curr_time + integration_delta
            self.times.append(curr_time)

        """
        figure(1)
        axis('tight')
        plot(self.times, self.truth, 'r-')
        plot(self.times, self.observations, 'b-')
        #plot(self.times, self.estimates, 'g-')
        grid()
        show()
        """

    def run_smoother(self, num_times):
        graph = self.graph
        self.time_partitions = list(linspace(0., self.max_time, num_times))
        self.alphas = zeros((graph.num_vert, num_times))
        self.betas = zeros((graph.num_vert, num_times))
        
        for n1_index in range(graph.num_vert):
            self.alphas[n1_index,0] = normal_val(graph.nodes[n1_index].x, array(init_state),
                                                 init_var)*normal_val(graph.nodes[n1_index].x, array(self.observations[0]), array(get_observation_var()))
            self.betas[n1_index,-1] = 1
        
        for curr_time in self.time_partitions[1:len(self.time_partitions)-2]:
            


if __name__ == "__main__":
        
    # patch kdtree so that we can pickle it
    kdtree.node = kdtree.KDTree.node
    kdtree.leafnode = kdtree.KDTree.leafnode
    kdtree.innernode = kdtree.KDTree.innernode

    if len(sys.argv) >=3:
        if sys.argv[1] == 'w':
            num_vert = int(sys.argv[2])
            tic = clock()
            smoother = Smoother(num_vert)
            smoother.run_smoother(100)

            # graph.print_graph()
            # pickle.dump(graph, open('graph.pkl','wb'))
            print clock() - tic, '[s]'
    else:
        graph = pickle.load(open('graph.pkl','rb'))

