#!/usr/bin/python

from time import *
from numpy import *
from random import *
from math import *
from scipy.spatial import KDTree
from pylab import *

NUM_DIM=1
zmin = [-0.5]
zmax = [1]
init_state = [1]
init_var = [10e-2]
def fdt(x, dt):
    return -x*dt
def get_process_var(dt):
    return [0.1*dt]
def observation_var():
    return 0.1
def get_holding_time(z,bowlr):
    dt = 1
    if NUM_DIM >=  2:
        return bowlr*bowlr/(bowlr*linalg.norm(fdt(z,dt)) + trace(array(get_process_var(dt))))
    else:
        return bowlr*bowlr/(bowlr*linalg.norm(fdt(z,dt)) + get_process_var(dt)[0])

def normal_val(x, mu, sigma):
    dim = len(mu)
    delx = x -mu;
    var = diag(sigma)*diag(sigma)
    det = sqrt(linalg.det(var))
    return 1/pow(2*pi,dim/2.0)/det*exp(-0.5*dot(delx,dot(linalg.inv(var),delx)))

class Node:
    x = []
    density = 0
    edge_probs = []
    edge_times = []

    def __init__(self, is_init=False):
        if is_init == False:
            self.x = array([zmin[i] + random()*(zmax[i]-zmin[i]) for i in range(NUM_DIM)])
        else:
            self.x = array(init_state)
        # self.density = normal_val(self.x, array(init_state), array(init_var))
        # print self.x
    
class Graph:

    num_vert  = 50
    bowlr = 2*pow(log(num_vert)/float(num_vert),1/float(NUM_DIM))
    delta = 0.1

    nodes = []
    point = []
    
    mydict = []

    def key(self, z):
        return [(z[i] - zmin[i])/(zmax[i]-zmin[i]) for i in range(NUM_DIM)]

    def __init__(self):
        self.nodes.append(Node(True))
        for i in range(self.num_vert-1):
            self.nodes.append(Node())
        self.points = [self.key(mynode.x) for mynode in self.nodes]
        
    def draw_edges(self, n1):
        # n1 is the key of mydict
        probs = []
        holding_time = get_holding_time(n1.x[1:], self.bowlr)
        for n2 in self.mydict[n1]:
            mu = n1.x + fdt(n1.x, holding_time)
            probs.append(normal_val(n2.x, mu, 
                            get_process_var(holding_time)))
        tot_prob = sum(probs)
        probs = probs/tot_prob
        
        # add self-transition and populate edge_probs
        self.mydict[n1].append(n1)
        pself = holding_time/(holding_time + self.delta)
        
        new_probs = [probs[i]*(1-pself) for i in range(len(probs))]
        new_probs.append(pself)
        new_holding_time = self.delta*holding_time/(self.delta+holding_time)
        
        n1.edge_probs = new_probs
        n1.edge_times = [new_holding_time for i in range(len(new_probs))]

    def connect(self):
        for n1 in self.nodes:
            neighbors = []
            for n2 in self.nodes:
                if n1 != n2:
                    if( linalg.norm(n1.x - n2.x) <= self.bowlr):
                        neighbors.append(n2)
            self.mydict.append((n1, neighbors))
        
        self.mydict = dict(self.mydict)
        
        for n1 in self.mydict.keys():
            self.draw_edges(n1)
        """ 
        for n1 in self.mydict.keys():
            print n1.x
            count = 0
            for n2 in self.mydict[n1]:
                print "\t", n2.x, " ", n1.edge_probs[count]
                count = count + 1
        """

    def simulate_trajectories(self):
        fig = figure(1)
        trajs = []
        for traj_index in range(10):
            traj_curr = []
            traj_time = 0
            node_curr = self.nodes[0]
            while traj_time < 100:
                cum_probs = cumsum(node_curr.edge_probs)
                coin_toss = random()
                self_index = len(self.mydict[node_curr]) -1
                next_index = 0
                for i in range(len(cum_probs)):
                    if coin_toss > cum_probs[i]:
                        next_index = next_index+1
                    else:
                        break
                
                #print len(self.mydict[node_curr]), next_index
                traj_curr.append(list(node_curr.x))
                node_curr = self.mydict[node_curr][next_index]
                if next_index == self_index:
                    traj_time = traj_time + self.delta
            
            to_put = [item for sublist in traj_curr for item in sublist]
            trajs.append(to_put)
            plot(to_put,'b-')
        
        grid()
        show()

if __name__ == "__main__":
    seed(0)
    
    tic = clock()
    graph = Graph()
    graph.connect()
    graph.simulate_trajectories()

    print clock() - tic, '[s]'
