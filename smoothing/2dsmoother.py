#!/usr/bin/python

"""
2D linear system smoother
"""

import sys
from math import *
import numpy as np
from pylab import *
from scipy.spatial import kdtree
import cPickle as pickle

NUM_DIM = 2
init_state = array([0.5, 0.5])
init_var = array([[1e-3,0.],[0., 1e-3]])
process_var = array([[1e-2, 0.],[0., 1e-2]])
observation_var = array([[1e-3, 0.],[0., 1e-3]])
zmin = [0.3, 0.3]
zmax = [0.6, 0.6]

A = array([[-1.,0.],[0.,-1.]])
C = array([[1.,0.],[0.,1.]])

def find_closest_index(mylist, myvar):
    tmp = [ abs(mylist[i] - myvar) for i in range(len(mylist))]
    return tmp.index(min(tmp))
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
def calci_moment(arrin, weights_in, m):
    arr = np.array(arrin, ndmin=1, copy=False)
    weights = np.array(weights_in, ndmin=1, dtype='f8', copy=False)
    nweights = weights/weights.sum()
    return ( nweights*(arr**m) ).sum()
def normalize_density(arr):
    arr = arr/sum(arr)

# system
def drift(z, dt=1):
    return dot(A,z)*dt
def get_observation(z):
    y = dot(C,z) + np.random.multivariate_normal([0,0], sqrt(observation_var))
    return y

class Node:
    def get_key(self, z):
        return [(z[i] - zmin[i])/(zmax[i]-zmin[i]) for i in range(NUM_DIM)]
    def __init__(self, index_in):
        self.z = array([zmin[i] + (zmax[i]-zmin[i])*np.random.rand() for i in range(NUM_DIM)])
        self.index = index_in
        self.key = self.get_key(self.z)
        self.htime = -1
    def get_htime(self, bowlr):
        h = bowlr*(zmax[0] - zmin[0])
        self.htime = h*h/(process_var[0,0] + h*np.linalg.norm(drift(self.z)))

class Approximate_chain:

    def __init__(self, num):
        self.num_vert = num
        self.bowlr = 2.1*pow(log(self.num_vert)/float(self.num_vert),1/float(NUM_DIM))
        self.nodes = []
        self.points = []
        node_tree = []
        for i in range(self.num_vert):
            n1 = Node(i)
            self.nodes.append(n1)
            self.points.append(n1.key)
            n1.get_htime(self.bowlr)
        self.node_tree = kdtree.KDTree(self.points)
        
        """
        self.P = zeros((self.num_vert, self.num_vert))
        self.Phdelta = zeros((self.num_vert, self.num_vert))
        for i in range(self.num_vert):
            n1 = self.nodes[i]
            mu = n1.z + drift(n1.z, n1.htime)
            var = process_var*n1.htime
            n1_neighbors = self.node_tree.query_ball_point(n1.key, self.bowlr)
            probs = []
            for n2_i in n1_neighbors:
                cprob = normal_val(self.nodes[n2_i].z, mu, var)
                probs.append(cprob)
            probs = array(probs)/sum(probs)
            count = 0
            for n2_i in n1_neighbors:
                self.P[n1.index][n2_i] = probs[count]
                count = count + 1
        """
        self.delta = 0.001  #0.89*min([n1.htime for n1 in self.nodes])
        
        """
        # get filtering one_step transition probabilities using matrices
        self.delta = min(edge_times)*0.999
        P1 = np.zeros((self.tot_vert, self.tot_vert))
        P0 = np.zeros((self.tot_vert, self.tot_vert))
        for i in range(self.tot_vert):
            pt = edge_times[i]/(self.delta + edge_times[i])
            if( (i!=0) and (i!=self.tot_vert-1)):
                P1[i,i-1] = edge_probs[i][0]*pt
                P1[i,i+1] = edge_probs[i][1]*pt
                P0[i,i-1] = edge_probs[i][0]*(1-pt)
                P0[i,i+1] = edge_probs[i][1]*(1-pt)
            elif(i ==0):
                P1[i,i+1] = edge_probs[i][1]*pt
                P0[i,i+1] = edge_probs[i][1]*(1-pt)
            else:
                P1[i,i-1] = edge_probs[i][0]*pt
                P0[i,i-1] = edge_probs[i][0]*(1-pt)
        Phdelta = np.linalg.inv(eye(self.tot_vert) - P0)*P1
        """
        """ 
        self.delta = min(edge_times)*0.899
        #self.delta = 0.0001
        #print 'min_htime: ', min(edge_times),' delta: ', self.delta
        if(min(edge_times) < self.delta):
            print "Add less nodes"
            sys.exit(0)
        # explicit method
        Phdelta = np.zeros((self.tot_vert, self.tot_vert))
        for i in range(self.tot_vert):
            ps = 1 - self.delta/edge_times[i]
            Phdelta[i,i] = ps 
            if( (i!=0) and (i!=self.tot_vert-1)):
                Phdelta[i,i+1] = edge_probs[i][1]*(1- ps)
                Phdelta[i,i-1] = edge_probs[i][0]*(1- ps)
            elif(i ==0):
                Phdelta[i,i+1] = edge_probs[i][1]*(1- ps)
            else:
                Phdelta[i,i-1] = edge_probs[i][0]*(1- ps)
        """
    """ 
    # use obs_curr on alphas[oi,:] to get alphas[oi+1,:]
    def update_alpha(self, oi, obs_curr, alphas):
        for i in range(len(states)):
            tprob = 0
            if (i!=0) and (i!= (len(states)-1)):
                tprob = alphas[oi,i-1]*Phdelta[i-1,i] + alphas[oi,i+1]*Phdelta[i+1,i] + alphas[oi,i]*Phdelta[i,i]
            elif (i==0):
                tprob = alphas[oi,i+1]*Phdelta[i+1,i] + alphas[oi,i]*Phdelta[i,i]
            else:
                tprob = alphas[oi,i-1]*Phdelta[i-1,i] + alphas[oi,i]*Phdelta[i,i]
            rep = normal_val(states[i], obs_curr, observation_var)
            alphas[oi+1,i] = tprob*rep
        normalize_density(alphas[oi+1,:])
 
    # use obs_curr on betas[oi,:] to get betas[oi-1,:]
    def update_beta(self, oi, obs_curr, betas):
        for i in range(len(states)):
            tprob = 0
            if (i!=0) and (i!= (len(states)-1)):
                tprob = betas[oi,i-1]*Phdelta[i,i-1]*normal_val(states[i-1], obs_curr, observation_var) + betas[oi,i+1]*Phdelta[i,i+1]*normal_val(states[i+1], obs_curr, observation_var)
                + betas[oi,i]*Phdelta[i,i]*normal_val(states[i], obs_curr, observation_var)
            elif (i==0):
                tprob = betas[oi,i+1]*Phdelta[i,i+1]*normal_val(states[i+1], obs_curr, observation_var) + betas[oi,i]*Phdelta[i,i]*normal_val(states[i], obs_curr, observation_var)
            else:
                tprob = betas[oi,i-1]*Phdelta[i,i-1]*normal_val(states[i-1], obs_curr, observation_var) + betas[oi,i]*Phdelta[i,i]*normal_val(states[i], obs_curr, observation_var)
            betas[oi-1,i] = tprob
        normalize_density(betas[oi-1,:])   
    """

    def propagate_system(self, max_time):
        times, truth, observations = [], [], []
        
        np.random.seed(10)
        if(max_time < self.delta):
            print "Increase max_time"
            sys.exit(0)
         
        curr_time = 0
        integration_delta = self.delta/2.0
        truth.append(init_state)
        observations.append(init_state)
        times.append(curr_time)
        while curr_time < max_time:
            next_state = truth[-1]
            runner_time = 0
            while np.fabs(runner_time - self.delta) > integration_delta/2.0:
                next_state = next_state + (drift(next_state)*integration_delta + np.random.multivariate_normal([0,0], process_var*integration_delta) )
                runner_time = runner_time + integration_delta
            
            next_obs_state = next_state + np.random.multivariate_normal([0,0], observation_var)
            
            curr_time = curr_time + self.delta
            truth.append(next_state)
            observations.append(next_obs_state)
            times.append(curr_time)
        
        truth = array(truth)
        observations = array(observations)
        """
        plot(times, truth[:,0], 'r-')
        plot(times, observations[:,0], 'b-')
        grid()
        show()
        """
        return times, truth, observations
    
    def run_kalman_smoother(self, times, truth, observations):
        kfestimates = []
        kfvars = []
        
        kf_var = init_var
        kfestimates.append(init_state)
        kfvars.append(kf_var)

        curr_time = 0
        integration_delta = self.delta/2.0
        oi = 1
        while oi < len(observations):
            estimate_state = kfestimates[-1]
            runner_time = 0
            while np.fabs(runner_time - self.delta) > integration_delta/2.0:
                estimate_state = estimate_state + (drift(estimate_state)*integration_delta + np.random.multivariate_normal([0,0], process_var*integration_delta) )
                runner_time = runner_time + integration_delta
            
            kf_var  = exp(-2*self.delta)*kf_var + process_var*self.delta
            S = observations[oi,:] - estimate_state
            gain = kf_var*np.linalg.inv(kf_var + observation_var)
            kf_var = (np.eye(2)-gain)*kf_var
            
            kfestimates.append(estimate_state + dot(gain,S))
            kfvars.append(kf_var)
            oi = oi + 1
        
        kfestimates = array(kfestimates)
        

        ksestimates = kfestimates.copy()
        backward_kf_vars = arange(len(kfvars)-1,1,-1)
        for oi in backward_kf_vars:
            x_k1_n = ksestimates[oi,:]
            x_k1_k = kfestimates[oi-1,:] + drift(kfestimates[oi-1,:], self.delta)
            kfvar_k1_k = exp(-2*self.delta)*kfvars[oi-1] + process_var*self.delta
            ak = dot(kfvars[oi-1]*exp(-1*self.delta), np.linalg.inv(kfvar_k1_k))
            x_k_n = kfestimates[oi-1] + dot(ak,x_k1_n - x_k1_k)
            ksestimates[oi-1,:] = x_k_n
        
        return kfestimates, ksestimates, kfvars
    
    """
    def run_hmmf_smoother(self):
        
        times, truth, observations = self.propagate_system(0.1)
        kfestimates, ksestimates, kfvars = self.run_kalman_smoother(times, truth, observations)

        alphas = zeros((len(times)+1, self.tot_vert))
        betas = zeros((len(times)+1, self.tot_vert))
        for si in range(len(states)):
            alphas[0,si] = normal_val(states[si], init_state, init_var)
            betas[-1,si] = 1.0
        normalize_density(alphas[0,:])
        
        forward_obs = arange(0, len(observations))
        for oi in forward_obs:
            self.update_alpha(oi, observations[oi], alphas)
        backward_obs = arange(len(observations)-1,-1,-1)
        for oi in backward_obs:
            self.update_beta(oi+1, observations[oi], betas)
        #print "updated both"
        
        sestimates = []
        density = zeros((len(times)+1, self.tot_vert))
        for oi in range(len(times)+1):
            density[oi,:] = alphas[oi,:]*betas[oi,:]
            normalize_density(density[oi,:])
            sestimates.append(calci_moment(states, density[oi,:], 1))
        sestimates.pop()
        
        festimates = []
        for oi in range(len(times)):
            festimates.append(calci_moment(states, alphas[oi+1,:], 1))
        
        figure(1)
        plot(truth, 'r-')
        #plot(observations, 'b-')
        plot(festimates, 'g-')
        plot(sestimates, 'g--')
        plot(kfestimates, 'c-')
        plot(ksestimates, 'c--')
        axis('tight')
        grid()
        show()
        return (np.linalg.norm(np.array(truth)-np.array(kfestimates))**2)/(0.1/self.delta), (np.linalg.norm(np.array(truth)-np.array(ksestimates))**2)/(0.1/self.delta), \
                (np.linalg.norm(np.array(truth)-np.array(festimates))**2)/(0.1/self.delta) ,(np.linalg.norm(np.array(truth)-np.array(sestimates))**2)/(0.1/self.delta)
    """

if __name__ == "__main__":
    
    # patch kdtree for pickling
    kdtree.node = kdtree.KDTree.node
    kdtree.leafnode = kdtree.KDTree.leafnode
    kdtree.innernode = kdtree.KDTree.innernode
    
    """
    if 0:
        for n in linspace(10, 600, 20):
            n = int(n)
            amc = Approximate_chain(n)
            kferr, kserr, ferr, serr = amc.run_hmmf_smoother()
            print n, kferr, kserr, ferr, serr
    else:
        n = int(sys.argv[1])
        amc = Approximate_chain(n)
        kferr, kserr, ferr, serr = amc.run_hmmf_smoother()
        print kferr, kserr, ferr, serr
    """
    
    if 1:
        n = int(sys.argv[1])
        amc = Approximate_chain(n)
        times, truth, observations = amc.propagate_system(0.5)
        amc.run_kalman_smoother(times, truth, observations)
        #pickle.dump(amc, open('amc.pkl','wb'))
    else:
        print "no"
        #amc = pickle.load(open('amc.pkl', 'rb'))
