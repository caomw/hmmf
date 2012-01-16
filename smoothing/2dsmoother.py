#!/usr/bin/python

"""
Smoothing using Markov chains
"""

import sys
from math import *
import numpy as np
from pylab import *
from scipy.spatial import kdtree
import cPickle as pickle
import scipy.linalg as slinalg
from matplotlib.font_manager import fontManager, FontProperties

font = FontProperties(size='medium')

NUM_DIM = 2
init_state = array([0.5, 0.5])
init_var = array([[1e-3,0.],[0., 1e-3]])
process_var = array([[1e-3, 0.],[0., 1e-3]])
observation_var = array([[1e-3, 0.],[0., 1e-3]])
zmin = [0.4, 0.4]
zmax = [0.9, 0.6]


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
    return dot(nweights, (arr**m))
def normalize_density(arr):
    arr = arr/sum(arr)

# vanderpol
def A(z=1):
    return array([[0., 1],[-1 - 4.0*z[0]*z[1], 2.0]])
def C(z=1):
    return array([0.,1.0])
def drift(z, dt=1):
    a = array([0., 0.])
    a[0] = z[1]
    a[1] = -z[0] + 2.0*z[1]*(1-z[0]*z[0])
    return a*dt
def get_observation(z):
    noise = np.random.multivariate_normal([0,0], sqrt(observation_var))
    y = z[0] + noise[0] 
    return y

"""
# linear
def A(z=1):
    return array([[-1., 0],[0, -1.0]])
def C(z=1):
    return array([[1.,0.0],[0.,1.0]])
def drift(z, dt=1):
    a = array([0., 0.])
    a[0] = -z[0]
    a[1] = -z[1]
    return a*dt
def get_observation(z):
    noise = np.random.multivariate_normal([0,0], sqrt(observation_var))
    y = z + noise
    return y
"""

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

    def __init__(self, num, to_dump):
        
        np.random.seed(10)
        self.num_vert = num
        self.bowlr = 2.2*pow(log(self.num_vert)/float(self.num_vert),1.0/float(NUM_DIM))
        self.nodes = []
        self.points = []
        node_tree = []
        self.states = []
        for i in range(self.num_vert):
            n1 = Node(i)
            self.states.append(n1.z)
            self.nodes.append(n1)
            self.points.append(n1.key)
            n1.get_htime(self.bowlr)
        self.node_tree = kdtree.KDTree(self.points)
        self.states = array(self.states)

        P = zeros((self.num_vert, self.num_vert))
        self.Phdelta = zeros((self.num_vert, self.num_vert))
        if to_dump:
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
                    P[n1.index][n2_i] = probs[count]
                    count = count + 1
            pickle.dump(P, open('p.pkl', 'wb'))
        else:
            P = pickle.load(open('p.pkl', 'rb'))

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
       
        min_htime = min([n1.htime for n1 in self.nodes])
        self.delta = 0.1*min_htime
        print 'min_htime: ', min_htime,' delta: ', self.delta
        # explicit method
        self.Phdelta = P.copy()
        for i in range(self.num_vert):
            ps = 1 - self.delta/self.nodes[i].htime
            self.Phdelta[i,:] = (1-ps)*P[i,:]
            self.Phdelta[i,i] = ps
                
    # use obs_curr on alphas[oi,:] to get alphas[oi+1,:]
    def update_alpha(self, oi, obs_curr, alphas):
        nodes = self.nodes
        Phdelta = self.Phdelta
        alphas[oi+1,:] = dot(alphas[oi,:], Phdelta)
        for n1 in nodes:
            rep = normal_val(n1.z, obs_curr, observation_var)
            alphas[oi+1,n1.index] = alphas[oi+1,n1.index]*rep
        normalize_density(alphas[oi+1,:])
 
    # use obs_curr on betas[oi,:] to get betas[oi-1,:]
    def update_beta(self, oi, obs_curr, betas):
        nodes = self.nodes
        Phdelta = self.Phdelta
        betac = betas[oi,:].copy()
        for n1 in nodes:
            betac[n1.index] = betac[n1.index]*normal_val(n1.z, obs_curr, observation_var)
            betas[oi-1,:] = dot(Phdelta, betac)
        normalize_density(betas[oi-1,:])

    def propagate_system(self, max_time):
        times, truth, observations = [], [], []
        
        np.random.seed(10)
        if(max_time < self.delta):
            print "Increase max_time", " max_time: ", max_time, " delta: ", self.delta
            sys.exit(0)
         
        curr_time = 0
        integration_delta = min(0.001, self.delta/2.0)
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
        return times, truth, observations
    
    def run_kalman_smoother(self, times, truth, observations):
        
        kfestimates = []
        kfvars = []
        
        kf_var = init_var
        kfestimates.append(init_state)
        kfvars.append(kf_var)

        curr_time = 0
        integration_delta = min(0.001, self.delta/2.0)
        oi = 1
        while oi < len(observations):
            estimate_state = kfestimates[-1]
            runner_time = 0
            while np.fabs(runner_time - self.delta) > integration_delta/2.0:
                estimate_state = estimate_state + (drift(estimate_state)*integration_delta + np.random.multivariate_normal([0,0], process_var*integration_delta) )
                runner_time = runner_time + integration_delta
            
            expA = slinalg.expm(A(kfestimates[-1])*self.delta)
            kf_var  = dot(dot(expA,kf_var),expA.T) + process_var*self.delta
            S = observations[oi,:] - estimate_state
            gain = dot(dot(kf_var, C(estimate_state).T), np.linalg.inv( dot(dot(C(estimate_state),kf_var),C(estimate_state).T) + observation_var))
            kf_var = dot(1.0 -dot(gain, C(estimate_state)), kf_var)
            
            kfestimates.append(estimate_state + dot(gain,S))
            kfvars.append(kf_var)
            oi = oi + 1
        
        kfestimates = array(kfestimates)
        

        ksestimates = kfestimates.copy()
        backward_kf_vars = arange(len(kfvars)-1,1,-1)
        for oi in backward_kf_vars:
            x_k1_n = ksestimates[oi,:]
            x_k1_k = kfestimates[oi-1,:] + drift(kfestimates[oi-1,:], self.delta)
            expA = slinalg.expm(A(kfestimates[oi-1,:])*self.delta)
            kfvar_k1_k = dot(dot(expA,kfvars[oi-1]),expA.T) + process_var*self.delta
            ak = dot(dot(kfvars[oi-1], expA.T), np.linalg.inv(kfvar_k1_k))
            x_k_n = kfestimates[oi-1] + dot(ak,x_k1_n - x_k1_k)
            ksestimates[oi-1,:] = x_k_n
        
        return kfestimates, ksestimates, kfvars
    
    def run_hmmf_smoother(self, max_time):
        nodes = self.nodes

        times, truth, observations = self.propagate_system(max_time)
        kfestimates, ksestimates, kfvars = self.run_kalman_smoother(times, truth, observations)

        sestimates = []
        alphas = zeros((len(times)+1, self.num_vert))
        betas = zeros((len(times)+1, self.num_vert))
        for n1 in nodes:
            alphas[0,n1.index] = normal_val(n1.z, init_state, init_var)
            betas[-1,n1.index] = 1.0
        normalize_density(alphas[0,:])
        
        forward_obs = arange(0, len(observations))
        for oi in forward_obs:
            self.update_alpha(oi, observations[oi], alphas)
        backward_obs = arange(len(observations)-1,-1,-1)
        for oi in backward_obs:
            self.update_beta(oi+1, observations[oi], betas)
        #print "updated both"
        
        density = zeros((len(times)+1, self.num_vert))
        for oi in range(len(times)+1):
            density[oi,:] = alphas[oi,:]*betas[oi,:]
            normalize_density(density[oi,:])
            sestimates.append(list(calci_moment(self.states, density[oi,:], 1)))
        sestimates.pop()
        festimates = []
        for oi in range(len(times)):
            festimates.append(list(calci_moment(self.states, alphas[oi+1,:], 1)))
        
        sestimates = array(sestimates)
        festimates = array(festimates)
        
        fig = figure(1)
        ax = fig.add_subplot(111, aspect='equal')
        plot(times, truth[:,0], 'r-', label='sys')
        #plot(times, observations, 'b-')
        plot(times, festimates[:,0], 'g--', label='hfilter')
        plot(times, sestimates[:,0], 'g-', label='hsmoothing')
        plot(times, kfestimates[:,0], 'c--', label='kfilter')
        plot(times, ksestimates[:,0], 'c-', label='ksmoothing')
        axis('tight')
        grid()
        xlabel('t [s]')
        legend(loc=2, prop=font)
        title('vanderpol_x')
        savefig('smooth_vanderpol_x'+str(self.num_vert)+'.pdf', bbox_inches='tight')
        
        fig = figure(2)
        ax = fig.add_subplot(111, aspect='equal')
        plot(times, truth[:,1], 'r-', label='sys')
        #plot(times, observations, 'b-')
        plot(times, festimates[:,1], 'g--', label='hfilter')
        plot(times, sestimates[:,1], 'g-', label='hsmoothing')
        plot(times, kfestimates[:,1], 'c--', label='kfilter')
        plot(times, ksestimates[:,1], 'c-', label='ksmoothing')
        axis('tight')
        grid()
        xlabel('t [s]')
        legend(loc=3, prop=font)
        title('vanderpol_x_dot')
        savefig('smooth_vanderpol_x_dot_'+str(self.num_vert)+'.pdf', bbox_inches='tight')
        
        show()

        return (np.linalg.norm(np.array(truth)-np.array(kfestimates))**2)/(max_time/self.delta), (np.linalg.norm(np.array(truth)-np.array(ksestimates))**2)/(max_time/self.delta), \
                (np.linalg.norm(np.array(truth)-np.array(festimates))**2)/(max_time/self.delta), 0#,(np.linalg.norm(np.array(truth)-np.array(sestimates))**2)/(max_time/self.delta)

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
        to_dump = int(sys.argv[2])
        if to_dump >= 1:
            amc = Approximate_chain(n, True)
            sys.exit(0)
        else:
            amc = Approximate_chain(n, False)
            kferr, kserr, ferr, serr = amc.run_hmmf_smoother(0.5)
            print kferr, kserr, ferr, serr
    else:
        print "no"
