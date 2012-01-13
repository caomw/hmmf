#!/usr/bin/python

import sys
from math import *
import numpy as np
from pylab import *

init_state = 0.5
init_var = 1e-3
process_var = 1e-2
observation_var = 1e-3
h = -1
zmin = 0.0
zmax = 1.0
states = []
edge_times = []
edge_probs = []     # (left, self, right)
Phdelta = []

def find_closest_index(mylist, myvar):
    tmp = [ abs(mylist[i] - myvar) for i in range(len(mylist))]
    return tmp.index(min(tmp))
def normal_val(x, mu, var):
    return 1/sqrt(2*pi*var)*exp((-0.5*(x-mu)**2)/var)
def calci_moment(arrin, weights_in, m):
    arr = np.array(arrin, ndmin=1, copy=False)
    weights = np.array(weights_in, ndmin=1, dtype='f8', copy=False)
    nweights = weights/weights.sum()
    return ( nweights*(arr**m) ).sum()
def normalize_density(arr):
    arr = arr/sum(arr)

class Approximate_chain:
    tot_vert = 0 
    def drift(self,s):
        k = 1
        return -k*s
    def get_observation(self,x):
        y = x + np.random.normal(0, sqrt(observation_var))
        if y > zmax:
            return zmax
        elif y< zmin:
            return zmin
        else:
            return y
    def __init__(self, num):
        global h, edge_times, edge_probs, states, Phdelta
        states, edge_probs, edge_times, Phdelta = [],[],[],[]
        
        h = (zmax-zmin)/float(num)
        self.tot_vert = num+1
        """
        for i in range(self.tot_vert):
            states.append( zmin + np.random.rand()*(zmax - zmin))
        sort(states)
        """
        states = np.linspace(zmin, zmax, self.tot_vert)
        for s in states:
            c = process_var + h*np.fabs(self.drift(s))
            htime = h*h/c
            edge_times.append(htime)
        for i in range(self.tot_vert):
            s = states[i]
            c = process_var + h*np.fabs(self.drift(s))
            if (i != 0) and (i != self.tot_vert-1):
                edge_probs.append([(process_var/2 + h*np.fabs(self.drift(s)))/c, process_var/2/c])
            elif (i == self.tot_vert -1):
                edge_probs.append([1.0, 0.])
            elif (i == 0):
                edge_probs.append([0., 1.0])
        
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
    def simulate_trajectories(self):
        figure(1)
        mc_states = []
        actual_states = []
        max_time = 1.0
        for i in range(50000):
            start_state = init_state #np.random.normal(init_state, sqrt(init_var))
            actual_state = init_state
            curr_index = find_closest_index(list(states), start_state)
            curr_time = 0
            while curr_time < max_time:
                cointoss = np.random.rand()
                actual_state = actual_state + self.drift(actual_state)*edge_times[curr_index] + normal(0, sqrt(process_var*edge_times[curr_index]))
                if(actual_state > zmax):
                    actual_state = zmax
                elif(actual_state < zmin):
                    actual_state = zmin
                curr_time = curr_time + edge_times[curr_index]
                prob_sum = cumsum(edge_probs[curr_index])
                #print prob_sum
                if (cointoss < prob_sum[0]):
                    curr_index = curr_index - 1
                else:
                    curr_index = curr_index +1

            mc_states.append(states[curr_index])
            actual_states.append(actual_state)

        print self.tot_vert, np.fabs(np.mean(mc_states) - np.mean(actual_states))
    
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
                next_state = next_state + (self.drift(next_state)*integration_delta + np.random.normal(0, sqrt(process_var*integration_delta)) )
                if(next_state > zmax):
                    next_state = zmax
                elif(next_state < zmin):
                    next_state = zmin
                runner_time = runner_time + integration_delta
            
            next_obs_state = next_state + normal(0, sqrt(observation_var))
            
            curr_time = curr_time + self.delta
            truth.append(next_state)
            observations.append(next_obs_state)
            times.append(curr_time)

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
                estimate_state = estimate_state + (self.drift(estimate_state)*integration_delta + np.random.normal(0, sqrt(process_var*integration_delta)) )
                runner_time = runner_time + integration_delta
            
            kf_var  = exp(-2*self.delta)*kf_var + process_var*self.delta
            S = observations[oi] - estimate_state
            gain = kf_var/(kf_var + observation_var)
            kf_var = (1-gain)*kf_var

            kfestimates.append( estimate_state + gain*S)
            kfvars.append(kf_var)
            oi = oi + 1
        
        ksestimates = kfestimates[:]
        backward_kf_vars = arange(len(kfvars)-1,1,-1)
        for oi in backward_kf_vars:
            x_k1_n = ksestimates[oi]
            x_k1_k = kfestimates[oi-1] + self.drift(kfestimates[oi-1])*self.delta
            kfvar_k1_k = exp(-2*self.delta)*kfvars[oi-1] + process_var*self.delta
            ak = kfvars[oi-1]*exp(-1*self.delta)/kfvar_k1_k
            x_k_n = kfestimates[oi-1] + ak*(x_k1_n - x_k1_k)
            ksestimates[oi-1] = x_k_n
        return kfestimates, ksestimates, kfvars

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
        
        """
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
        """
        return (np.linalg.norm(np.array(truth)-np.array(kfestimates))**2)/(0.1/self.delta), (np.linalg.norm(np.array(truth)-np.array(ksestimates))**2)/(0.1/self.delta), (np.linalg.norm(np.array(truth)-np.array(festimates))**2)/(0.1/self.delta) ,(np.linalg.norm(np.array(truth)-np.array(sestimates))**2)/(0.1/self.delta)

if __name__ == "__main__":
    
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

