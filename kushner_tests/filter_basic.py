#!/usr/bin/python

import sys
from math import *
import numpy as np
from pylab import *

init_state = 0.5
init_var = 0.001
process_var = 0.0001
observation_var = 0.0001
h = -1
zmin = 0.48
zmax = 0.52
states = []
edge_times = []
edge_probs = []     # (left, self, right)
state_density = []
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
def normalize_density():
    global state_density
    state_density = state_density/sum(np.array(state_density))

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
        global h, edge_times, edge_probs, state_density, states, Phdelta
        states, state_density, edge_probs, edge_times, Phdelta = [],[],[],[],[]
        
        h = (zmax-zmin)/float(num)
        self.tot_vert = num+1
        states = np.linspace(zmin, zmax, self.tot_vert)
        for s in states:
            c = process_var + h*np.fabs(self.drift(s))
            htime = h*h/c
            edge_times.append(htime)
            state_density.append( normal_val(s, init_state, init_var) )
        for i in range(self.tot_vert):
            s = states[i]
            c = process_var + h*np.fabs(self.drift(s))
            if (i != 0) and (i != self.tot_vert-1):
                edge_probs.append([(process_var/2 + h*np.fabs(self.drift(s)))/c, process_var/2/c])
            elif (i == self.tot_vert -1):
                edge_probs.append([1.0, 0.])
            elif (i == 0):
                edge_probs.append([0., 1.0])
        
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

        self.delta = min(edge_times)*0.999
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

    def update_conditional_density(self, obs_delta, obs_curr):
        global state_density, Phdelta
        for i in range(len(states)):
            tprob = 0
            if (i!=0) and (i!= (len(states)-1)):
                tprob = state_density[i-1]*Phdelta[i-1,i] + state_density[i+1]*Phdelta[i+1,i] + state_density[i]*Phdelta[i,i]
            elif (i==0):
                tprob = state_density[i+1]*Phdelta[i+1,i] + state_density[i]*Phdelta[i,i]
            else:
                tprob = state_density[i-1]*Phdelta[i-1,i] + state_density[i]*Phdelta[i,i]
            
            rep = exp(states[i]*obs_delta/observation_var - 0.5*states[i]*states[i]*self.delta/observation_var)
            #rep = normal_val(states[i], obs_curr, observation_var)
            state_density[i] = tprob*rep
        normalize_density()

    def run_filter(self):
        global state_density, Phdelta
        normalize_density()
        truth, observations, estimates, kfestimates = [], [], [], []
        max_time = 0.01
        curr_time = 0
        integration_delta = self.delta/2.0
        
        kf_var = init_var
        np.random.seed(10)
        truth.append(init_state)
        mean = calci_moment(states, state_density,1)
        estimates.append(mean)
        kfestimates.append(init_state)
        observations.append(init_state)
        while curr_time < max_time:
            next_state = truth[-1]
            next_obs_state = observations[-1]
            runner_time = 0
            while np.fabs(runner_time - self.delta) > integration_delta/2.0:
                next_state = next_state + (self.drift(next_state)*integration_delta + np.random.normal(0, sqrt(process_var*integration_delta)) )
                if(next_state > zmax):
                    next_state = zmax
                elif(next_state < zmin):
                    next_state = zmin
                next_obs_state = next_obs_state + (next_state*integration_delta + np.random.normal(0, sqrt(observation_var*integration_delta)) )
                if(next_obs_state > zmax):
                    next_obs_state = zmax
                elif(next_obs_state < zmin):
                    next_obs_state = zmin    
                runner_time = runner_time + integration_delta
            
            kf_var  = exp(-2*self.delta)*kf_var + process_var*self.delta
            S = next_obs_state - (observations[-1] + kfestimates[-1]*self.delta)        # this is the next observation, noiseless (ydot = x)
            gain = kf_var*self.delta/(kf_var*self.delta**2 + observation_var*self.delta)
            kfestimates.append( kfestimates[-1] + self.drift(kfestimates[-1])*self.delta + gain*S)
            kf_var = (1-gain*self.delta)*kf_var

            truth.append(next_state)
            observations.append( next_obs_state)
            self.update_conditional_density(observations[-1] - observations[-2], observations[-1])
            mean = calci_moment(states, state_density,1)
            estimates.append(mean)
            curr_time = curr_time + self.delta
        """
        figure(1)
        plot(truth, 'r-')
        plot(observations, 'b-')
        plot(estimates, 'g-')
        plot(kfestimates, 'c-')
        grid()
        show()
        """
        return fabs(np.linalg.norm(np.array(truth)-np.array(estimates))**2 - np.linalg.norm(np.array(truth)-np.array(kfestimates))**2)


if __name__ == "__main__":
    
    """
    d = 10
    nodes = linspace(d,100,(100)/d)
    #nodes = [800]
    for n in nodes:
        amc = Approximate_chain(int(n))
        amc.simulate_trajectories()
    """

    
    np.random.seed(10)
    first, last, step = 1000, 1001, 1
    if(len(sys.argv) > 1):
        first = int(sys.argv[1])
        last = int(sys.argv[2])
        step = int(sys.argv[3])
    nodes = np.arange(first, last, step)
    for n in nodes:
        err = []
        for tries in range(1):
            amc = Approximate_chain(n)
            cerr = amc.run_filter()
            #print n,cerr
            err.append(cerr)
        print n, np.mean(err)
    
