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

def find_closest_index(mylist, myvar):
    tmp = [ abs(mylist[i] - myvar) for i in range(len(mylist))]
    return tmp.index(min(tmp))
def normal_val(x, mu, var):
    return 1/sqrt(2*pi*var)*exp((-0.5*(x-mu)**2)/var)
def calci_moment(arrin, weights_in, m):
    arr = array(arrin, ndmin=1, copy=False)
    weights = array(weights_in, ndmin=1, dtype='f8', copy=False)
    nweights = weights/weights.sum()
    return ( nweights*(arr**m) ).sum()
def normalize_density():
    global state_density
    state_density = state_density/sum(array(state_density))

init_state = 0.7
init_var = 0.001
process_var = 0.001
observation_var = 0.001
h = -1
zmin = 0.2
zmax = 1
states = linspace(zmin, zmax, (zmax-zmin)/h)
edge_times = []
edge_probs = []     # (left, self, right)
state_density = []

class Approximate_chain:
    tot_vert = 0 
    k = 3.0
    def drift(self,s):
        #return cos(2*pi*s)
        return -self.k*s
    def get_observation(self,x):
        return x + normal(0, sqrt(observation_var))
    def __init__(self, num):
        global h, edge_times, edge_probs, state_density, states
        states, state_density, edge_probs, edge_times = [],[],[], []
        
        self.tot_vert = num
        h = (zmax-zmin)/float(num)
        states = linspace(zmin, zmax, (num+1))
        delta = h*h/(h*fabs(self.drift(zmin)) + process_var)
        for s in states:
            c = process_var + h*fabs(self.drift(s))
            htime = h*h/c
            # print "htime: ", htime
            pself = 0*htime/delta
            state_density.append( normal_val(s, init_state, init_var) )
            edge_times.append(htime*(1-pself))
            if (s < zmax) and (s > zmin):
                edge_probs.append([(process_var/2 + h*fabs(self.drift(s)))/c*(1-pself), pself, process_var/2/c*(1-pself)])
            elif (s >= zmax):
                edge_probs.append([(1-pself), pself, 0.])
            elif (s <= zmin):
                edge_probs.append([0., pself, 1-pself])
    
    def simulate_trajectories(self):
        figure(1)
        trajs = []
        traj_times = []
        traj_probs = []
        traj_states = []
        max_time = 0.1
        for i in range(1000):
            traj_state = []
            traj_time = []
            start_state = normal(init_state, sqrt(init_var))
            curr_index = find_closest_index(list(states), start_state)
            prob = normal_val(init_state, init_var, states[curr_index]) 
            curr_time = 0
            traj_state.append(states[curr_index])
            traj_time.append(curr_time)

            while curr_time < max_time:
                cointoss = random()
                curr_time = curr_time + edge_times[curr_index]
                prob_sum = cumsum(edge_probs[curr_index])
                #print prob_sum
                if (cointoss < prob_sum[0]):
                    prob = prob*edge_probs[curr_index][0]
                    curr_index = curr_index - 1
                elif (cointoss < prob_sum[1]):
                    prob = prob*edge_probs[curr_index][1]
                    curr_index = curr_index
                else:
                    prob = prob*edge_probs[curr_index][2]
                    curr_index = curr_index +1

                #print prob
                traj_state.append(states[curr_index])
                traj_time.append(curr_time)

            if prob < 10e-300:
                print prob

            traj_probs.append(prob)
            traj_states.append(traj_state)
            traj_times.append(traj_time)
            #plot(traj_time, traj_state, 'b--')

            #if (i+1)%1000 == 0:
            #    print i
    
        # print traj_probs
        refine = max_time/10 # min(edge_times)
        time_array = linspace(0,max_time,(max_time+refine)/refine)
        state_array = zeros((len(traj_states), len(time_array)))
        for ti in range(len(time_array)):
            t = time_array[ti]
            for si in range(len(traj_states)):
                state_index = find_closest_index(traj_times[si], t)
                state_array[si,ti] = traj_states[si][state_index]

        state_array_avg = array([calci_moment(state_array[:,i],traj_probs,1) for i in range(len(state_array[0,:]))])
        state_array_m2 = array([calci_moment(state_array[:,i],traj_probs,2) for i in range(len(state_array[0,:]))])
        cont_avg = init_state*exp(-self.k*time_array)
        cont_m2 = init_var/self.k*(1.0 + (self.k-1)*exp(-self.k*time_array)) + cont_avg**2

        """
        plot(time_array, state_array_avg, 'b-')
        plot(time_array, state_array_avg+state_array_std, 'b--')
        plot(time_array, state_array_avg-state_array_std, 'b--')
        plot(time_array, init_state*exp(-k*time_array), 'r-', label='cont. mean')
        grid()
        show()
        """
        print self.tot_vert, fabs(cont_avg[-1] - state_array_avg[-1]), fabs(cont_m2[-1] - state_array_m2[-1])

    def update_conditional_density(self,curr_observation):
        global state_density
        for i in range(len(states)):
            tprob = 0
            if (i!=0) and (i!= (len(states)-1)):
                tprob = state_density[i-1]*edge_probs[i-1][2] + state_density[i+1]*edge_probs[i+1][0]
            elif (i==0):
                tprob = state_density[i+1]*edge_probs[i+1][0]
            else:
                tprob = state_density[i-1]*edge_probs[i-1][2]
            tprob = tprob + (state_density[i])*edge_probs[i][1]
            state_density[i] = tprob*normal_val(states[i], curr_observation, observation_var)
        normalize_density()

    def run_filter(self):
        global state_density
        normalize_density()
        truth, observations, estimates = [], [], []
        max_time = 0.5
        curr_time = 0
        integration_delta = 1e-3

        truth.append(init_state)
        mean_std = wstd(states, state_density)
        estimates.append(mean_std[0])
        observations.append( get_observation(self,init_state))
        while curr_time < max_time:
            curr_state = truth[-1]
            next_state = curr_state
            runner_time = 0
            while runner_time < delta:
                next_state = next_state + (self.drift(next_state)*integration_delta + normal(0, sqrt(process_var*integration_delta)) )
                runner_time = runner_time + integration_delta
            truth.append(next_state)
            observations.append( get_observation(self,next_state))
            update_conditional_density(observations[-1])
            mean_std = wstd(states, state_density)
            estimates.append(mean_std[0])
            curr_time = curr_time + delta
        return norm(array(truth)-array(estimates))

        """
        figure(1)
        plot(truth, 'r-')
        # plot(observations, 'b-')
        plot(estimates, 'g-')
        grid()
        show()
        """


if __name__ == "__main__":
    
    d = 10
    nodes = linspace(d,100,(100)/d)
    #nodes = [800]
    for n in nodes:
        amc = Approximate_chain(n)
        amc.simulate_trajectories()
    
    #err = array([run_filter() for i in range(10)])
    #print mean(err)
