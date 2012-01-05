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

init_state = 0.0
init_var = 0.001
process_var = 0.001
observation_var = 0.001
h = 0.01
delta = 0.001
zmin = 0
zmax = 1
states = linspace(zmin, zmax+h, (zmax-zmin)/h)
edge_times = []
edge_probs = []     # (left, self, right)
state_density = []
def drift(s):
    return cos(2*pi*s)

for s in states:
    c = process_var*process_var + h*drift(s)
    htime = h*h/c
    pself = htime/(htime+delta)
    state_density.append( normal_val(s, init_state, init_var) )
    edge_times.append(htime*delta/(htime+delta))
    if (s < zmax) and (s > zmin):
        edge_probs.append([(process_var*process_var/2 + h*drift(s))/c*(1-pself), pself, process_var*process_var/2/c*(1-pself)])
    elif (s >= zmax):
        edge_probs.append([(1-pself), pself, 0.])
    elif (s <= zmin):
        edge_probs.append([0., pself, 1-pself])

def get_observation(x):
    return x + normal(0, sqrt(observation_var))
def wstd(arrin, weights_in):
    arr = array(arrin, ndmin=1, copy=False)
    weights = array(weights_in, ndmin=1, dtype='f8', copy=False)
  
    wtot = weights.sum()
    wmean = ( weights*arr ).sum()/wtot

    wvar = ( weights*(arr-wmean)**2 ).sum()/wtot
    wsdev = sqrt(wvar)
    return wmean,wsdev
def normalize_density():
    global state_density
    state_density = state_density/sum(array(state_density))
def simulate_trajectories():
    figure(1)
    trajs = []
    traj_times = []
    traj_probs = []
    traj_states = []
    max_time = 1
    for i in range(10):
        traj_state = []
        traj_time = []
        prob = 1
        curr_index = find_closest_index(list(states),init_state)
        curr_time = 0
        traj_state.append(states[curr_index])
        traj_time.append(curr_time)

        while curr_time < max_time:
            cointoss = random()
            curr_time = curr_time + edge_times[curr_index]
            prob_sum = cumsum(edge_probs[curr_index])
            if (cointoss < prob_sum[0]):
                prob = prob*edge_probs[curr_index][0]
                curr_index = curr_index - 1
            elif (cointoss < prob_sum[1]):
                prob = prob*edge_probs[curr_index][1]
                curr_index = curr_index
            else:
                prob = prob*edge_probs[curr_index][2]
                curr_index = curr_index +1

            traj_state.append(states[curr_index])
            traj_time.append(curr_time)
        
        traj_probs.append(prob)
        traj_states.append(traj_state)
        traj_times.append(traj_time)
        # plot(traj_time, traj_state, 'b--')
    
    traj_states = array(traj_states)
    
    traj_state_avg = array([wstd(traj_states[:,i],traj_probs)[0] for i in range(len(traj_states[0,:]))])
    traj_state_std = array([wstd(traj_states[:,i],traj_probs)[1] for i in range(len(traj_states[0,:]))])
    plot(traj_times[0], traj_state_avg, 'b-')
    plot(traj_times[0], traj_state_avg+traj_state_std, 'b--')
    plot(traj_times[0], traj_state_avg-traj_state_std, 'b--')
    plot(linspace(0,max_time,1000), init_state*exp(-linspace(0,max_time,1000)), 'r-', label='cont. mean')
    grid()
    show()

def update_conditional_density(curr_observation):
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

def run_filter():
    global state_density

    normalize_density()
    truth, observations, estimates = [], [], []
    max_time = 1
    curr_time = 0
    integration_delta = 1e-3

    truth.append(init_state)
    mean_std = wstd(states, state_density)
    estimates.append(mean_std[0])
    observations.append( get_observation(init_state))
    while curr_time < max_time:
        curr_state = truth[-1]
        next_state = curr_state

        runner_time = 0
        while runner_time < delta:
            next_state = next_state + (drift(next_state)*integration_delta + normal(0, sqrt(process_var*integration_delta)) )
            runner_time = runner_time + integration_delta

        truth.append(next_state)
        observations.append( get_observation(next_state))

        update_conditional_density(observations[-1])

        mean_std = wstd(states, state_density)
        estimates.append(mean_std[0])

        curr_time = curr_time + delta

    figure(1)
    plot(truth, 'r-')
    # plot(observations, 'b-')
    plot(estimates, 'g-')
    grid()
    show()

# simulate_trajectories()
run_filter()
