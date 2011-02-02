#!/usr/bin/python

from numpy import *
from numpy.random import *
from pylab import *

sigma_process = 1
def system(x):
    #return x/2.0 + 25*x/(1 + x*x) + 1.2*cos(x)
    return x/(1 + pow(x, 2))

sigma_obs = 1.0
def obs_model(predict, reading):
    return (1 /sqrt(sigma_obs*2*pi))* exp((-(reading-predict)**2)/2/sigma_obs)

def resample(x, w):
    N = len(x)
    new_x = empty(N)
    new_w = ones(N, dtype=float)/N
    
    c = cumsum(w)
    u = rand()/float(N)
    i = 0
    for j in range(N):
        uj = u + j/float(N)
        while uj > c[i]:
            i += 1
        new_x[j] = x[i]
    
    return new_x, new_w

def particle_filter(y, N, init):
    x0 = init
    xpart = zeros(N, float)
    w = zeros(N, float)
    xhat = zeros(len(y), float)
    xhat[0] = x0

    for i in range(N):                                  # init N particles
        xpart[i] = x0 + normal(0, 10)
        w[i] = 1.0/N

    for k in range(len(y)):
        if k > 0:
            for i in range(N):
                xpart[i] = system(xpart[i])                 # predict
                w[i] = w[i]*obs_model(xpart[i], y[k])       # calci weights
            
            if 1/sum(w) < 2*N/3:
                #resample
                w /= sum(w)     
                xpart, w = resample(xpart, w)    
            xhat[k] = sum(xpart*w)
    
    return xhat

if __name__=="__main__":

    t = 50
    seed(0)
    x = zeros(t, float)
    x[0] = 1.0
    y = zeros(t, float)
    y[0] = 1.0

    for i in range(len(x) - 1):
        x[i+1] = system(x[i])
        tmp = x[i+1] + random()             # process noise
        y[i+1] = tmp*tmp/2 + 2*random()     # obs noise

    xhat = particle_filter(y, 50, 1.0)
    
    plot(x, 'r-')
    plot(y, 'b-')
    plot(xhat, 'g-')
    legend( ('System', 'Obs', 'Filtered'), loc='best')
    grid()
    show()
