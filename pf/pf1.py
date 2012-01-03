#!/usr/bin/python

from numpy import *
from numpy.random import *
from pylab import *

which_noise = 0
sigma_obs = 1
def system(x):
    if which_noise == 0:
        return 1 + sin(1*pi*x) + x/2 + 0.01*normal(1*fabs(x), 1*x*x)
    else:
        return 1 + sin(1*pi*x) + x/2 + 0.01*gamma(1,fabs(x))

def normal_prob(predict, reading, sigma):
    return (1 /sqrt(sigma*sigma*2*pi))* exp((-(reading-predict)**2)/2/sigma/sigma)

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
    
    # init N particles
    for i in range(N):
        #xpart[i] = (random() - 0.5)*100
        xpart[i] = x0 + normal(0, 10)
        w[i] = 1.0/N
    #print xpart

    for k in range(len(y)):
        for i in range(N):
            # predict
            xpart[i] = system(xpart[i])
            # update weights
            w[i] = w[i]*normal_prob(xpart[i], y[k], sigma_obs)

        ct = std(w)**2/mean(w)**2
        w /= sum(w)     
        if 1/(1 + ct) < 0.75:
            #resample
            xpart, w = resample(xpart, w)
        
        xhat[k] = sum(xpart*w)

    return xhat


def run_pf():
    global which_noise

    t = 50
    # seed(0)
    x = zeros(t, float)
    x[0] = 1.0
    y = zeros(t, float)
    y[0] = 1.0
    
    which_noise = 1
    for i in range(len(x) - 1):
        x[i+1] = system(x[i])
        tmp = x[i+1]
        y[i+1] = tmp + normal(0, sigma_obs)
    
    which_noise = 0
    xhat = particle_filter(y, 200, 1.0)
    print "normal -- mean: ", mean(x-xhat), " stddev: ", std(x-xhat)
    
    which_noise = 1
    xhat = particle_filter(y, 200, 1.0)
    print "gamma -- mean: ", mean(x-xhat), " stddev: ", std(x-xhat)
    
    print "max x: ", max(x)

if __name__=="__main__":
    
    for i in range(10):
        run_pf()
        print

    """
    plot(x, 'r-')
    #plot(y, 'b-')
    plot(xhat, 'g-')
    #legend( ('System', 'Obs', 'Filtered'), loc='best')
    grid()
    show()
    """
