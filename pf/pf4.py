#!/usr/bin/python

from numpy import *
from numpy.random import *
from pylab import *

sigobs = 1.0

def system(x):
    return x

def obs(x):
    return x

def obsx(x):
    return 1.0

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

def implicit_filter(y, N):

    yhat = zeros(len(y), float)
    yhat[0] = y[0]
    xpart = zeros(N, float)
    w = zeros(N, float)
    for i in range(N):
        xpart[i] = y[0]
        w[i] = 1.0/N

    k = 1
    while k < len(y):
        for i in range(N):
            zeta = normal(0, 1.0)
            xj = 1.
            xj1 = 0.0
            phij = 0.
            sigj = 1.0
            mbarj = 0.
            while(abs(xj1 - xj) > 0.001):
                xj = xj1
                
                sigj = 1/(1. + obsx(xj)*obsx(xj))
                zj = y[k] - obs(xj) + obsx(xj)
                mbarj = sigj*(xpart[i] + system(xpart[i]) + obsx(xj)*zj)
                kj = obsx(xj)*obsx(xj) + 1.
                phij = pow( (zj - obsx(xj)*(xpart[i] + system(xpart[i]))), 2)/kj/2

                xj1 = mbarj + sqrt(sigj)*zeta
                
            J = zeta*sigj/(xj1 - mbarj)
            
            xpart[i] = xj1
            w[i] = exp(-1*phij)*abs(J)

        #xpart, w = resample(xpart, w)
        """ 
        # updated all particles, now resample
        ct = std(w)**2/mean(w)**2
        w /= sum(w)
        if 1/(1 + ct) < 0.75:
            #resample
            xpart, w = resample(xpart, w)
        """
        w /= sum(w)
        yhat[k] = sum(xpart*w)
        #print y[k], yhat[k]
        #temp = raw_input()
        
        k = k+1
    return yhat

if __name__=="__main__":
    
    t = 1000
    dt = 0.001
    x = zeros(t, float)
    y = zeros(t, float)
     
    x[0] = 0.
    y[0] = x[0]
    i = 1
    while i< t:
        x[i] = x[i-1] + system(x[i-1])*dt + sqrt(dt)*normal(0, 1.)
        y[i] = obs(x[i]) + normal(0, sigobs)
        i += 1
        
    yhat = implicit_filter(y, 50)    
    
    plot(x, 'r-')
    plot(yhat, 'g-')
    grid()
    show()
