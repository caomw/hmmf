#!/usr/bin/python

from numpy import *
from numpy.random import *
from pylab import *

dt = 0.01
Q = 0.01

def system(x):
    return -10*x*(x*x - 1)

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
            xj = 0.5
            xj1 = 0
            phij = 0.
            sigj = 0.0
            mbarj = 0.
            Fn = system(xpart[i])*dt
            Gn2 = dt
            Q2 = Q*Q

            while(abs(xj1 - xj) > 0.001):
                xj = xj1
                
                sigj = 1/(1/Gn2 + obsx(xj)/Q2*obsx(xj))
                zj = y[k] - obs(xj) + obsx(xj)*xj
                mbarj = sigj*( 1/Gn2*(xpart[i] + Fn) + obsx(xj)*zj/Q2)
                kj = Gn2*obsx(xj)*obsx(xj) + Q2
                phij = pow( (zj - obsx(xj)*(xpart[i] + Fn)), 2)/kj/2

                xj1 = mbarj + sqrt(sigj)*zeta
                
                #print y[k], xj1, sigj, mbarj
                #raw_input()

            J = zeta*sigj/(xj1 - mbarj)
            
            xpart[i] = xj1
            w[i] = exp(-1*phij)*abs(J)
        
        #print
        w /= sum(w)
        xpart, w = resample(xpart, w)
        """ 
        # updated all particles, now resample
        ct = std(w)**2/mean(w)**2
        w /= sum(w)
        if 1/(1 + ct) < 0.75:
            #resample
            xpart, w = resample(xpart, w)
        """
        yhat[k] = sum(xpart*w)
        k = k+1
    return yhat

if __name__=="__main__":
    
    t = 100
    x = zeros(t, float)
    y = zeros(t, float)
     
    x[0] = 1.0
    y[0] = x[0]
    i = 1
    while i< t:
        x[i] = x[i-1] + system(x[i-1])*dt + normal(0, dt)
        y[i] = obs(x[i]) + Q*normal(0, 1.0)
        i += 1
   
    fig1 = figure(1)
    plot(x, 'r-')
    grid()
    fig1.show()

    yhat = implicit_filter(y, 100)    
    
    figure(2)
    plot(x, 'r-')
    plot(y, 'b-')
    plot(yhat, 'g-', lw=0.5)
    grid()
    show()

