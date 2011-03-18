#!/usr/bin/python

from numpy import *
from numpy.random import *
from pylab import *
set_printoptions(precision=3)
set_printoptions(suppress=True)

XMAX = 2.0
XMIN = -2.0

class node:
    public:
    private:

def system(x1):
    return 0.8*x1

def obs(x1):
    return x1

if __name__=="__main__":
    
    dt = 0.01
    N = 100
    x = zeros(N, float)
    x[0] = 1.0;
    time = arange(0, N*dt, dt)

    y = zeros(N, float)
    yc = zeros(N, float)

    for i in arange(1, N):
        x[i] = system(x[i-1])
    
    figure()
    plot(time, x, 'r--')
    grid()
    show()
