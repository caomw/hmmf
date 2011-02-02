#!/usr/bin/python

from numpy import *
from numpy.random import *
from pylab import *
import mpl_toolkits.mplot3d.axes3d as p3

sigma_process = array([[1.0, 0., 0.], [0., 1., 0.], [0., 0., 1.]])/50
sigma_predict = array([[1.0, 0., 0.], [0., 1., 0.], [0., 0., 1.]])/50
sigma_init = array([[1.0, 0., 0.], [0., 1., 0.], [0., 0., 1.]])/10

dt = 0.01
def system(x):
    new = [0., 0., 0.]
    new[0] = dt*8.0*(x[1] - x[0]) + x[0]
    new[1] = dt*x[0]*(25. - x[2]) - dt*x[1] + x[1]
    new[2] = dt*(x[0]*x[1] - 8/3.0*x[2]) + x[2]
    
    return new

def normal_prob(predict, reading, sigma):
    delta = reading - predict
    temp = 1 /( ((2*pi)**3/2)*sqrt(det(sigma)) )* exp(-0.5*dot(delta, dot(inv(sigma), delta)))
    return temp


def resample(x, w):
    N = len(x[:,0])
    new_x = zeros( (N, 3), float)
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
    t = len(y[:, 0])
    x0 = init
    xpart = zeros( (N, 3), float)
    w = zeros(N, float)
    xhat = zeros( (t, 3), float)
    
    # init N particles
    for i in range(N):
        xpart[i,:] = x0 + multivariate_normal([0., 0., 0.], sigma_init)
        w[i] = 1.0/N
    
    for k in range(t):
        for i in range(N):
            # predict
            xpart[i] = system(xpart[i])
            # update weights
            w[i] = w[i]*normal_prob(xpart[i], y[k], sigma_predict)
        
        ct = std(w)**2/mean(w)**2
        w /= sum(w)
        if 1/(1 + ct) < 0.95:
            #resample
            xpart, w = resample(xpart, w)
        
        for i in range(N):
            xhat[k] += xpart[i]*w[i]

    return xhat


if __name__=="__main__":

    t = 500
    seed(0)
    x = zeros((t, 3), float)
    x[0] = array([1.0, 1.0, 1.0])

    y = zeros((t, 3), float)
    y[0] = x[0]

    for i in range(t - 1):
        x[i+1] = system(x[i])
        tmp = x[i+1]
        y[i+1] = tmp
        #+ multivariate_normal([0., 0., 0.], sigma_process)
    
    xhat = particle_filter(y, 150, x[0])
   
    
    fig = figure()
    ax = p3.Axes3D(fig)
    x, y, xhat = x.T, y.T, xhat.T
    ax.plot3D(x[0], x[1], x[2], 'r-')
    ax.plot3D(y[0], y[1], y[2], 'b-')
    ax.plot3D(xhat[0], xhat[1], xhat[2], 'g-')

    show()

