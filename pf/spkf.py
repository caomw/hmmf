#!/usr/bin/python

from numpy import *
import numpy.matlib as nm
from numpy.random import *
from pylab import *
import mpl_toolkits.mplot3d.axes3d as p3
set_printoptions(precision=3)
set_printoptions(suppress=True)

sigma_process = nm.matrix([[1.0, 0., 0.], [0., 1., 0.], [0., 0., 1.]])*0
sigma_predict = nm.matrix([[1.0, 0., 0.], [0., 1., 0.], [0., 0., 1.]])
sigma_init = nm.matrix([[1.0, 0., 0.], [0., 1., 0.], [0., 0., 1.]])/10
sigma_good = nm.matrix([[1.0, 0., 0.], [0., 1., 0.], [0., 0., 1.]])/100

dt = 0.005
def system(x, w):
    x = x.T    
    x0 = x[0]
    x1 = x[1]
    x2 = x[2]
    new = nm.matrix([[0., 0., 0.]]).T
    
    new[0] = dt*10.0*(x1 - x0) + x0
    new[1] = dt*x0*(28. - x2) - dt*x1 + x1
    new[2] = dt*(x0*x1 - 8/3.0*x2) + x2
    """
    new[0] = cos(x2)*dt + x0
    new[1] = sin(x2)*dt + x1
    new[2] = 1.0*dt + x2
    """
    return new.T + w

def obs(x, w):
    return x+w

def sigma_point_filter(y, ygood, init):
    init = init.T


    t = len(y[:, 0])
    x0 = init
    yhat = nm.zeros( (t, 3))
    xbar = nm.matrix([x0[0,0], x0[1,0], x0[2,0]])
    yhat[0] = xbar
    print "xbar: ", xbar

    na = 3
    alpha = 1e-4
    beta = 2.0
    kappa = 0.0
    lamda = alpha*alpha*(na + kappa) - na
    xsig = nm.zeros( (2*na + 1, na))
    ysig = nm.zeros( (2*na + 1, 3))
    wsigm = nm.zeros( (2*na+1, 1))
    wsigc = nm.zeros( (2*na+1, 1))
    
    Px = nm.zeros( (3, 3))
    Px = sigma_init
    for k in arange(1, t, 1):
        sqP = linalg.cholesky((na + lamda)*Px)
        
        # 1. calculate sigma points
        for i in range(2*na + 1):
            if i == 0:
                wsigm[i] = lamda/(na + lamda)
                wsigc[i] = wsigm[i]
                xsig[i] = xbar
            elif (i <= na):
                xsig[i] = xbar + sqP[i-1,:]
                wsigm[i] = 1.0/2/(na + lamda)
                wsigc[i] = wsigm[i]
            elif (i > na):
                xsig[i] = xbar - sqP[i-1-na,:]
                wsigm[i] = 1.0/2/(na + lamda)
                wsigc[i] = wsigm[i]
            #print xsig[i]

        # 2. time-update equations
        xbar = 0*xbar
        for i in range(2*na + 1):
            xsig[i] = system(xsig[i], 0*xsig[i])
            xbar += xsig[i]*wsigm[i,0]

        Px = nm.zeros( (3, 3))
        for i in range(2*na + 1):
            Px += wsigc[i,0]*(nm.matrix(xsig[i]) - xbar).T*(nm.matrix(xsig[i]) - xbar)
        Px += sigma_process
        #print "P: ", P
        #print "xbar: ", xbar

        ybar = nm.matrix([0., 0., 0.])
        # obs equation
        for i in range(2*na + 1):
            ysig[i] = obs(xsig[i], 0*xsig[i])
            ybar += ysig[i]*wsigm[i, 0]
        
        # 3. measurement equations
        Pyy = nm.zeros( (3, 3))
        Pxy = nm.zeros( (3, 3))
        for i in range(2*na + 1):
            Pyy += wsigc[i,0]*nm.matrix(ysig[i] - ybar).T*(nm.matrix(ysig[i] - ybar))
            Pxy += wsigc[i,0]*nm.matrix(xsig[i] - xbar).T*(nm.matrix(ysig[i] - ybar))

        Pyy += sigma_predict
        Kt = Pxy*inv(Pyy)
        
        if (k%50 == 0):
            nu = nm.matrix(ygood[k] - ybar)
        else:
            nu = nm.matrix(y[k] - ybar)
        xbar = xbar + (Kt*nu.T).T
        Px = Px - Kt*Pyy*Kt.T
        
        yhat[k] = xbar
        #if( k%100 == 0):
        #    print k
    return yhat


if __name__=="__main__":


    N = 1000
    x = nm.zeros((N, 3))
    x[0] = [1.0, 1.0, 0.0]
    time = arange(0, N*dt, dt)

    y = nm.zeros((N, 3))
    ygood = nm.zeros((N, 3))

    for i in range(N - 1):
        x[i+1] = system(x[i], multivariate_normal([0., 0., 0.], sigma_process))
    
    for i in range(N):
        y[i] = obs(x[i], multivariate_normal([0., 0., 0.], sigma_predict))
        ygood[i] = obs(x[i], multivariate_normal([0., 0., 0.], sigma_good))

    print "x: ", x[0] 
    
    for i in range(10):
        yhat = sigma_point_filter(y, ygood, x[0] + multivariate_normal([0., 0., 0.], sigma_init)) 
    
        diffo, diffn = 0, 0
        for i in range(N):
            diffn += dot(yhat[i] - x[i], (yhat[i] - x[i]).T)
            diffo += dot(y[i] - x[i], (y[i] - x[i]).T)
        diffn = sqrt(diffn)
        diffo = sqrt(diffo)
        print "diff was: %.2f"  %diffo
        print "diff new: %.2f"  %diffn

    """
    x, y = x.T.tolist(), y.T.tolist()
    yhat = yhat.T.tolist()
    figure(1)
    plot(x[0], x[1], 'r-', y[0], y[1], 'b+', yhat[0], yhat[1], 'g-')
    grid()
    
    figure(2)
    plot(time, x[1], 'r-', time, y[1], 'b-', time, yhat[1], 'g-')
    grid()
    figure(3)
    plot(time, x[2], 'r-', time, y[2], 'b-', time, yhat[2], 'g-')
    grid()
    
    show()
    """
