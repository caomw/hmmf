#!/usr/bin/python

from numpy import *
from numpy.random import *
from pylab import *
import mpl_toolkits.mplot3d.axes3d as p3

sigma_process = array([[1.0, 0., 0.], [0., 1., 0.], [0., 0., 1.]])/100
sigma_predict = array([[1.0, 0., 0.], [0., 1., 0.], [0., 0., 1.]])/50
sigma_init = array([[1.0, 0., 0.], [0., 1., 0.], [0., 0., 1.]])/10

dt = 0.01
def system(x, w):
    
    x0 = x[0]
    x1 = x[1]
    x2 = x[2]
    new = array([0., 0., 0.])
    """
    new[0] = dt*8.0*(x1 - x0) + x0
    new[1] = dt*x0*(25. - x2) - dt*x1 + x1
    new[2] = dt*(x0*x1 - 8/3.0*x2) + x2
    """
    new[0] = cos(x2)*dt + x0
    new[1] = sin(x2)*dt + x1
    new[2] = 1.0*dt + x2
    return new + w

def obs(x, w):
    return x+w

def sigma_point_filter(y, init):
    
    na = 9
    alpha = 0.5
    k = 0.0
    beta = 2.0
    #lamda = alpha*alpha*(na + k) - na
    lamda = 3 - na

    t = len(y[:, 0])
    x0 = init
    P0 = sigma_init
    Pa0 = zeros( (9, 9), float)
    Pa0[0:3, 0:3] = sigma_init
    Pa0[3:6, 3:6] = sigma_process
    Pa0[6:9, 6:9] = sigma_predict
    P = Pa0
    yhat = zeros( (t, 3), float)
    xbar = array([[x0[0], x0[1], x0[2], 0, 0, 0, 0, 0, 0]])

    xsig = zeros( (2*na + 1, na), float)
    ysig = zeros( (2*na + 1, 3), float)
    wsigm = zeros( (2*na+1, 1), float)
    wsigc = zeros( (2*na+1, 1), float)
    
    for k in arange(0, t, 1):
        sqP = linalg.cholesky((na + lamda)*P)
        
        # 1. calculate sigma points
        for i in range(2*na + 1):
            if i == 0:
                wsigm[i] = lamda/(na + lamda)
                wsigc[i] = wsigm[i] + (1 -alpha*alpha + beta)
                xsig[i] = xbar
            elif (i <= na):
                xsig[i] = xbar + sqP[:,i-1]
                wsigm[i] = 1.0/2/(na + lamda)
                wsigc[i] = wsigm[i]
            elif (i > na):
                xsig[i] = xbar - sqP[:, i-1-na]
                wsigm[i] = 1.0/2/(na + lamda)
                wsigc[i] = wsigm[i]
            #print xsig[i]

        # 2. time-update equations
        xbar = 0*xbar
        for i in range(2*na + 1):
            xsig[i,0:3] = system(xsig[i, 0:3], 0*xsig[i, 3:6])
            xbar += xsig[i]*wsigm[i]

        P = 0*P
        for i in range(2*na + 1):
            P += wsigc[i]*dot( array(xsig[i] - xbar).T, (xsig[i] - xbar))
        #print "P: ", P
        #print "xbar: ", xbar
        
        ybar = array([[0., 0., 0.]])
        # obs equation
        for i in range(2*na + 1):
            ysig[i] = obs(xsig[i, 0:3], 0*xsig[i, 6:9])
            ybar += ysig[i]*wsigm[i]
        
        # 3. measurement equations
        Pyy = zeros( (3, 3), float)
        Pxy = zeros( (9, 3), float)
        for i in range(2*na + 1):
            Pyy += wsigc[i]*dot( array(ysig[i] - ybar).T, (ysig[i] - ybar))
            Pxy += wsigc[i]*dot( array(xsig[i] - xbar).T, (ysig[i] - ybar))
         
        Kt = dot(Pxy, inv(Pyy))
       
        nu = y[k] - ybar
        xbar = xbar + dot(Kt, nu[0])
        P = P - dot(Kt, dot(Pyy, Kt.T))
        #print "P: ", P

        yhat[k] = xbar[0, 0:3]
        #print y[k], yhat[k]
        #raw_input()
        print k

    return yhat



if __name__=="__main__":

    t = 500
    seed(0)
    x = zeros((t, 3), float)
    x[0] = array([[1.0, 1.0, 1.0]])

    y = zeros((t, 3), float)
    y[0] = x[0]

    for i in range(t - 1):
        x[i+1] = system(x[i], multivariate_normal([0., 0., 0.], sigma_process))
    
    for i in range(t):
        y[i] = obs(x[i], multivariate_normal([0., 0., 0.], sigma_predict))
   
    
    yhat = sigma_point_filter(y, x[0] + multivariate_normal([0., 0., 0.], sigma_init)) 
    
    #ax = p3.Axes3D(fig)
    x, y = x.T, y.T
    yhat = yhat.T
    #ax.plot3D(x[0], x[1], x[2], 'r-')
    #ax.plot3D(y[0], y[1], y[2], 'b-')
    #ax.plot3D(yhat[0], yhat[1], yhat[2], 'g-')
    figure(1)
    plot(x[0], 'r-', y[0], 'b-', yhat[0], 'g-')
    grid()
    figure(2)
    plot(x[1], 'r-', y[1], 'b-', yhat[1], 'g-')
    grid()
    figure(3)
    plot(x[2], 'r-', y[2], 'b-', yhat[2], 'g-')
    grid()

    show()

