#!/usr/bin/python

from numpy import *
from pylab import *
from scipy import *


data = loadtxt('err.dat', delimiter=' ')

figure(1)
plot(data[:,0], data[:,1], 'b-')
pcoeff = polyfit(data[:,0], data[:,1], 10)
yfit = polyval(pcoeff, data[:,0])
plot(data[:,0], yfit, 'r-', lw=1.5, label='poly. fit')
pcoeff = polyfit(data[:,0], data[:,1], 1)
yfit = polyval(pcoeff, data[:,0])
plot(data[:,0], yfit, 'g-', lw=1.5, label='linear fit')
grid()
xlabel('No. of samples')
ylabel('Err. in first moment')
legend()
savefig('err_convergence_random_m1.pdf', bbox_inches='tight')

figure(2)
plot(data[:,0], data[:,2], 'b-')
pcoeff = polyfit(data[:,0], data[:,2], 10)
yfit = polyval(pcoeff, data[:,0])
plot(data[:,0], yfit, 'r-', lw=1.5, label='poly. fit')
pcoeff = polyfit(data[:,0], data[:,2], 1)
yfit = polyval(pcoeff, data[:,0])
plot(data[:,0], yfit, 'g-', lw=1.5, label='linear fit')
grid()
xlabel('No. of samples')
ylabel('Err. in second moment')
legend()
savefig('err_convergence_random_m2.pdf', bbox_inches='tight')

show()
