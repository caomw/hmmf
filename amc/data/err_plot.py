#!/usr/bin/python

import sys
from numpy import *
from pylab import *
from scipy import *

to_save = False
save_name1 = 'none.pdf'
save_name2 = 'none.pdf'
if len(sys.argv) > 1:
    to_save = bool(sys.argv[1])
    save_name1 = sys.argv[2]
    save_name2 = sys.argv[3]

data = loadtxt('err.dat', delimiter=' ')

figure(1)
semilogx(data[:,0], data[:,1], 'bo-')
semilogx(data[:,0], pow(data[0,0],1)*data[0,1]/pow(data[:,0],1), 'ro-', label='y = c/n')
semilogx(data[:,0], pow(data[0,0],2)*data[0,1]/pow(data[:,0],2), 'ro-', label='y = c/n')
"""
pcoeff = polyfit(data[:,0], data[:,1], 10)
yfit = polyval(pcoeff, data[:,0])
semilogx(data[:,0], yfit, 'r-', lw=1.5, label='poly. fit')
pcoeff = polyfit(data[:,0], data[:,1], 1)
yfit = polyval(pcoeff, data[:,0])
loglog(data[:,0], yfit, 'g-', lw=1.5, label='linear fit')
"""
axis('tight')
grid(which='minor')
#xticks(linspace(data[0,0], data[-1,0], 20))
xlabel('No. of samples')
ylabel('Err. in first moment')
legend()
if to_save:
    savefig(save_name1, bbox_inches='tight')

figure(2)
loglog(data[:,0], data[:,2], 'bo-')
#pcoeff = polyfit(data[:,0], data[:,2], 10)
#yfit = polyval(pcoeff, data[:,0])
#plot(data[:,0], yfit, 'r-', lw=1.5, label='poly. fit')
#pcoeff = polyfit(data[:,0], data[:,2], 1)
#yfit = polyval(pcoeff, data[:,0])
#plot(data[:,0], yfit, 'g-', lw=1.5, label='linear fit')

axis('tight')
grid(which='minor')
xlabel('No. of samples')
ylabel('Err. in second moment')
legend()
if to_save:
    savefig(save_name2, bbox_inches='tight')

show()
