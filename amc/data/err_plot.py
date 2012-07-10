#!/usr/bin/python

import sys
from numpy import *
from pylab import *
from matplotlib.ticker import MultipleLocator, FixedLocator, FormatStrFormatter
from scipy import *

rc('font', family='serif')
rc('text', usetex='True')

to_save = False
save_name1 = 'none.pdf'
save_name2 = 'none.pdf'
if len(sys.argv) > 1:
    to_save = bool(sys.argv[1])
    save_name1 = sys.argv[2]
    save_name2 = sys.argv[3]

data = loadtxt('err_conv_batch.dat', delimiter=' ')
n = data[:,0]
m1 = data[:,1]
m2 = data[:,2]

"""
fig = figure(1)
ax = fig.add_subplot(111, aspect='equal')
setp(ax.get_xticklabels(), fontsize=20)
setp(ax.get_yticklabels(), fontsize=20)
semilogx(n, m1, 'bo-')
axis('tight')
grid(which='majorminor')
xlabel(r'Number of samples')
title(r'$1^{st}$ moment')
if to_save:
    savefig(save_name1, bbox_inches='tight')

fig = figure(2)
ax = fig.add_subplot(111, aspect='equal')
setp(ax.get_xticklabels(), fontsize=20)
setp(ax.get_yticklabels(), fontsize=20)
semilogx(data[:,0], data[:,2], 'bo-')
axis('tight')
grid(which='majorminor')
xlabel(r'Number of samples')
title(r'$2^{nd}$ moment')
if to_save:
    savefig(save_name2, bbox_inches='tight')
"""

fig = figure(3)
#ax = fig.add_subplot(111, aspect='equal')
#setp(ax.get_xticklabels(), fontsize=20)
#setp(ax.get_yticklabels(), fontsize=20)
plot(n, m1/pow(log(n)/n, 1))
grid()

show()
