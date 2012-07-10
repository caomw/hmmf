#!/usr/bin/python

import sys
from numpy import *
from pylab import *
from scipy import *
import matplotlib.ticker as ticker

rc('font', family='serif')
rc('text', usetex='True')

data = loadtxt('errd2.dat')

to_save = False
save_name = 'none.pdf'
if len(sys.argv) > 1:
    to_save = bool(sys.argv[1])
    save_name = sys.argv[2]

fig = figure(1)
ax = fig.add_subplot(111, aspect='equal')
#ax.yaxis.set_major_locator(ticker.LogLocator(base=0.0001))
setp(ax.get_xticklabels(), fontsize=20)
setp(ax.get_yticklabels(), fontsize=20)
semilogx(data[:,0], data[:,2], 'bo-')
pcoeff = polyfit(data[:,0], data[:,2], 10)
yfit = polyval(pcoeff, arange(data[0,0], data[-1,0]))
semilogx(arange(data[0,0], data[-1,0]), yfit, 'r-', lw=1.5, label='poly. fit')
axis('tight')
grid(which='majorminor')
xlabel(r'$n$')
ylabel(r'$\int ||\theta(t) - x(t)||_2^2\ dt$')
legend()

if to_save:
    savefig(save_name, bbox_inches='tight')

show()
