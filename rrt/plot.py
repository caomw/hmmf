#!/usr/bin/python

import sys
from numpy import *
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Circle


if len(sys.argv) <= 1:
    print "usage:", sys.argv[0], "<obs file>"
    sys.exit() 

obsfile = sys.argv[1]
obs = open( obsfile, 'r');

points= open('points.dat', 'r');
bot = open('input/bot.txt', 'r');
goal = open('input/goal.txt', 'r');
box = open('input/box.txt', 'r');

num_obs = 0
obsx = []
obsy = []
obsr = []

fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')

if obs:
    l = obs.readline()
    while l:
        s = l.split(',')
        obsx.append( float(s[0]) )
        obsy.append( float(s[1]) )
        obsr.append( float(s[2]) )
        num_obs += 1
        l = obs.readline()
    obs.close()

if bot:
    l = bot.readline()
    if l:
        s = l.split(',')
        botx = float(s[0])
        boty = float(s[1])
        botr = float(s[2])
    bot.close()

if goal:
    l = goal.readline()
    if l:
        s = l.split(',')
        goalx = float(s[0])
        goaly = float(s[1])
        goalr = float(s[2])
    goal.close()

if box:
    l = box.readline()
    if l:
        s = l.split(',')
        box_cx = float(s[0])
        box_cy = float(s[1])
        box_xsize = float(s[2])
        box_ysize = float(s[3])
    box.close()

if points:
    
    found_optpath = 0
    # Read cost on the first line
    l = points.readline()
    s = l.split('\n')
    if len(s) == 2:
        print "Path Cost: ", float(s[0])
    
    while l:
        
        l = points.readline()
        s = l.split(' ')
        #print s

        if (len(s) == 3) and (s[0] == 'Duration:'):
            print "Duration: ", float(s[1]), '[ms]'
        
        else:
            if (len(s) == 3) and (s[0] == 'optpath_cost:'):
                print "optpath_cost: ", float(s[1])
            else: 
                if len(s) == 3:
                    px = []
                    py = []
                    px.append( float( s[0]))
                    py.append( float( s[1]))

                    l = points.readline()
                    s = l.split(' ')
                    px.append( float( s[0]))
                    py.append( float( s[1]))
        
                    if found_optpath == 1:
                        plt.plot(px, py, 'ko-', linewidth=1.5)
                    else:
                        plt.plot(px, py, 'y-', linewidth=1.5)
        
            if (len(s) == 1) and (s[0] == 'optpath\n'):
                print "Found optpath"
                found_optpath = 1

       
box_minx = box_cx - box_xsize/2
box_maxx = box_cx + box_xsize/2

box_miny = box_cy - box_ysize/2
box_maxy = box_cy + box_ysize/2


for i in range(num_obs):
    circle = Circle( (obsx[i], obsy[i]), obsr[i], fc='blue', alpha=.2)
    ax.add_patch(circle)

circle = Circle( (botx, boty), botr, fc='red', alpha=.5)
ax.add_patch(circle)

circle = Circle( (goalx, goaly), goalr, fc='green', alpha=.5)
ax.add_patch(circle)

ax.set_xlim(box_minx, box_maxx)
ax.set_ylim(box_miny, box_maxy)

plt.grid()
plt.show()

