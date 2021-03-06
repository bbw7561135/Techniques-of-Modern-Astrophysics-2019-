from nbody import *
from random import *
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d
import matplotlib.animation as animation

from sys import argv
from time import sleep

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.axis('off')
ax.set_aspect('equal')

s = System.read(argv[1])
x = s.all_x()
y = s.all_y()
z = s.all_z()

line, = ax.plot(x, y, z,  color="black", marker="o", markersize=1, alpha=1, linestyle='none')
#ax.set_xlim(-1, 1)
#ax.set_ylim(-1, 1)


def update(fname):
    #sleep(1)
    s = System.read(fname)
    x = s.all_x()
    y = s.all_y()
    z = s.all_z()
    
    #print "animating", fname
    
    line.set_data(x, y)
    line.set_3d_properties(zs=z)
    

if len(argv) > 2:
    ani = animation.FuncAnimation(fig, update, frames=argv[2:], interval=20)

plt.show()

