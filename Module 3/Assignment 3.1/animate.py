from nbody import *
from random import *
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d
import matplotlib.animation as animation

from sys import argv
from time import sleep

fig = plt.figure()
ax = fig.add_subplot(111) #, projection='3d')
ax.axis('off')
ax.set_aspect('equal')

s = System.read(argv[1])
# x = s.all_x()
# y = s.all_y()
# z = s.all_z()
halox = []
haloy = []
haloz = []
diskx = []
disky = []
diskz = []

def separateParts():
    thalox = []
    thaloy = []
    thaloz = []
    tdiskx = []
    tdisky = []
    tdiskz = []
    for p in s.parts:
            if p.ID == 1:
                thalox.append(p.pos[0])
                thaloy.append(p.pos[1])
                thaloz.append(p.pos[2])
            else:
                tdiskx.append(p.pos[0])
                tdisky.append(p.pos[1])
                tdiskz.append(p.pos[2])
    return thalox, thaloy, thaloz, tdiskx, tdisky, tdiskz

halox, haloy, haloz, diskx, disky, diskz = separateParts()

# line, = ax.plot(halox, haloy, haloz,  color="#FF8000", marker="o", markersize=1, alpha=1, linestyle='none')
# line2, = ax.plot(diskx, disky, diskz,  color="#FFFF00", marker="o", markersize=1, alpha=1, linestyle='none')
line, = ax.plot(halox, haloy,  color="#FF8000", marker="o", markersize=1, alpha=1, linestyle='none')
line2, = ax.plot(diskx, disky,  color="#FFFF00", marker="o", markersize=1, alpha=1, linestyle='none')
#ax.set_xlim(-1, 1)
#ax.set_ylim(-1, 1)


def update(fname):
    #sleep(1)
    s = System.read(fname)
    # x = s.all_x()
    # y = s.all_y()
    # z = s.all_z()

    halox, haloy, haloz, diskx, disky, diskz = separateParts()

    #print "animating", fname

    line.set_data(halox, haloy)
    line2.set_data(diskx, disky)
    # line.set_data(x, y)
    # line.set_3d_properties(zs=z)


if len(argv) > 2:
    writer = animation.FFMpegWriter()
    ani = animation.FuncAnimation(fig, update, frames=argv[2:], interval=20)
    ani.save("antennae-collision.mpg", fps=30, extra_args=['-vcodec','libx264'])

plt.show()
