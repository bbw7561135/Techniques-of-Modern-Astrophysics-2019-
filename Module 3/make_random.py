from nbody import *
from random import random

N = 1000
m = 1.0 / N
t = 0.0

parts = []
for i in range(N):
    
    pos = [ 1.0 - 2.0 * random(), 1.0 - 2.0 * random(), 1.0 - 2.0 * random() ]
    while pos[0]**2 + pos[1]**2 + pos[2]**2 > 1.0:
        pos = [ 1.0 - 2.0 * random(), 1.0 - 2.0 * random(), 1.0 - 2.0 * random() ]
    
    vel =  [ (1.0 - 2.0 * random()), (1.0 - 2.0 * random()), (1.0 - 2.0 * random()) ]
    
    p = particle(m, pos, vel)
    parts.append(p)
    
s = System(N, t, parts)

s.write("random.dat")     

