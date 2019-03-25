from numpy import sqrt
from nbody import particle

comet = particle(1.0, [1,1,1], [5489,83,92])
baseball = particle(1e-8, [-7,-4,8], [3.14,7,38])
sun = particle(0.9, [-0.1,0,0], [0,0.133,0])
earth = particle(0.1, [0.9,0,0], [0,-1.2,0])

#print(comet.distance(baseball))
#print(baseball.distance(comet))

G = 1
dt = 0.01
t = 0
t_end = 500
count = 0

while t < t_end:
	dist3 = earth.distance(sun)**3.0
	for i in range(3):
		earth.acc[i] = -sun.m * (earth.pos[i] - sun.pos[i]) / dist3
		sun.acc[i] = -earth.m * (sun.pos[i] - earth.pos[i]) / dist3
		earth.vel[i] += earth.acc[i] * dt
		sun.vel[i] += sun.acc[i] * dt
		earth.pos[i] += earth.vel[i] * dt
		sun.pos[i] += sun.vel[i] * dt
	ke_earth = 0.5 * earth.m * earth.v2()
	ke_sun = 0.5 * sun.m * sun.v2()
	pe = - G * earth.m * sun.m  / earth.distance(sun)
	E = ke_earth + ke_sun + pe

	t += dt

	count += 1

	if count % 100 == 0:
		print(t, earth.pos[0], earth.pos[1], earth.pos[2], sun.pos[0], sun.pos[1], sun.pos[2], E)
