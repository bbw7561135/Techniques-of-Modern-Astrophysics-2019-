from nbody import *
import matplotlib.pyplot as plt
import numpy as np

solar = System.read("solar_system.dat")

solar.G = 0.0001184
m_increase = 1.0

t_end = 500.0
dt = 1e-3
dt_out = 0.1
dt_log = 0.01

r = []

sun = solar.parts[0]
mercury = solar.parts[1]
venus = solar.parts[2]
earth = solar.parts[3]
mars = solar.parts[4]
jupiter = solar.parts[5]
saturn = solar.parts[6]
neptune = solar.parts[7]
uranus = solar.parts[8]

hillradii = []

for i in range(1,9):
	solar.parts[i].m *= m_increase
	hillradii.append(sun.hillradius(solar.parts[i]))

tl = []

d_mercury = []
d_venus = []
d_earth = []
d_mars = []
d_jupiter = []
d_saturn = []
d_neptune = []
d_uranus = []

px_doom = 12.0
py_doom = 0
pz_doom = 45.0
vx_doom = 0
vy_doom = 0
vz_doom = -5.0
m_doom = sun.m
d_doom = []

doom = particle(m_doom, [px_doom, py_doom, pz_doom], [vx_doom, vy_doom, vz_doom])

solar.parts.append(doom)
solar.N += 1

solar.accel_calc()

v = np.sqrt(vx_doom**2.0 + vy_doom**2.0 + vz_doom**2.0)
b = np.sqrt(px_doom**2.0 + py_doom**2.0 + pz_doom**2.0)

# with open("temp.txt", 'w') as f:
# 	for part in solar.parts:
# 		f.write(part.__str__())
# 		f.write("\n")

t = 0.0
t_out = 0.0
t_log = 0.0
step = 0

while t < t_end:
	for part in solar.parts:
		part.pos[0] += part.vel[0]*dt + 0.5*part.acc[0]*dt**2.0
		part.pos[1] += part.vel[1]*dt + 0.5*part.acc[1]*dt**2.0
		part.pos[2] += part.vel[2]*dt + 0.5*part.acc[2]*dt**2.0

		part.vel[0] += 0.5*part.acc[0]*dt
		part.vel[1] += 0.5*part.acc[1]*dt
		part.vel[2] += 0.5*part.acc[2]*dt

	solar.accel_calc()

	for part in solar.parts:

		part.vel[0] += 0.5*part.acc[0]*dt
		part.vel[1] += 0.5*part.acc[1]*dt
		part.vel[2] += 0.5*part.acc[2]*dt

	t += dt
	step += 1
	solar.t = t

	if t >= t_out:
		filename = f"solar/run_{step:06d}"
		solar.write(filename)
		t_out += dt_out
		d_mercury.append(sun.distance(mercury))
		d_venus.append(sun.distance(venus))
		d_earth.append(sun.distance(earth))
		d_mars.append(sun.distance(mars))
		d_jupiter.append(sun.distance(jupiter))
		d_saturn.append(sun.distance(saturn))
		d_neptune.append(sun.distance(neptune))
		d_uranus.append(sun.distance(uranus))
		d_doom.append(sun.distance(doom))
		tl.append(t)

	if t >= t_log:
		print(f"step:{step} t = {t}")
		t_log += dt_log

labels = ["Mercury","Venus","Earth","Mars","Jupiter","Saturn","Neptune","Uranus","Nemesis"]

fig = plt.figure(figsize=(12,7))
ax = fig.add_subplot(111)
ax.plot(tl,d_mercury, label=labels[0])
ax.plot(tl,d_venus, label=labels[1])
ax.plot(tl,d_earth, label=labels[2])
ax.plot(tl,d_mars, label=labels[3])
ax.plot(tl,d_jupiter, label=labels[4])
ax.plot(tl,d_saturn, label=labels[5])
ax.plot(tl,d_neptune, label=labels[6])
ax.plot(tl,d_uranus, label=labels[7])
ax.plot(tl,d_doom, label=labels[8])
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.set_ylim(-5,60)
# ax.set_xlim(-10,150)
ax.set_xlabel('Time [yr]')
ax.set_ylabel('Distance')
plt.savefig('doom-distance-{}v-{}b.pdf'.format(v,b))
plt.show()

print(hillradii)
# np.savetxt("hillradii.txt", hillradii)
