from nbody import *

solar = System.read("random.txt")

t_end = 5.0
dt = 1e-3
dt_out = 0.01
dt_log = 0.01


solar.accel_calc()

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
		filename = f"run_{step:06d}"
		solar.write(filename)
		t_out += dt_out

	if t >= t_log:
		print(f"step:{step} t = {t}")
		t_log += dt_log
