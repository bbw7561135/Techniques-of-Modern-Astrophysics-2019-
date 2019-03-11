
#!/usr/bin/env python3
from math import sqrt
from struct import unpack, pack
from array import array

class particle:

	def __init__(self, m, pos, vel, ID=0):
		self.m = m
		self.pos = pos
		self.vel = vel
		self.ID = ID
		self.acc = [0,0,0]
	
	def __str__(self):
		return f"{self.m:6e} {self.pos[0]:6e} {self.pos[1]:6e} {self.pos[2]:6e} {self.vel[0]:6e} {self.vel[1]:6e} {self.vel[2]:6e}"
	
	@classmethod
	def fromstring(cls, part):
		words = part.split()
		m = float(words[0])
		x = float(words[1])
		y = float(words[2])
		z = float(words[3])
		vx = float(words[4])
		vy = float(words[5])
		vz = float(words[6])
		return cls(m,[x,y,z], [vx,vy,vz])
	
	def v2(self):
		return self.vel[0]**2 + self.vel[1]**2 + self.vel[2]**2

	def radius(self):
		return sqrt(self.pos[0]**2.0 + self.pos[1]**2.0 + self.pos[2]**2.0)

	def distance(self, part):
		return sqrt((self.pos[0]-part.pos[0])**2.0 + (self.pos[1]-part.pos[1])**2.0 + (self.pos[2]-part.pos[2])**2.0)



class System:

	def __init__(self, N, t, parts, G=1):
		self.N = N
		self.t = t
		self.parts = parts
		self.G = G
		self.KE = 0
		self.PE = 0
		return
	
	def __str__(self):
		ret = f"{self.N}\n{self.t}\n{self.G}\n"
		for item in self.parts:
			ret += f"{str(item)}\n"
		return ret

	def write_ascii(self, filename):
		with open(filename, "w") as f:
			f.write(str(self))
		return

	def write(self, filename):
		with open(filename, "wb") as f:
			Nb = pack('<i', self.N)
			tb = pack('<d', self.t)
			f.write(Nb)
			f.write(tb)
			
			arr = array('f')
			for i in range(self.N):
				arr.fromlist([self.parts[i].m, self.parts[i].pos[0], self.parts[i].pos[1], self.parts[i].pos[2], self.parts[i].vel[0], self.parts[i].vel[1], self.parts[i].vel[2]])
			arr.tofile(f)

	@classmethod
	def read_ascii(cls, filename):
		with open(filename, "r") as f:
			N = int(f.readline())
			t = float(f.readline())
			G = float(f.readline())
			parts = []
			for i in range(N):
				parts.append(particle.fromstring(f.readline()))
		return cls(N, t, parts, G);
		
	@classmethod
	def read(cls, filename):
		with open(filename, "rb") as f:
			N = unpack('<i', f.read(4))[0]
			t = unpack('<d', f.read(8))[0]
			
			arr = array('f')
			arr.fromfile(f, N * 7)
			
			parts = []
			for i in range(N):
				m = arr[7 * i + 0]
				x = arr[7 * i + 1]
				y = arr[7 * i + 2]
				z = arr[7 * i + 3]
				vx = arr[7 * i + 4]
				vy = arr[7 * i + 5]
				vz = arr[7 * i + 6]
				
				parts.append(particle(m, [x, y, z], [vx, vy, vz]))
		return cls(N, t, parts)
		
				

	def accel_calc(self):
		for i in range(self.N):
			self.parts[i].acc = [0,0,0]
		for i in range(self.N):
			for j in range(i+1, self.N):
				pi = self.parts[i]
				pj = self.parts[j]
				xd = pi.pos[0]-pj.pos[0]
				yd = pi.pos[1]-pj.pos[1]
				zd = pi.pos[2]-pj.pos[2]
				r3 = pi.distance(pj)**3.0
				pi.acc[0] += -self.G*pj.m*xd/r3
				pj.acc[0] += self.G*pi.m*xd/r3

				pi.acc[1] += -self.G*pj.m*yd/r3
				pj.acc[1] += self.G*pi.m*yd/r3

				pi.acc[2] += -self.G*pj.m*zd/r3
				pj.acc[2] += self.G*pi.m*zd/r3
		return

	def all_x(self):
		x = []
		for i in range(self.N):
			x.append(self.parts[i].pos[0])
		return x

	def all_y(self):
		x = []
		for i in range(self.N):
			x.append(self.parts[i].pos[1])
		return x

	def all_z(self):
		x = []
		for i in range(self.N):
			x.append(self.parts[i].pos[2])
		return x






if __name__=="__main__":
	comet = particle(1.0, [1,1,1], [5489,83,92])
	baseball = particle(1e-8, [-7,-4,8], [3.14,7,38])
	sun = particle(0.9, [-0.1,0,0], [0,0.133,0])
	earth = particle(0.1, [0.9,0,0], [0,-1.2,0])
	solar = System(4,0,[comet,baseball, sun, earth])
	solar.write("solar.txt")
	solar.accel_calc()
	solar2 = System.read("solar.txt")
	solar2.write("solar2.txt")
