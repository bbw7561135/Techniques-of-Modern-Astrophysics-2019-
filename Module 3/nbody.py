from math import sqrt
from struct import unpack, pack

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
			
			for p in self.parts:
				m = pack('<f', p.m)
				x = pack('<f', p.pos[0])
				y = pack('<f', p.pos[1])
				z = pack('<f', p.pos[2])
				vx = pack('<f', p.vel[0])
				vy = pack('<f', p.vel[1])
				vz = pack('<f', p.vel[2])
				ID = pack('<i', p.ID)
				f.write(m)
				f.write(x)
				f.write(y)
				f.write(z)
				f.write(vx)
				f.write(vy)
				f.write(vz)
				f.write(ID)
				
				

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
			
			parts = []
			for i in range(N):
				m = unpack('<f', f.read(4))[0]
				x = unpack('<f', f.read(4))[0]
				y = unpack('<f', f.read(4))[0]
				z = unpack('<f', f.read(4))[0]
				vx = unpack('<f', f.read(4))[0]
				vy = unpack('<f', f.read(4))[0]
				vz = unpack('<f', f.read(4))[0]
				ID = unpack('<i', f.read(4))[0]
				
				parts.append(particle(m, [x, y, z], [vx, vy, vz], ID))
				
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

	test = System.read("random.dat")
	print(test.parts[0].pos[2])
	test.write("random2.dat")
