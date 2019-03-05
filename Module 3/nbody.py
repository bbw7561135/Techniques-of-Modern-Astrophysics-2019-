from numpy import sqrt
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

	def write(self, filename):
		with open(filename, "w") as f:
			f.write(str(self))
		return
	@classmethod
	def read(cls, filename):
		with open(filename, "r") as f:
 			N = int(f.readline())
			t = float(f.readline())
			G = float(f.readline())
			parts = []
			for i in range(N):


		return
if __name__=="__main__":
	comet = particle(1.0, [1,1,1], [5489,83,92])
	baseball = particle(1e-8, [-7,-4,8], [3.14,7,38])
	sun = particle(0.9, [-0.1,0,0], [0,0.133,0])
	earth = particle(0.1, [0.9,0,0], [0,-1.2,0])
	solar = System(4,0,[comet,baseball, sun, earth])
	solar.write("solar.txt")
