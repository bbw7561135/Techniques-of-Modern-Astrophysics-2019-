import numpy as np
import matplotlib.pyplot as plt

t, ex, ey, ez, sx, sy, sz, E = np.loadtxt('data1.txt', unpack=True)

#plt.plot(t, ex)
#plt.plot(t,sx)
#plt.plot(ex,ey)
#plt.plot(sx,sy)
plt.plot(t,E)
plt.show()
