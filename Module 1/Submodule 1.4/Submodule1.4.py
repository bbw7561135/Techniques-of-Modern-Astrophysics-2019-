import numpy as np
import numpy.random
import matplotlib.pyplot as plt
import time
from mpl_toolkits.mplot3d import Axes3D

# final x,y coordinates lists
xf = []
yf = []

def Move_Photon(zMax =1, tauMax = 10):

    N = 10**3
    x = np.zeros(N)
    y = np.zeros(N)
    z = np.zeros(N)

    scatters = 0

    flag = True

    while flag:

        scatters = 0
        #print('Start')
        for i in range(N-1):
            chi,xi,nu = np.random.random(3)

            cosTheta = 1-2*chi
            sinTheta = np.sqrt(1-cosTheta**2)

            phi = 2*np.pi*xi


            tau = -np.log(1 - nu)

            s = tau/tauMax

            x[i+1] = x[i]+ s*sinTheta*np.cos(phi)
            y[i+1] = y[i]+ s*sinTheta*np.sin(phi)
            z[i+1] = z[i]+ s*cosTheta

            if z[i+1] < 0:
                scatters = 0
                x = np.zeros(N)
                y = np.zeros(N)
                z = np.zeros(N)
                #print('Reset')
                break

            if z[i+1] > zMax:
                flag = False
                #print('Escaped')
                break
            scatters +=1

    return cosTheta, phi, x[scatters],y[scatters]

starttime = time.time()

for i in range(1000):
    cosTheta, phi, xt, yt = Move_Photon()
    xf.append(xt)
    yf.append(yt)

print(time.time()-starttime,'s')

# print(scatters)

fig = plt.figure()
ax = fig.add_subplot(111)#, projection='3d')

plt.scatter(xf,yf, s=4.0, alpha=0.5)
plt.show()
