import numpy as np
import numpy.random
import matplotlib
import matplotlib.pyplot as plt
import time
import phy4910
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rc
rc('text', usetex=True)
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size':14})
plt.rc('axes', labelsize=14)
plt.rcParams.update({'figure.autolayout': True})

filename = "assignment2b1.pdf"
filename2 = "assignment2b2.pdf"

# final x,y coordinates lists
xf = []
yf = []
# final cosTheta value list
mu = []
N = 10 # number of bins
n = int(100000) # number of photons


def Move_Photon(zMax=1, tauMax=10):

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

for i in range(n):
    cosTheta, phi, xt, yt = Move_Photon()
    xf.append(xt)
    yf.append(yt)
    mu.append(cosTheta)


bins = np.arange(0,1,1/N)
bins += (1/N)/2.0

bin = np.zeros(N)

for angle in mu:
    pos = int(angle * N)
    bin[pos] += 1

# fig = plt.figure()
# ax = fig.add_subplot(111)#, projection='3d')

bin = bin / (2 * bins * n / N)

bins = np.arccos(bins)

itheory = 1.0/2.0 + 3.0/4.0 * np.cos(bins) # theoretical scatter angles

print(time.time()-starttime,'s')

phy4910.plot_print(bins, [bin, itheory], filename, "Angle [$^\\circ$]", "Normalized Counts", labels=("calculated","theoretical"),noline=True)

fig = plt.figure()
ax = fig.add_subplot(111)
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.scatter(xf, yf, s=3, alpha=0.5)
ax.grid(True)
plt.savefig(filename2)
plt.show()
