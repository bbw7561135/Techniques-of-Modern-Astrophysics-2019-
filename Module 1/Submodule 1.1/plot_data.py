import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc

(x, f) = np.loadtxt("data.txt", unpack=True)

# a quick plot
#plt.plot(x, f)
#plt.show()

# now a longer, but better, plot

rc('text', usetex=True)
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size':14})
plt.rc('axes', labelsize=16)
plt.rcParams.update({'figure.autolayout': True})

fig = plt.figure(figsize=(6, 4))
ax = fig.add_subplot(1,1,1)
#ax.grid(True)
ax.plot(x, f, color="black")
ax.set_xlabel('$x$')
ax.set_ylabel('$f(x)$')

plt.savefig("plot_data.pdf")
#plt.show()
