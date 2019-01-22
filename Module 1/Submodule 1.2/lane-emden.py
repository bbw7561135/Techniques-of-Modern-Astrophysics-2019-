import numpy as np
import scipy as sp
import scipy.integrate
import phy4910
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size':14})
plt.rc('axes', labelsize=16)
plt.rcParams.update({'figure.autolayout': True})

stepsize = 0.001
filename = "assignment1b1.pdf"

n = 3.0

def f(x, y, z):
    return z;

def g(x, y, z):
    return -2.0 / x * z - np.power(y, n)

eta1, rho1, tmp = phy4910.ode_euler(0.00001, 10, stepsize, 1.0, 0.0, f, g)

eta2, rho2, tmp = phy4910.ode_rk4(0.00001, 10, stepsize, 1.0, 0.0, f, g)

# fig = plt.figure(figsize=(6, 4))
# ax = fig.add_subplot(1,1,1)
# ax.grid(True)

#ax.plot(eta1, rho1, color="black")
# ax.plot(eta2, rho2, color="red")
#ax.plot(eta1, np.sin(eta1)/eta1, color="blue", linewidth=0.7)
#ax.plot(eta2, rho2-rho1, color="red")

# plt.savefig("assignment1a1.pdf")

# plt.show()

phy4910.plot_print(eta2, rho2, filename, "$\\eta$", "$\\rho(\\eta)$")

varrho = np.power(rho2[rho2 > 0],n) * np.power(eta2[rho2 > 0],2)

m = sp.integrate.trapz(varrho,eta2[rho2 > 0])
print("Eta_s: {}".format(max(eta2[rho2 >0])))
print("Mass: {}".format(m))
