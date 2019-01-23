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
filename2 = "assignment1b1-2.pdf"
filename3 = "assignment1c3.pdf"

n = 3.0

# k = 3.166e12
k = 4.936e14
G = 6.674e-8
M_sun = 1.989e33
M = M_sun
delta = 0.1 # integration stop point

def f(x, y, z):
    return z

def g(x, y, z):
    return -2.0 / x * z - np.power(y, n)

def A(x):
    return -(5.0/9.0) * np.power(x,-4.0/3.0) * (1 + np.power(x,2.0/3.0)) ** (-0.5) - (2.0/3.0) * np.power(x,-2.0/3.0) * (1 + np.power(x,2.0/3.0)) ** (-1.5) + (1.0/3.0) * (1+np.power(x,2.0/3.0))**(-5.0/2.0)

def B(x):
    return 5.0/3.0 * np.power(x,-1.0/3.0) * (1+np.power(x,2.0/3.0))**(-0.5) - 1.0/3.0 * np.power(x,1.0/3.0) * (1+np.power(x,2.0/3.0))**(-3.0/2.0)

def rnr_g(x,y,z):
    return -2.0/x * z - A(y)/B(y) * z**2 - 1/B(y) * y

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
print("Eta_s: {}".format(max(eta2[rho2 > 0])))
print("Mass: {}".format(m))

rho_c = (M / (4 * np.pi * m)) ** 2 * ((5 * k) / (8 * np.pi * G)) ** -3
print("Rho_c: {}".format(rho_c))

lambdan = (((n+1) * k * rho_c ** ((1-n)/n))/(4*np.pi*G)) ** 0.5
r = (lambdan*eta2[rho2 > 0]) / 695e8
phy4910.plot_print(r, rho_c * rho2[rho2 > 0] / 1e6, filename2, "Radius (solar radii)", "Density ($10^6g/cm^3$)")
print("Radius: {}".format(max(r)))

eta2, rho2, tmp = phy4910.ode_rk4(0.00001, 10, stepsize, 1.0, 0.0, f, rnr_g)

phy4910.plot_print(eta2, rho2, filename3, "$\\eta$", "$\\varrho(\\eta)$")
