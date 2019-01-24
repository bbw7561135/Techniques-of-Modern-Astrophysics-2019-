import numpy as np
import scipy as sp
<<<<<<< HEAD
import scipy.integrate
=======
>>>>>>> e97fbb2816f534478310c0ce82e35ef80aea28ca
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
N = 25
rho0 = 3.789e6

<<<<<<< HEAD
# k = 3.166e12
k = 4.936e14
G = 6.674e-8
M_sun = 1.989e33
M = M_sun
delta = 0.1 # integration stop point
=======

n = 1.5
>>>>>>> e97fbb2816f534478310c0ce82e35ef80aea28ca

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

<<<<<<< HEAD
# fig = plt.figure(figsize=(6, 4))
# ax = fig.add_subplot(1,1,1)
# ax.grid(True)

#ax.plot(eta1, rho1, color="black")
# ax.plot(eta2, rho2, color="red")
#ax.plot(eta1, np.sin(eta1)/eta1, color="blue", linewidth=0.7)
#ax.plot(eta2, rho2-rho1, color="red")

# plt.savefig("assignment1a1.pdf")

# plt.show()

# phy4910.plot_print(eta2, rho2, filename, "$\\eta$", "$\\rho(\\eta)$")
#
# varrho = np.power(rho2[rho2 > 0],n) * np.power(eta2[rho2 > 0],2)
#
# m = sp.integrate.trapz(varrho,eta2[rho2 > 0])
# print("Eta_s: {}".format(max(eta2[rho2 > 0])))
# print("Mass: {}".format(m))
#
# rho_c = (M / (4 * np.pi * m)) ** 2 * ((5 * k) / (8 * np.pi * G)) ** -3
# print("Rho_c: {}".format(rho_c))
#
# lambdan = (((n+1) * k * rho_c ** ((1-n)/n))/(4*np.pi*G)) ** 0.5
# r = (lambdan*eta2[rho2 > 0]) / 695e8
# phy4910.plot_print(r, rho_c * rho2[rho2 > 0] / 1e6, filename2, "Radius (solar radii)", "Density ($10^6g/cm^3$)")
# print("Radius: {}".format(max(r)))

eta2, rho2, tmp = phy4910.ode_rk4(0.00001, 10, stepsize, 1.0, 0.0, f, rnr_g)

phy4910.plot_print(eta2, rho2, filename3, "$\\eta$", "$\\varrho(\\eta)$")
varrho = np.power(rho2[rho2 > 0],n) * np.power(eta2[rho2 > 0],2)

m = sp.integrate.trapz(varrho,eta2[rho2 > 0])
print("Eta_s: {}".format(max(eta2[rho2 > 0])))
print("Mass: {}".format(m))

rho_c = (M / (4 * np.pi * m)) ** 2 * ((5 * k) / (8 * np.pi * G)) ** -3
print("Rho_c: {}".format(rho_c))

rhoCR = np.logspace(4,12,N)
#eta_array = np.empty(N)
#varrho_array = np.empty(N)
#sigma_array = np.empty(N)
mass = np.empty(N)
eta_array = []
varrho_array = []
sigma_array = []
radii = []
lambda2 = np.sqrt(k/(4*np.pi*G*np.power(rho0,1.0/3.0)))

for i in range(N):
    x,y,z = phy4910.ode_rk4(0.00001, 10, stepsize, rhoCR[i] / rho0, 0.0, f, rnr_g)
    eta_array.append(x)
    varrho_array.append(y)
    sigma_array.append(z)
    varrho = np.power(varrho_array[i],n) * np.power(eta_array[i],2)
    mass[i] = sp.integrate.trapz(varrho,eta_array[i])/2.0
    radii.append(max(y) * lambda2 / 695e12)

proper_mass = (mass[i] * (rhoCR/rho0)**(0.5) * 4 * np.pi * (5 * k /(8*np.pi*G) ) ** (3.0/2.0)) * 1e-16
=======
#ax.plot(eta1, rho1, color="black")
ax.plot(eta2, rho2, color="red")
#ax.plot(eta1, np.sin(eta1)/eta1, color="blue", linewidth=0.7)
#ax.plot(eta2, rho2-rho1, color="red")

plt.savefig("assignment1a1.pdf")

#plt.show()

varrho = np.power(rho2[rho2 > 0],n) * np.power(eta2[rho2 > 0],2)

m = sp.integrate.trapz(varrho,x)
print(m)
>>>>>>> e97fbb2816f534478310c0ce82e35ef80aea28ca

phy4910.plot(rhoCR, proper_mass/2e30, "Central Density", "Solar Masses")
print("Stable mass: {}".format(max(proper_mass/2e30)))

phy4910.plot(rhoCR, radii, "Central Density", "Solar Radii/100")

phy4910.plot(proper_mass/2e30,radii, "Solar Masses", "Solar Radii/100")
