#!/usr/bin/env python3

import numpy as np
import random as random
from random import randint
import phy4910
import matplotlib
import matplotlib.pyplot as plt
import scipy as sp
from scipy.integrate import simps


def f(x, y, z):
    return z;
    
def g(x, y, z):
    a = (-5.0/9.0) * y**(-4.0/3.0) * (1.0 + y**(2.0/3.0))**(-1.0/2.0) - (2.0/3.0) * y**(-2.0/3.0) * (1.0 + y**(2.0/3.0))**(-3.0/2.0) + (1.0/3.0) * (1.0 + y**(2.0/3.0))**(-5.0/2.0);
    b =  (5.0/3.0) * y**(-1.0/3.0) * (1.0 + y**(2.0/3.0))**(-1.0/2.0) - (1.0/3.0) * y**(1.0/3.0) * (1.0 + y**(2.0/3.0))**(-3.0/2.0);
    #print(x, y, z) 
    return (-2.0 / x * z) - (a / b) * z**2.0 - (1.0 / b) * y;

N = 25
rho_c = np.logspace(4, 12, num=N)
dm = np.zeros(N)
rho_0 = 3.789E6
d_rho_c = rho_c/rho_0 #This is to get dimensionless central density
mass = np.zeros(N)
radius = np.zeros(N)
k_nr = 3.166E12
k_r = 4.936E14
G = 6.674E-8
smconv = 1.989E33 # g
lam = ((k_nr * rho_0**(-1.0 / 3.0)) / (4.0 * sp.pi * G))**(1.0/2.0)
srconv = 6.957E10 #cm
	


for i in range(N):
	eta2, varrho2, tmp = phy4910.ode_rk4(0.00001, 20, 0.0001, d_rho_c[i], 0.0, f, g)
	varrho_min = varrho2[varrho2> d_rho_c[i] * 0.001]
	eta2_min= eta2[varrho2> d_rho_c[i] * 0.001]
	dm[i] = simps(varrho_min * eta2_min**2.0, eta2_min)
	mass[i] = (sp.pi * 4.0 * rho_0 * dm[i] * lam**3.0) / smconv
	radius[i] = (lam * eta2_min[-1]) / srconv


"""print(f"Dimensionless mass equals {dm}")
print (lam)
print (mass)
print (radius)"""


fig = plt.figure(figsize=(6, 4))
ax = fig.add_subplot(1,1,1)
ax.grid(True)


ax.semilogx(d_rho_c, mass, '.', color="red")
plt.ticklabel_format(style = 'sci', axis = 'y')
plt.ylabel('mass')
plt.xlabel(r'log($\rho_c$)')
plt.savefig('bakhtiar_assignment1_C_4_mass-den.pdf')

fig = plt.figure(figsize=(6, 4))
ax = fig.add_subplot(1,1,1)
ax.grid(True)

ax.semilogx(d_rho_c, radius, '.', color="blue")
plt.ticklabel_format(style = 'sci', axis = 'y')
plt.xlabel(r'log($\rho_c$)')
plt.ylabel('radius')
plt.savefig('bakhtiar_assignment1_C_4_radius-den.pdf')

fig = plt.figure(figsize=(6, 4))
ax = fig.add_subplot(1,1,1)
ax.grid(True)

ax.plot(mass, radius, '.', color="green")
plt.ticklabel_format(style = 'sci', axis = 'y')
plt.xlabel('mass')
plt.ylabel('radius')
plt.savefig('bakhtiar_assignment1_C_4_radius-mass.pdf')
plt.show()
