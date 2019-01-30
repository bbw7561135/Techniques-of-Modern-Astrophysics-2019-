#!/usr/bin/env python3

import numpy as np
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
                                 #x_start, x_end, h, y0, z0,  f, g   
eta2, varrho2, tmp = phy4910.ode_rk4(0.00001, 6, 0.001, 1.0, 0.0, f, g)

y0 = 1.0
varrho_min = varrho2[varrho2> y0*0.001]
eta2_min = eta2[varrho2> y0*0.001]
G = 6.674E-8

print (eta2_min[-1])

k_nr = 3.166E12
k_r = 4.936E14
rho_0 = 3.789E6
lam = ((k_nr * rho_0**(-1.0/3.0))/(4.0 * sp.pi * G))**(1.0/2.0)

dm = simps(varrho_min * eta2_min**2.0, eta2_min)
print(f"Dimensionless mass equals {dm}")


mass = sp.pi * 4.0 * rho_0 * dm * lam**3.0 
solar_mass_conv = 1.989E33
solar_mass = mass / solar_mass_conv
radius = lam * eta2_min[-1]

#print (lam)
#print (solar_mass)
#print (radius)


fig = plt.figure(figsize=(6, 4))
ax = fig.add_subplot(1,1,1)
ax.grid(True)

ax.plot(eta2_min, varrho_min, color="red")
plt.ticklabel_format(style = 'sci', axis = 'y')
plt.xlabel(r'$\eta$')
plt.ylabel(r'$\varrho$')
plt.savefig('bakhtiar_assignment1_C_3.pdf')
plt.show()
