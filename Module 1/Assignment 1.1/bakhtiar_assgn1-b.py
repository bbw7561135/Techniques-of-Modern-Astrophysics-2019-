#!/usr/bin/env python3

import numpy as np
import phy4910
import matplotlib
import matplotlib.pyplot as plt
from scipy.integrate import simps


n = 3.0

def f(x, y, z):
    return z;
    
def g(x, y, z):
    #print(x, y, z, n)
    return -2.0 / x * z - np.power(y, n);

eta2, rho2, tmp = phy4910.ode_rk4(0.00001, 6.90, 0.001, 1.0, 0.0, f, g)

varrho = np.power(rho2,n)
dm = simps(varrho * eta2**2, eta2)
#print(f"Dimensionless mass equals {dm}")

rhoc= 5.33E7
lam = 4.885199E10* rhoc**(-1/3)

fig = plt.figure(figsize=(6, 4))
ax = fig.add_subplot(1,1,1)
ax.grid(True)

ax.plot(eta2, rho2, color="red")
plt.xlabel(r'$\eta (Dimensionless radius)$')
plt.ylabel(r'$\varrho$ (Dimensionless Density')

plt.savefig('bakhtiar_assignment1_B_1.pdf')


fig = plt.figure(figsize=(6, 4))
ax = fig.add_subplot(1,1,1)
ax.grid(True)
ax.plot(eta2*lam, rho2**n *rhoc, color='blue')
plt.xlabel('radius')
plt.ylabel(r'$\rho$ (Density)')
plt.savefig('Part_B_2.pdf')
plt.show()
