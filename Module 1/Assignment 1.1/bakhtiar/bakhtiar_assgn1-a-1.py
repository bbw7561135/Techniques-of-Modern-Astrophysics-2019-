#!/usr/bin/env python3

import numpy as np
import phy4910
import matplotlib
import matplotlib.pyplot as plt
from scipy.integrate import simps


n = 1.5

def f(x, y, z):
    return z;
    
def g(x, y, z):
    #print(x, y, z, n)
    return -2.0 / x * z - np.power(y, n);

eta2, rho2, tmp = phy4910.ode_rk4(0.00001, 3.65401, 0.001, 1.0, 0.0, f, g)

varrho = np.power(rho2,n)
dm = simps(varrho * eta2**2, eta2)
print(f"Dimensionless mass equals {dm}")


fig = plt.figure(figsize=(6, 4))
ax = fig.add_subplot(1,1,1)
ax.grid(True)

ax.plot(eta2, rho2, color="red")

plt.xlabel(r'$\eta$ (Distance from center to surface)')
plt.ylabel(r'$\varrho$ (Dimensionless Density)')

plt.savefig('bakhtiar_assignment1_A_1.pdf')
plt.show()
