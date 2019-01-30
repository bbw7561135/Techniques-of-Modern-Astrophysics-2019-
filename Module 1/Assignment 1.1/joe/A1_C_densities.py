import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import phy4910

rho_0 = 3.789e6  # cgs units
lam =   1.557e8  # cm
M_sun = 1.989e33 # grams
R_sun = 6.957e10 # cm
G = 6.6743e-8  # cgs units
  
  
# Let's build a few white dwarfs for fun -- I want to compare their densities directly.
r1, rho1, R1, M1 = phy4910.build_a_white_dwarf(1e4)
r2, rho2, R2, M2 = phy4910.build_a_white_dwarf(1e6)
r3, rho3, R3, M3 = phy4910.build_a_white_dwarf(1e8)
r4, rho4, R4, M4 = phy4910.build_a_white_dwarf(1e10)

rc('text', usetex=True)
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size':14})
plt.rc('axes', labelsize=16)
plt.rcParams.update({'figure.autolayout': True})

fig = plt.figure(figsize=(6, 4))
ax = fig.add_subplot(1,1,1)
ax.grid(True)
ax.set_xlim(8e-5, 4e-2)
ax.set_ylim(8e-4, 2e4)
ax.loglog(r1, rho1, color="black", label=f"{M1:.3f} Msun")
ax.loglog(r2, rho2, color="red", label=f"{M2:.3f} Msun")
ax.loglog(r3, rho3, color="blue", label=f"{M3:.3f} Msun")
ax.loglog(r4, rho4, color="green", label=f"{M4:.3f} Msun")
ax.legend()
ax.set_xlabel(r"$r$ (R$_\odot$)")
ax.set_ylabel(r"$\rho$ ($10^6$ g/cm$^3$)")
plt.savefig("A1_C_densities.pdf")

