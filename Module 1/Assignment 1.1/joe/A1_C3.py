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
    
#
#  Question C 3 -- with varrho_c = 1 (so rho_c = rho_0)
#
r1, rho1, R1, M1 = phy4910.build_a_white_dwarf(rho_0)
rc('text', usetex=True)
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size':14})
plt.rc('axes', labelsize=16)
plt.rcParams.update({'figure.autolayout': True})

fig = plt.figure(figsize=(6, 4))
ax = fig.add_subplot(1,1,1)
ax.grid(True)
ax.plot(r1, rho1, color="black")
ax.set_xlabel(r"$r$ (R$_\odot$)")
ax.set_ylabel(r"$\rho$ ($10^6$ g/cm$^3$)")
plt.savefig("A1_C3_rho.pdf")
