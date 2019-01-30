import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc


(rho_c, R, M) = np.loadtxt("A1_C4.txt", unpack=True)

rc('text', usetex=True)
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size':14})
plt.rc('axes', labelsize=16)
plt.rcParams.update({'figure.autolayout': True})

fig = plt.figure(figsize=(6, 4))
ax = fig.add_subplot(1,1,1)
ax.grid(True)
ax.semilogx(rho_c, M, color="black")
ax.set_xlabel(r"$\rho_c$ (g/cm$^3$)")
ax.set_ylabel(r"$M$ (M$_\odot$)")
ax.semilogx(rho_c, 9.887e29 * np.sqrt(rho_c) / 1.989e33, "--", color="black")
ax.semilogx(rho_c, np.full(len(rho_c), 1.458), "--", color="black")
ax.set_ylim(0.0, 1.55)


ax2 = ax.twinx()
ax2.semilogx(rho_c, R, color="red")
ax2.set_ylabel(r"$R$ (R$_\odot$)", color="red")
ax2.tick_params('y', colors='red')

plt.savefig("A1_C4_plot_rho_c.pdf")

# now plot R and M

fig = plt.figure(figsize=(6, 4))
ax = fig.add_subplot(1,1,1)
ax.grid(True)
ax.plot(M, R, color="black")
ax.set_ylabel(r"$R$ (R$_\odot$)")
ax.set_xlabel(r"$M$ (M$_\odot$)")

plt.savefig("A1_C4_plot_M_R.pdf")

