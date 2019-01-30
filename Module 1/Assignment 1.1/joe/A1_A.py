import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import phy4910

rc('text', usetex=True)
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size':14})
plt.rc('axes', labelsize=16)
plt.rcParams.update({'figure.autolayout': True})



n = 1.5

k = 3.166e12  # cgs units
G = 6.6743e-8  # cgs units
M_sun = 1.989e33 # grams
R_sun = 6.957e10 # cm

def f(eta, varrho, sigma):
    return sigma
    
def g(eta, varrho, sigma):
    return -2.0/eta * sigma - pow(varrho, n)
    
# solve the Lane-Emden equation using Runge-Kutta

eta, varrho, sigma = phy4910.ode_rk4(0.000000001, 4.0, 0.00001, 1.0, 0.0, f, g)

phy4910.plot_print(eta, varrho, "A1_A1.pdf", r"$\eta$", r"$\varrho$")

# where is the surface?  We can find it by keeping only the positive values of varrho.
# Warning: this assumes that the density function doesn't later on go positive again (i.e. oscillate)
condition = varrho > 0.0
eta = eta[condition]
varrho = varrho[condition]

surface = eta[-1]
print(f"Surface at {surface:.3f}")

# now do the integration
m = np.trapz(np.power(varrho, n) * eta**2, eta)
print(f"The dimensionless mass is {m:.3f}")


# now plot the correct density plot for a 1.0 solar mass white dwarf
lam = 2434 # km
r = lam * eta
rho_c = 4.045 # 10^6 g/cc
rho = rho_c * np.power(varrho, n)

print(f"The surface of a 1.0 solar mass white dwarf is at {r[-1]:.0f} km.")

phy4910.plot_print(r, rho, "A1_A2.pdf", r"$r$ (km)", r"$\rho$ ($10^6$ g/cm$^3$)")

# print out the data for later
np.savetxt("A1_A1.dat", np.column_stack((eta, varrho)))
np.savetxt("A1_A2.dat", np.column_stack((r, rho)))


