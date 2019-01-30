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



n = 3.0

k = 4.936e14  # cgs units
G = 6.6743e-8  # cgs units
M_sun = 1.989e33 # grams
R_sun = 6.957e10 # cm

def f(eta, varrho, sigma):
    return sigma
    
def g(eta, varrho, sigma):
    return -2.0/eta * sigma - pow(varrho, n)
    
def stop(x, y, z, data):
    return False
    
# solve the Lane-Emden equation using Runge-Kutta

eta, varrho, sigma = phy4910.ode_rk4(0.000000001, 8.0, 0.001, 1.0, 0.0, f, g, stop, None)

phy4910.plot_print(eta, varrho, "A1_B1.pdf", r"$\eta$", r"$\varrho$")

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


# now plot the correct density plot for a white dwarf with radius 8890 km
r_s = 8.89e8 # cm
rho_c = np.power( r_s / surface / 4.853e10, -3.0 ) / 1e6  # in 10^6 g/cc
print(f"The central density for a star with r_s = 8890 km is {rho_c:.3f} x 10^6 g/cc")
rho = rho_c * np.power(varrho, n)
r = r_s / surface * eta / 1e5 # in km

phy4910.plot_print(r, rho, "A1_B2.pdf", r"$r$ (km)", r"$\rho$ ($10^6$ g/cm$^3$)")

# print out the data for later
np.savetxt("A1_B1.dat", np.column_stack((eta, varrho)))
np.savetxt("A1_B2.dat", np.column_stack((r, rho)))


