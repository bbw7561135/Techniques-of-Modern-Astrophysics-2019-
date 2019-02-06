import numpy as np
import scipy as sp
import scipy.integrate
import phy4910
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size':14})
plt.rc('axes', labelsize=14)
plt.rcParams.update({'figure.autolayout': True})

filename = "assignment2a2.pdf"

G = 6.674e-8
M_sun = 1.989e33
delta = 0.01 # integration stop point
gamma = 5.0/3.0
R_sun = 6.95508e8
mu_s = 1e-5
mu_star = 4.5
d_eta0 = 0
eta0 = 44.5
stepsize = 0.01
M = 4.5 * M_sun
R = 44.5 * R_sun
# tau units are 23 minutes, 1e3 ~ 10 days

Peq = mu_star * mu_s / eta0 ** 4
# print(Peq)

def f(x, y, z, p):
    return z

def g(x, y, z, p):
    return -mu_star / (y ** 2) + (y ** 2 * p) / mu_s

tau, eta, nu, p = phy4910.ode_rk4(0, 1e3, stepsize, eta0, d_eta0, f, g, delta, Peq)
tau2, eta2, nu2, p2 = phy4910.ode_rk4(0, 1e3, stepsize, eta0, d_eta0, f, g, delta, Peq*1.1)
tau3, eta3, nu3, p3 = phy4910.ode_rk4(0, 1e3, stepsize, eta0, d_eta0, f, g, delta, Peq*1.2)
tau4, eta4, nu4, p4 = phy4910.ode_rk4(0, 1e3, stepsize, eta0, d_eta0, f, g, delta, Peq*1.3)
tau5, eta5, nu5, p5 = phy4910.ode_rk4(0, 1e3, stepsize, eta0, d_eta0, f, g, delta, Peq*0.9)

labels = []
labels.append("{:0.2e}".format(Peq))
labels.append("{:0.2e}".format(Peq*1.1))
labels.append("{:0.2e}".format(Peq*1.2))
labels.append("{:0.2e}".format(Peq*1.3))
labels.append("{:0.2e}".format(Peq*0.9))
tlabels = tuple(labels)

phy4910.plot_print(tau, [eta, eta2, eta3, eta4, eta5], filename, "$\\tau [23 min]$", "$\\eta [R_\\odot]$", tlabels, title="p0")
