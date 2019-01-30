import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc

rc('text', usetex=True)
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size':14})
plt.rc('axes', labelsize=16)
plt.rcParams.update({'figure.autolayout': True})

"""
Plot with reasonable settings
"""
def plot(x, f, xlabel="", ylabel=""):
        
    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(1,1,1)
    ax.grid(True)
    ax.plot(x, f, color="black")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.show()
    
"""
Print the plot with reasonable settings
"""
def plot_print(x, f, filename, xlabel="", ylabel=""):
        
    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(1,1,1)
    ax.plot(x, f, color="black")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.savefig(filename)
    
    
"""
Solves a coupled pair of ODEs using euler
 
Takes as arguments:
  x_start - starting point for independent coordinate
  x_end - ending point
  h - difference between x_i and x_i+1 (i.e., delat x)
  y0 - initial value for first variable y(x_start)
  z_0 - initial value for second variable z(x_start)
  f - function for derivative of first variable (i.e., f = dy/dx)
  g - function for deriviative of second variable (i.e., g = dz/dx)
    
  returns three arrays, x[0,N-1], y[0,N-1], and z[0,N-1].
"""
def ode_euler(x_start, x_end, h, y0, z0, f, g):
    
    x = np.arange(x_start, x_end, h)
    N = len(x)
    y = np.zeros(N)
    y[0] = y0
    z = np.zeros(N)
    z[0] = z0
    
    for i in range(0, N-1):
        k1 = h * f(x[i], y[i], z[i])
        l1 = h * g(x[i], y[i], z[i])
            
        y[i+1] = y[i] + k1
        z[i+1] = z[i] + l1
                
    return x, y, z
    

""" 

Solves a coupled pair of ODEs using runge kutta.
 
Takes as arguments:
  x_start - starting point for independent coordinate
  x_end - ending point
  h - difference between x_i and x_i+1 (i.e., delat x)
  y0 - initial value for first variable y(x_start)
  z_0 - initial value for second variable z(x_start)
  f - function for derivative of first variable (i.e., f = dy/dx)
  g - function for deriviative of second variable (i.e., g = dz/dx)
  stop - a function for a stopping criteria.  Function must take three number (xi, yi, zi).
  data - any data that needs to be passed to the stop function in addition to xi, yi, zi.
  
  returns three arrays, x[0,N-1], y[0,N-1], and z[0,N-1].
"""

def ode_rk4(x_start, x_end, h, y0, z0, f, g, stop, data):
    
    x = np.arange(x_start, x_end, h)
    N = len(x)
    y= np.zeros(N)
    y[0] = y0
    z = np.zeros(N)
    z[0] = z0
    
    for i in range(0, N-1):
        k1 = h * f(x[i], y[i], z[i])
        l1 = h * g(x[i], y[i], z[i])
        k2 = h * f(x[i] + 0.5 * h, y[i] + 0.5 * k1, z[i] + 0.5 * l1)
        l2 = h * g(x[i] + 0.5 * h, y[i] + 0.5 * k1, z[i] + 0.5 * l1)
        k3 = h * f(x[i] + 0.5 * h, y[i] + 0.5 * k2, z[i] + 0.5 * l2)
        l3 = h * g(x[i] + 0.5 * h, y[i] + 0.5 * k2, z[i] + 0.5 * l2)
        k4 = h * f(x[i] + h, y[i] + k3, z[i] + l3)
        l4 = h * g(x[i] + h, y[i] + k3, z[i] + l3)
    
        y[i+1] = y[i] + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
        z[i+1] = z[i] + (l1 + 2.0 * l2 + 2.0 * l3 + l4) / 6.0
        
        if stop(x[i+1], y[i+1], z[i+1], data):
            return x[0:i], y[0:i], z[0:i]
        
    return x, y, z
    
    

'''
Build a white dwarf with a central density specified by rho_c (in units of g/cm^3).
    rho_0 is the density scale, in cgs units
    lam is the distance scale, in cm
    M_sun and R_sun is the mass and radius of the sun, in cgs
'''
def build_a_white_dwarf(rho_c, rho_0 = 3.789e6, lam = 1.557e8, M_sun = 1.989e33, R_sun = 6.957e10):

    print(f"Let's build a white dwarf!")

    varrho_c = rho_c / rho_0
    
    def A(x):
        return -(5.0/9.0) * np.power(x, -4.0/3.0) * np.power(1.0 + np.power(x, 2.0/3.0), -0.5) - (2.0/3.0) * np.power(x, -2.0/3.0) * np.power(1.0 + np.power(x, 2.0/3.0), -1.5) + (1.0/3.0) * np.power(1.0 + np.power(x, 2.0/3.0), -2.5)

    def B(x):
        return (5.0/3.0) * np.power(x, -1.0/3.0) * np.power(1.0 + np.power(x, 2.0/3.0), -0.5) - (1.0/3.0) * np.power(x, 1.0/3.0) * np.power(1.0 + np.power(x, 2.0/3.0), -1.5)

    def f(eta, varrho, sigma):
        return sigma
    
    def g(eta, varrho, sigma):
        return -2.0/eta * sigma - A(varrho) / B(varrho) * sigma*sigma - 1.0/B(varrho) * varrho
    
    def stop(eta, varrho, sigma, data):
        varrho_c = data[0]
        rel_error = data[1]
        if varrho/varrho_c < rel_error:
            print(f"### Stopping ODE solving at radius eta = {eta}")
            return True
        else:
            return False
    
    # solve the Lane-Emden equation using Runge-Kutta

    eta, varrho, sigma = ode_rk4(1e-5, 20.0, 1e-4, varrho_c, 0.0, f, g, stop, [varrho_c, 1e-6])

    # where is the surface?  We can find it by keeping only the positive values of varrho.
    # Warning: this assumes that the density function doesn't later on go positive again (i.e. oscillate)
    condition = varrho > 0.0
    eta = eta[condition]
    varrho = varrho[condition]

    # calculate the dimensionless mass, too.
    m = np.trapz(varrho * eta**2, eta)

    # this is all dimensionless stuff, so let's add the units in.
    # I prefer km to cm, though, for the radius, and 10^ g/cm^3:
    r = eta * lam / R_sun
    rho = varrho * rho_0 / 1e6
    M = (4.0 * np.pi * rho_0 * lam**3 * m) / M_sun

    print(f"  We built a white dwarf:")
    print(f"  * Central density rho_c = {varrho_c * rho_0 / 1e6}*10^6 g/cm^3")
    print(f"  * Radius R_s = {r[-1]:.3f} R_sun")
    print(f"  * Mass M = {M:3f} M_sun")
    
    return r, rho, r[-1], M
    
