
# coding: utf-8

# In[2]:

import numpy as np
import matplotlib
import matplotlib.pyplot as plt


#Commented this out, it was causing errors
#from matplotlib import rc

#rc('text', usetex=True)
#matplotlib.rcParams['mathtext.fontset'] = 'cm'
#matplotlib.rcParams['font.family'] = 'STIXGeneral'
#plt.rcParams.update({'font.size':14})
#plt.rc('axes', labelsize=16)
#plt.rcParams.update({'figure.autolayout': True})

"""
Plot with reasonable settings

-Added logx and logy arguments to allow logarithmic plots

"""
def plot(x, f, xlabel="", ylabel="", logx= False, logy=False):
        
    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(1,1,1)
    ax.grid(True)
    
    if logx==True:
        plt.semilogx()
    if logy==True:
        plt.semilogy()
    
    ax.plot(x, f, color="black")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.show()
    
"""
Print the plot with reasonable settings

-Added logx and logy arguments to allow logarithmic plots


"""
def plot_print(x, f, filename, xlabel="", ylabel="", logx= False, logy=False):
        
    fig = plt.figure(figsize=(4, 3), dpi=200)
    ax = fig.add_subplot(1,1,1)
    
    if logx==True:
        plt.semilogx()
    if logy==True:
        plt.semilogy()
    
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
Empty function used as default in ode_rk4
"""

def Stop(xi,yi,zi):
    return
    
    
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
  
  return_xf - returns the final x value used in the integration
  
"""




def ode_rk4(x_start, x_end, h, y0, z0, f, g, stop = Stop , return_xf = False):
    
    x = np.arange(x_start, x_end, h)
    N = len(x)
    y= np.zeros(N)
    y[0] = y0
    z = np.zeros(N)
    z[0] = z0
    
    for i in range(0, N-1):
        
        if stop(x[i],y[i],z[i]) == True:
            xf = x[i]
            break
        
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
        
        xf = x[i]
        
        
    Final = [x,y,z]
        
    if return_xf==True:
        Final.append(xf)
        
    return Final

"""
Numerically integrates a function of x using the trapezoid method

Takes as arguments:
    x_start - starting point for independent coordinate
    x_end - ending point
    dx - difference between x_i and x_i+1 (i.e., delat x)
    f - function to be integrated
    
    returns two arrays, x[0,n-1], F[0,n-1]
"""

def fint(x_start, x_end, dx, f):
    F = []
    x= np.arange(x_start, x_end, dx)
    n= len(x)
    F.append(dx/2*f(x_start))
    a = 0
    for i in x[1:-1]:
        a = a + 1
        F.append(F[a-1] + dx*f(x_start+a*dx))
    F.append(F[n-2] + dx/2*f(x_end))
    return x,F


"""
Alternative Numerical integration of a function of x using the trapezoid method

Takes as arguments:
    M - an array of values already generated for the function of x
    dx - step size to use
    returns one array, F[0,len(M)-1]
    x values must already be generated
"""

def trapInt (M , dx):
    F = []
    F.append((dx/2)*M[0])
    a = 0
    for i in M[1:-1]:
        a = a + 1
        F.append(F[a-1] + dx*i)
    F.append(F[a] + (dx/2)*M[-1])
    return F






