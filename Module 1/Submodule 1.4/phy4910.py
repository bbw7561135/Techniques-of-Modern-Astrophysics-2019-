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

def plot(x, f, xlabel="", ylabel="", labels=None, title=None,  axv=None, noline=False):
    """
    Plot with reasonable settings
    """

    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(1,1,1)
    ax.grid(True)
    for item in f:
        if noline: ax.plot(x, item, ".")
        else: ax.plot(x,item)
    if labels is not None:
        if title is None:
            ax.legend(labels, loc='center left', bbox_to_anchor=(1, 0.5))
        else:
            ax.legend(labels, title=title, loc='center left', bbox_to_anchor=(1, 0.5))
    if axv is not None:
        for item in axv:
            plt.axvline(x=item, linestyle='--')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.show()
    plt.close('all')

def plot_print(x, f, filename, xlabel="", ylabel="", labels=None, title=None, axv=None, dpi=600, noline=False):
    """
    Print the plot with reasonable settings
    """

    fig = plt.figure(figsize=(11, 4), dpi=dpi)
    ax = fig.add_subplot(1,1,1)
    ax.grid(True)
    for item in f:
        if noline: ax.plot(x, item, ".")
        else: ax.plot(x,item)
    if labels is not None:
        if title is None:
            ax.legend(labels, loc='center left', bbox_to_anchor=(1, 0.5))
        else:
            ax.legend(labels, title=title, loc='center left', bbox_to_anchor=(1, 0.5))
    if axv is not None:
        for item in axv:
            plt.axvline(x=item, linestyle='--')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.savefig(filename)
    plt.close('all')

def ode_euler(x_start, x_end, h, y0, z0, f, g):
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



def ode_rk4(x_start, x_end, h, y0, z0, f, g, stop=0, p0=None):
    """
    Solves a coupled pair of ODEs using runge kutta.
    ode_rk4(x_start, x_end, h, y0, z0, f, g, stop=0, p0=None)
    Takes as arguments:
      x_start - starting point for independent coordinate
      x_end - ending point
      h - difference between x_i and x_i+1 (i.e., delat x)
      y0 - initial value for first variable y(x_start)
      z0 - initial value for second variable z(x_start)
      f - function for derivative of first variable (i.e., f = dy/dx)
      g - function for deriviative of second variable (i.e., g = dz/dx)
      stop - a function for a stopping criteria.  Function must take three number (xi, yi, zi).
      data - any data that needs to be passed to the stop function in addition to xi, yi, zi.
      p0 - initial pressure value for fourth variable (optional)

      returns three arrays, x[0,N-1], y[0,N-1], and z[0,N-1].
    """

    x = np.arange(x_start, x_end, h)
    N = len(x)
    y= np.zeros(N)
    y[0] = y0
    z = np.zeros(N)
    z[0] = z0
    if p0 is None:
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
    else:
        p = np.zeros(N)
        p[0] = p0
        for i in range(0, N-1):
            k1 = h * f(x[i], y[i], z[i], p[i])
            l1 = h * g(x[i], y[i], z[i], p[i])
            k2 = h * f(x[i] + 0.5 * h, y[i] + 0.5 * k1, z[i] + 0.5 * l1, p[i])
            l2 = h * g(x[i] + 0.5 * h, y[i] + 0.5 * k1, z[i] + 0.5 * l1, p[i])
            k3 = h * f(x[i] + 0.5 * h, y[i] + 0.5 * k2, z[i] + 0.5 * l2, p[i])
            l3 = h * g(x[i] + 0.5 * h, y[i] + 0.5 * k2, z[i] + 0.5 * l2, p[i])
            k4 = h * f(x[i] + h, y[i] + k3, z[i] + l3, p[i])
            l4 = h * g(x[i] + h, y[i] + k3, z[i] + l3, p[i])

            y[i+1] = y[i] + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
            z[i+1] = z[i] + (l1 + 2.0 * l2 + 2.0 * l3 + l4) / 6.0
            p[i+1] = p[i] * (y[i] / y[i+1]) ** 5

        if(y[i+1] < stop * y0):
            if p0 is None:
                return x[y > 0],y[y > 0],z[y > 0]
            else:
                return x[y > 0],y[y > 0],z[y > 0], p[y > 0]
    if p0 is None:
        return x, y, z
    else:
        return x, y, z, p
