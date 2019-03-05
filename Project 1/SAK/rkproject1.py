#!/usr/bin/env python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from decimal import Decimal


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

    fig = plt.figure(figsize=(4, 3), dpi=200)
    ax = fig.add_subplot(1,1,1)
    ax.plot(x, f, color="black")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.savefig(filename)



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

def ode_rk4(x_start, x_end, h, P0, M0, L0, T0, f, g, c, t):

	x = np.arange(x_start, x_end, h)
	N = len(x)
	P = np.zeros(N)
	P[0] = P0
	M = np.zeros(N)
	M[0] = M0
	L = np.zeros(N)
	L[0] = L0
	T = np.zeros(N)
	T[0] = T0


	for i in range(0, N-1):
		k1 = h * f(x[i], P[i], M[i], L[i], T[i])
		l1 = h * g(x[i], P[i], M[i], L[i], T[i])
		m1 = h * c(x[i], P[i], M[i], L[i], T[i])
		n1 = h * t(x[i], P[i], M[i], L[i], T[i])

		k2 = h * f(x[i] + 0.5 * h, P[i] + 0.5 * l1, M[i] + 0.5 * k1, L[i] + 0.5 * m1, T[i] + 0.5 * n1)
		l2 = h * g(x[i] + 0.5 * h, P[i] + 0.5 * l1, M[i] + 0.5 * k1, L[i] + 0.5 * m1, T[i] + 0.5 * n1)
		m2 = h * c(x[i] + 0.5 * h, P[i] + 0.5 * l1, M[i] + 0.5 * k1, L[i] + 0.5 * m1, T[i] + 0.5 * n1)
		n2 = h * t(x[i] + 0.5 * h, P[i] + 0.5 * l1, M[i] + 0.5 * k1, L[i] + 0.5 * m1, T[i] + 0.5 * n1)

		k3 = h * f(x[i] + 0.5 * h, P[i] + 0.5 * l2, M[i] + 0.5 * k2, L[i] + 0.5 * m2, T[i] + 0.5 * n2)
		l3 = h * g(x[i] + 0.5 * h, P[i] + 0.5 * l2, M[i] + 0.5 * k2, L[i] + 0.5 * m2, T[i] + 0.5 * n2)
		m3 = h * c(x[i] + 0.5 * h, P[i] + 0.5 * l2, M[i] + 0.5 * k2, L[i] + 0.5 * m2, T[i] + 0.5 * n2)
		n3 = h * t(x[i] + 0.5 * h, P[i] + 0.5 * l2, M[i] + 0.5 * k2, L[i] + 0.5 * m2, T[i] + 0.5 * n2)

		k4 = h * f(x[i] + h, P[i] + l3, M[i] + k3, L[i] + m3, T[i] + n3)
		l4 = h * g(x[i] + h, P[i] + l3, M[i] + k3, L[i] + m3, T[i] + n3)
		m4 = h * c(x[i] + h, P[i] + l3, M[i] + k3, L[i] + m3, T[i] + n3)
		n4 = h * t(x[i] + h, P[i] + l3, M[i] + k3, L[i] + m3, T[i] + n3)

		P[i+1] = P[i] + (l1 + 2.0 * l2 + 2.0 * l3 + l4) / 6.0
		M[i+1] = M[i] + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
		L[i+1] = L[i] + (m1 + 2.0 * m2 + 2.0 * m3 + m4) / 6.0
		T[i+1] = T[i] + (n1 + 2.0 * n2 + 2.0 * n3 + n4) / 6.0


	return x, P, M, L, T
