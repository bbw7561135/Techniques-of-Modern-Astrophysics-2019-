#!/usr/bin/env python3

import numpy as np
import rkproject1
from cons_funcs import *
import matplotlib
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")


X = 0.71 #mass fraction of hydrogen
Y = 0.271 # mass fraction of helium
T_0 = 1.56e7
P_0 = 2.477e16
rad = 1.0

Z = 1 - (Y + X) #mass fractionheavier elements
s_rad = rad * sun_rad

radius, Pressure, Mass, Lum, Temp = rkproject1.ode_rk4(0.001, s_rad , 10000, P_0, 0.0, 0.0, T_0, f, g, c, t)

T_eff = (max(Lum)/(4 * np.pi * sigma * sun_rad**2))**0.25

np.savetxt("star0.txt", np.stack((T_eff, max(Lum))).transpose())
