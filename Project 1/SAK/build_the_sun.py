#!/usr/bin/env python3

import numpy as np
import rkproject1
from cons_funcs import *
import matplotlib
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

radius, Pressure, Mass, Lum, Temp = rkproject1.ode_rk4(0.001, sun_rad, 10000, 2.477e16, 0.0, 0.0, 1.493969953443072643135e7, f, g, c, t)

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)

ax1.plot(radius/sun_rad, Mass/sun_mass, color="red", label= 'Mass')
ax1.set_xlabel('radius')
ax1.set_ylabel('Mass')
ax1.set_title ('Radius vs Mass')
ax1.grid(True)

ax2.plot(radius/sun_rad, Pressure/sun_p, color="blue", label = 'Pressure')
ax2.set_xlabel('radius')
ax2.set_ylabel('Pressure')
ax2.set_title ('Radius vs Pressure')
ax2.grid(True)

ax3.plot(radius/sun_rad, Lum/sun_lum, color="green", label = 'Luminosity')
ax3.set_xlabel('radius')
ax3.set_ylabel('Luminosity')
ax3.set_title ('Radius vs Luminosity')
ax3.grid(True)

ax4.plot(radius/sun_rad, Temp/sun_temp, color="black", label = 'Temperature')
ax4.set_xlabel('radius')
ax4.set_ylabel('Temperature')
ax4.set_title ('Radius vs Temperature')
ax4.grid(True)

plt.tight_layout(pad = 0.2, w_pad = 0.2, h_pad = 0.2)
plt.savefig('sun.pdf')
plt.show()
