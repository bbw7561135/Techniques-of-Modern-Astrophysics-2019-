#!/usr/bin/env python3

import numpy as np
from rkproject1 import *

G = 6.67408e-11 #gravitation constant (m^3/Kg*s^2)
k_b = 1.38e-23 #boltzmans constant (J/K or Kg*m^2/s^2*K)
C = 3.0e8 # speed of light (m/s)
a = 7.565767e-16 #stephen-boltzmann constant (J/m^3*K^4 or Kg/s^2*m*K^4)
sigma = 5.670373e-8 #Wm^-2K^-4
m_H = 1.6737236e-27 # in kg
fpp = 1.0 #assume its 1. PP chain screening factor

sun_p = 2.477e16 #pressure in pascals (Kg/ms^2)
sun_rad = 6.96e8 #meters
sun_mass = 1.99e30 #kg
sun_lum = 3.828e26 # total radiated energy - solar luminosity (kg*m^2/s^3)
sun_temp = 1.56e7 #Kelvin

X = 0.71 #mass fraction of hydrogen
Y = 0.271 # mass fraction of helium
T_0 = 1.56e7
P_0 = 2.477e16
rad = 1.0

Z = 1 - (Y + X) #mass fractionheavier elements
s_rad = rad * sun_rad

mu = 1.0 / ((2.0 * X) + (0.75 * Y) + (0.5 * Z)) #molecular weight mu

def d(x, P, M, L, T):
	return (P - ((1.0/3.0)*a*T**4.0))*(mu*(m_H)/(k_b*T));

def kappa_bar(x, P, M, L, T):
	kappa_bf = 4.34e21*(1.0/0.708)*(1.0/(d(x, P, M, L, T)*(1.0+X))**0.2)*Z*(1.0+X)*d(x, P, M, L, T)/(T**3.5) # m^2/kg
	kappa_ff = 3.68e18*(1.0-Z)*(1.0+X)*d(x, P, M, L, T)/(T**3.5) # m^2/kg
	kappa_es = 0.02*(1.0+X) # m^2/kg
	kappa_H = 0.0
	if T <= 6000 and T >= 3000 and d(x, P, M, L, T)>= 1e-7 and d(x, P, M, L, T) <= 1e-2 and z <= 0.03 and z >= 0.001:
		kappa_H = 7.9e-34*(Z/0.02)*((d(x, P, M, L, T))**0.5)*(T**9.0) # m^2/kg
	return (kappa_bf + kappa_ff + kappa_es + kappa_H)/4.0;

def eps(x, P, M, L, T):
    psi_pp = 1.0 + 1.412e8 * (1.0/X - 1.0) * np.exp(-49.98 * ((T/1.0e6)**(-1.0/3.0)))
    cpp = 1.0 + 0.0123 * ((T/1.0e6)**(1.0/3.0)) + 0.0109 * ((T/1.0e6)**(2.0/3.0)) + (0.000938 * (T/1.0e6))
    cno = 1.0 + 0.0027 * (T/1.0e6)**(1.0/3.0) - (0.00778 * ((T/1.0e6)**(2.0/3.0))) - (0.00149 * (T/1.0e6))
    epp = 0.21 * d(x, P, M, L, T) * X**2.0 * fpp * psi_pp * cpp * ((T/1.0e6)**(-2.0/3.0)) * np.exp(-33.80 * ((T/1.0e6)**(-1.0/3.0)))
    ecno = 8.67e20 * d(x, P, M, L, T) * X * Z/2.0 * cno * ((T/1.0e6)**(-2.0/3.0)) * np.exp(-152.28 * ((T/1.0e6)**(-1.0/3.0)))
    e3a = 50.9 * ((d(x, P, M, L, T))**2.0) * (Y**3.0) * ((T/1.0e8)**(-3.0)) *  np.exp(-44.027 * ((T/1.0e8)**(-1.0)))
    return epp + ecno + e3a;

def f(x, P, M, L, T): # dM/dr
	return 4.0*np.pi*(x**2.0)*d(x, P, M, L, T);

def g(x, P, M, L, T): # dP/dr
	return (-G*M*d(x, P, M, L, T))/(x**2.0);

def c(x, P, M, L, T): # dL/dr
	return 4.0*np.pi*(x**2.0)*d(x, P, M, L, T)*eps(x, P, M, L, T);

def t(x, P, M, L, T): # dT/dr
	return -3.0*kappa_bar(x, P, M, L, T)*d(x, P, M, L, T)*L/((T**3.0)*16.0*a*C*np.pi*x**2.0);
