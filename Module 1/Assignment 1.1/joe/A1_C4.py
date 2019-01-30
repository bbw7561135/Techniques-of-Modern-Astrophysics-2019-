import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import phy4910

rho_0 = 3.789e6  # cgs units
lam =   1.557e8  # cm
M_sun = 1.989e33 # grams
R_sun = 6.957e10 # cm
G = 6.6743e-8  # cgs units
    
#
#  Question C 4 -- with rho_c = 10^4 to 10^12 g/cm^3
#

N = 25
rho_c = np.logspace(4, 12, N)
R = np.zeros(N)
M = np.zeros(N)
for i in range(N):
    r, rho, R[i], M[i] = phy4910.build_a_white_dwarf(rho_c[i])
    
# print out the data, this takes too long to run all the time.
np.savetxt("A1_C4.txt", np.column_stack((rho_c, R, M)))
