import numpy as np
import matplotlib.pyplot as plt

# moves one photon through a slab of gas of thickness z_max and maxmimum vertical optical depth tau_max
def move_photon(tau_max = 10.0, z_max = 1.0):
    
    # create photon at origin
    x = 0
    y = 0
    z = 0
    
    while True:
    
        # pick a direction that it scatters into
        # theta goes from 0 to pi, but not uniformly
        cos_theta = 1.0 - 2.0 * np.random.random()
        sin_theta = np.sqrt(1.0 - cos_theta**2)
    
        # phi goes from 0 to 2pi uniformly
        phi = 2.0 * np.pi * np.random.random()
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)

        # move photon in that direction a distance s = tau/t_max
        tau = -np.log(1.0 - np.random.random())
        s = tau/tau_max
    
        # new photon position is
        x += s * sin_theta * cos_phi
        y += s * sin_theta * sin_phi
        z += s * cos_theta
    
        #print(f"{x} {y} {z}")
    
        # did the photon go out the bottom of the slab?
        if z < 0.0:
            x = 0
            y = 0
            z = 0
            #print()
        
        # did the photon leave the top of the slab?
        if z > z_max:
            return x, y, z, cos_theta, sin_theta, cos_phi, sin_phi
            

Nphotons = 10000
Nbins = 20
bins = np.zeros(Nbins)
for i in range(Nphotons):
    x, y, z, cos_theta, sin_theta, cos_phi, sin_phi = move_photon()
    #print i, cos_theta
    pos = int(cos_theta * Nbins)
    bins[pos] += 1
    
angles = np.zeros(Nbins)
intensity = np.zeros(Nbins)
for i in range(Nbins):
    mu = 0.5 * 1.0 / Nbins + i * 1.0/Nbins
    angles[i] = np.arccos(mu)
    intensity[i] = bins[i] / (2.0 * Nphotons * np.cos(angles[i])) * Nbins
    
    print( angles[i], intensity[i])
    
angles = angles * 180.0 / np.pi

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(angles, intensity, color="black", marker='o', linestyle='none')
ax.set_xlabel(r'$\theta$')
ax.set_ylabel(r'Normalized Intensity')


I_theory = 0.5 + 0.75 * np.cos(angles / 180.0 * np.pi)
ax.plot(angles, I_theory, '-')
plt.show()
#plt.savefig("rad_trans_2.pdf")
