import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.visualization import MinMaxInterval, LogStretch, make_lupton_rgb, PercentileInterval
from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder, IRAFStarFinder, CircularAperture, aperture_photometry

# first open the "detection" image so we can locate all the stars

image_detection = fits.open("hst_12911_01_wfc3_uvis_f775w_sci.fits")['SCI'].data
mean, median, std = sigma_clipped_stats(image_detection, sigma=3.0) 
print(mean, median, std)

iraffind = IRAFStarFinder(fwhm=3.0, threshold=10*std, exclude_border=True)
sources = iraffind(image_detection - median) 

positions = (sources['xcentroid'], sources['ycentroid'])
apertures = CircularAperture(positions, r=4.)

# check!
transform = LogStretch() + PercentileInterval(99.9)
plt.imshow(transform(image_detection), cmap='gray_r', origin='lower')
plt.colorbar()
apertures.plot(color='blue', lw=1.5, alpha=0.5)
plt.show()

# now we can do aperture photometry on our two filtered images
image_u = fits.open("hst_12911_01_wfc3_uvis_f775w_sci.fits")['SCI'].data
image_i = fits.open("hst_12911_01_wfc3_uvis_f467m_sci.fits")['SCI'].data

phot_table_u = aperture_photometry(image_u, apertures)
phot_table_i = aperture_photometry(image_i, apertures)
print(phot_table_u['aperture_sum'][0], phot_table_i['aperture_sum'][0])

m_u = -2.5 * np.log10(phot_table_u['aperture_sum'])
m_i = -2.5 * np.log10(phot_table_i['aperture_sum'])

print(m_u, m_i)

import matplotlib
from matplotlib import rc
rc('text', usetex=True)
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size':14})
plt.rc('axes', labelsize=16)
plt.rcParams.update({'figure.autolayout': True})

fig = plt.figure(figsize=(6, 4))
ax = fig.add_subplot(1,1,1)
ax.grid(True)
ax.invert_yaxis()
ax.plot(m_u - m_i, m_u, ',')
ax.set_xlabel("U - I")
ax.set_ylabel("U")
plt.savefig("colour_mag.pdf")

