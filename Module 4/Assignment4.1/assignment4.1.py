import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.visualization import MinMaxInterval, LogStretch, make_lupton_rgb, PercentileInterval, AsinhStretch, SinhStretch
from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder, IRAFStarFinder, CircularAperture, aperture_photometry
from scipy import ndimage
from astropy import wcs
from astropy.table import Table, hstack
from astropy import units as u
from astropy.coordinates import SkyCoord

# 814, 275, astropy WCS
# transform = LogStretch() + PercentileInterval(99.9)
# transform = AsinhStretch() + PercentileInterval(99.9)
transform = SinhStretch() + PercentileInterval(99.9)

rot = -57
upperb = 2200
lowerb = 5800

image_detection = fits.open("hst_10188_10_acs_wfc_total_drz.fits")['SCI'].data

# image_detection = ndimage.rotate(image_detection, rot)#[lowerb:upperb,lowerb:upperb]

mean, median, std = sigma_clipped_stats(image_detection, sigma=3.0)

image_detection = transform(image_detection) - mean

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)
# ax.imshow(transform(image_detection), cmap='gray_r', origin='lower')
ax.imshow(image_detection, cmap='gray_r', origin='lower')
# plt.show()
plt.close()

image_detection_u = fits.open("hst_10188_10_acs_wfc_f435w_drz.fits")['SCI'].data
image_detection_b = fits.open("hst_10188_10_acs_wfc_f550m_drz.fits")['SCI'].data
# image_detection_b = fits.open("hst_10188_10_acs_wfc_f658n_drz.fits")['SCI'].data
image_detection_y = fits.open("hst_10188_10_acs_wfc_f814w_drz.fits")['SCI'].data

# image_detection_b = ndimage.rotate(image_detection_b, rot)#[lowerb: upperb,lowerb:upperb]
# image_detection_y = ndimage.rotate(image_detection_y, rot)#[lowerb: upperb,lowerb:upperb]
# image_detection_u = ndimage.rotate(image_detection_u, rot)#[lowerb: upperb,lowerb:upperb]

fig = plt.figure(figsize=(12,7))
ax1 = fig.add_subplot(131)
ax1.imshow(transform(image_detection_b), cmap="gray_r", origin='lower')
ax2 = fig.add_subplot(132)
ax2.imshow(transform(image_detection_y), cmap="gray_r", origin='lower')
ax3 = fig.add_subplot(133)
ax3.imshow(transform(image_detection_u), cmap="gray_r", origin='lower')
plt.savefig("assignment4.1a-grayscale.pdf")
# plt.colorbar()
# plt.show()
plt.close()

# LUPTON RGB
image_detection_u = transform(image_detection_u)
image_detection_b = transform(image_detection_b)
image_detection_y = transform(image_detection_y)
mean_u, median_u, std_u = sigma_clipped_stats(image_detection_u, sigma=3.0)
mean_b, median_b, std_b = sigma_clipped_stats(image_detection_b, sigma=3.0)
mean_y, median_y, std_y = sigma_clipped_stats(image_detection_y, sigma=3.0)

image_detection_um = image_detection_u - median_u
image_detection_ym = image_detection_y - median_y
image_detection_bm = image_detection_b - median_b

image = make_lupton_rgb(image_detection_ym, image_detection_bm, image_detection_um, minimum=0.1e-3, Q=-1, stretch=0.5)#, filename="antennae.jpg")

wlist =  fits.open("hst_10188_10_acs_wfc_total_drz.fits")['SCI']
w = wcs.WCS(wlist.header)

radio_north = fits.open("Antennae_Band7_CalibratedData/Antennae_North.CO3_2Line.Clean.pcal1.image.mom0.fits")['PRIMARY']
radio_south = fits.open("Antennae_Band7_CalibratedData/Antennae_South.CO3_2Line.Clean.pcal1.image.mom0.fits")['PRIMARY']

radio_north_wcs = wcs.WCS(radio_north.header)
radio_south_wcs = wcs.WCS(radio_south.header)

radio_north_wcs = radio_north_wcs.sub(2)
radio_south_wcs = radio_south_wcs.sub(2)

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111, projection=w)
ax.imshow(image)
ax.contour(radio_north.data[0][0], levels=[-9, 5, 20, 35, 50, 65, 80, 95], cmap='Reds', transform=ax.get_transform(radio_north_wcs))
ax.contour(radio_south.data[0][0], levels=[-9, 10, 29, 48], cmap='Reds', transform=ax.get_transform(radio_south_wcs))
plt.savefig("antennae_large.pdf")
plt.show()
plt.close()
