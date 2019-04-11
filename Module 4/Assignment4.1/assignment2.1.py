import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.visualization import MinMaxInterval, LogStretch, make_lupton_rgb, PercentileInterval, AsinhStretch
from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder, IRAFStarFinder, CircularAperture, aperture_photometry
from scipy import ndimage
from astropy import wcs
from astropy.table import Table, hstack
from astropy import units as u
from astropy.coordinates import SkyCoord

# 814, 275, astropy WCS
transform = LogStretch() + PercentileInterval(99.9)
# transform = AsinhStretch() + PercentileInterval(99.9)

image_detection_u = fits.open("hst_12193_04_wfc3_uvis_f395n_drz.fits")['SCI'].data
image_detection_b = fits.open("hst_12193_04_wfc3_uvis_f467m_drz.fits")['SCI'].data
image_detection_y = fits.open("hst_12193_04_wfc3_uvis_f547m_drz.fits")['SCI'].data

rot = 43.75
upperb = 6300
lowerb = 2400
# image_detection_b = ndimage.rotate(image_detection_b, rot)[lowerb: upperb,lowerb:upperb]
# image_detection_y = ndimage.rotate(image_detection_y, rot)[lowerb: upperb,lowerb:upperb]
# image_detection_u = ndimage.rotate(image_detection_u, rot)[lowerb: upperb,lowerb:upperb]

# fig = plt.figure(figsize=(12,7))
# ax1 = fig.add_subplot(131)
# ax1.imshow(transform(image_detection_b), cmap="gray_r", origin='lower')
# ax2 = fig.add_subplot(132)
# ax2.imshow(transform(image_detection_y), cmap="gray_r", origin='lower')
# ax3 = fig.add_subplot(133)
# ax3.imshow(transform(image_detection_u), cmap="gray_r", origin='lower')
# plt.savefig("assignment2.1a-grayscale.pdf")
# plt.colorbar()
# plt.show()

# LUPTON RGB
# image_detection_u = transform(image_detection_u)
# image_detection_b = transform(image_detection_b)
# image_detection_y = transform(image_detection_y)
# mean_u, median_u, std_u = sigma_clipped_stats(image_detection_u, sigma=3.0)
# mean_b, median_b, std_b = sigma_clipped_stats(image_detection_b, sigma=3.0)
# mean_y, median_y, std_y = sigma_clipped_stats(image_detection_y, sigma=3.0)
#
# image_detection_um = image_detection_u - median_u
# image_detection_ym = image_detection_y - median_y
# image_detection_bm = image_detection_b - median_b

# image = make_lupton_rgb(image_detection_ym, image_detection_bm, image_detection_um, minimum=0.002, Q=0, stretch=0.3, filename="m-22.jpg")

image_detection = fits.open("hst_12311_08_wfc3_uvis_total_drz.fits")['SCI'].data

wlist = fits.open("hst_12311_08_wfc3_uvis_total_drz.fits")['SCI']

w = wcs.WCS(wlist.header)

#
# image_detection = ndimage.rotate(image_detection, rot)[lowerb: upperb,lowerb:upperb]
mean, median, std = sigma_clipped_stats(image_detection, sigma=3.0)

iraffind = IRAFStarFinder(fwhm=3.0, threshold=5*std, exclude_border=True)
sources = iraffind(image_detection - median)

positions = (sources['xcentroid'], sources['ycentroid'])
apertures = CircularAperture(positions, r = 4.0)

pc = w.all_world2pix(279.09975,-23.90475,1)

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111, projection = w)
ax.imshow(transform(image_detection), cmap='gray_r', origin='lower')
# ax.colorbar()
apertures.plot(color='blue', lw=1.5, alpha=0.5)
plt.savefig("apertures.pdf")
# print(ax.get_ylim(),ax.get_xlim())
# print(ax.get_yticks(),ax.get_xticks())
# print(ax.get_yaxis_transform(),ax.get_xaxis_transform())
plt.xlabel('Right Ascension')
plt.ylabel('Declination')
# ax.axhline(y=3960, linestyle='--', color='black')
# ax.axvline(x=3300, linestyle='--', color='black')
ax.axhline(y=pc[1], linestyle='--', color='black')
ax.axvline(x=pc[0], linestyle='--', color='black')
plt.savefig('assignment2.1b.pdf')
# plt.show()
plt.close()

image_u = fits.open("hst_12311_08_wfc3_uvis_f275w_drz.fits")['SCI'].data
image_i = fits.open("hst_12311_08_wfc3_uvis_f814w_drz.fits")['SCI'].data

phot_table_u = aperture_photometry(image_u, apertures)
phot_table_i = aperture_photometry(image_i, apertures)
# print(phot_table_u['aperture_sum'][0], phot_table_i['aperture_sum'][0])

m_u = -2.5 * np.log10(phot_table_u['aperture_sum'])
m_i = -2.5 * np.log10(phot_table_i['aperture_sum'])
m_i_u = m_i - m_u

fig = plt.figure(figsize=(6,4))
ax = fig.add_subplot(111)
ax.grid(True)
ax.invert_yaxis()
ax.set_xlabel("U - I")
ax.set_ylabel("U")
ax.plot(m_u - m_i, m_u, ',')
plt.savefig("assignment2-colour_mag.pdf")
plt.close()

pos = np.swapaxes(np.asarray(positions),0,1)
# print(pos.shape)

radec = w.all_pix2world(pos, 1)
# print(radec)
astro_table = hstack([pos, radec, m_u, m_i, m_i_u])
astro_table[0].name = 'X'
astro_table[1].name = 'Y'
astro_table[2].name = 'RA'
astro_table[3].name = 'DEC'
astro_table[4].name = 'm_u'
astro_table[5].name = 'm_i'
astro_table[6].name = 'm_i-m_u'

astro_table.write("assignment2-table.txt", format='ascii.fixed_width', overwrite=True)

# image_detection_u = fits.open("hst_12193_04_wfc3_uvis_f395n_drz.fits")['SCI'].data
# image_detection_b = fits.open("hst_12193_04_wfc3_uvis_f467m_drz.fits")['SCI'].data
# image_detection_y = fits.open("hst_12193_04_wfc3_uvis_f547m_drz.fits")['SCI'].data
#

b1 = w.all_world2pix(279.062738,-23.9038194,1)
b2 = w.all_world2pix(279.0992667, -23.9093056,1)

b1_aperture = CircularAperture(b1, r = 4.0)
b2_aperture = CircularAperture(b2, r = 4.0)

image_detection = fits.open("hst_12311_08_wfc3_uvis_total_drz.fits")['SCI'].data

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111, projection = w)
ax.imshow(transform(image_detection), cmap='gray_r', origin='lower')
b1_aperture.plot(color='blue', lw=1.5, alpha=0.5)
b2_aperture.plot(color='blue', lw=1.5, alpha=0.5)
plt.show()
