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

image_detection = fits.open("Antennae_Band7_CalibratedData/Antennae_South.CO3_2Line.Clean.pcal1.image.fits", velocity=True)['PRIMARY'].data[0]

for i in range(len(image_detection)):
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    # ax.imshow(transform(image_detection), cmap='gray_r', origin='lower')
    ax.imshow(image_detection[i], origin='lower')
    # plt.show()
    plt.savefig("antennae_south_{}.pdf".format(i))
    plt.close()
