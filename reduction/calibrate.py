#! /usr/bin/env python
import numpy as np
from astropy.io import fits
import astropy.units as u
import matplotlib.pyplot as plt
T_AMB = 31
T_SYS = 5
T_CAL = 23
DATA = 6


hdulist = fits.open('AGBT10B_045_01.raw.dcr.fits')
hdulist.close()

hdulist = fits.open('AGBT10B_045_01.raw.acs.fits')
hdu = hdulist[1]
img = []

data = hdu.data
img = np.array(data[0][DATA])

print (type(img))
plt.imshow(img)

hdulist.close()
