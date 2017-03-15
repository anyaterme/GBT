#! /usr/bin/env python
import numpy as np
#from astropy.io import fits
import pyfits as fits
import astropy.units as u
import matplotlib.pyplot as plt
T_AMB = 31
T_SYS = 5
T_CAL = 23
DATA = 6

DATA_PATH = "../data/AGBT10B_045_01_SDFITS/"
filename = "AGBT10B_045_01.raw"



hdulist = fits.open("%s/%s.dcr.fits" % (DATA_PATH, filename))
hdulist.close()

hdulist = fits.open("%s/%s.acs.fits" % (DATA_PATH, filename))
hdu = hdulist[1]
columns = hdu.columns
data = hdu.data
hdulist.close()

scans = data.field('SCAN')
scans_keys,counts=np.unique(scans, return_counts=True)

summary = []
for scan in scans_keys:
    #scans_aux = data[np.where(data.field('SCAN') == 11)]
    objects_in_scan = data[np.where(data.field('SCAN') == scan)]
    objects_key, counts = np.unique(objects_in_scan.field('OBJECT'), return_counts=True)
    for object_key in objects_key:
        summary.append("%d\t\t%s" % (scan, object_key))

for line in summary:
    print line

#objects = data.field('OBJECT')
#unique,counts=np.unique(objects, return_counts=True)
#objects = dict(zip(unique, counts))
#for object_name in objects.keys():
#    print object_name




#for index in range(len(data)):
#    freq_resolution = data[index][35] * u.Hz
#    center_freq = data[index][34] * u.Hz
#    init_freq = center_freq - len(data[index][DATA])/2 * freq_resolution
#    end_freq = center_freq + len(data[index][DATA])/2 * freq_resolution
#    x_axis = np.linspace(init_freq.to(u.GHz).value, end_freq.to(u.GHz).value, len(data[index][DATA]))
#    img = (data[index][DATA])
#    plt.plot(x_axis, img)
#plt.xlabel ("GHz")
#plt.ylabel ("Flux")
#plt.show()
