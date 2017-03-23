#! /usr/bin/env python
import numpy as np
import time
from astropy.io import fits
import astropy.units as u
import matplotlib.pyplot as plt
T_AMB = 31
T_SYS = 5
T_CAL = 23
DATA = 6

DATA_PATH = "../data/AGBT10B_045_01_SDFITS/"
filename = "AGBT10B_045_01.raw"

class fit_file:
    def __init__(self, path):
        print "Loading file %s..." % path
        start_time = time.time()
        self.path = path
        hdulist = fits.open(path)
        hdu = hdulist[1]
        self.data = hdu.data
        hdulist.close()

        scans = self.data.field('SCAN')
        scans_keys,counts=np.unique(scans, return_counts=True)
        self.summary = []
        for scan in scans_keys:
            objects_in_scan = self.data[np.where(self.data.field('SCAN') == scan)]
            first_row = objects_in_scan[0]
            velocity = first_row.field('VELOCITY') * (u.m / u.s)
            rest_freq = first_row.field('RESTFREQ') * u.Hz
            ifnum,counts=np.unique(objects_in_scan.field('IFNUM'), return_counts=True)
            fdnum,counts=np.unique(objects_in_scan.field('FDNUM'), return_counts=True)
            nint = "%d" % len(np.unique(objects_in_scan.field('DATE-OBS')))
            azimuth = "%.2lf" % first_row.field('AZIMUTH')
            elevation = "%.2lf" % first_row.field('ELEVATIO')
            self.summary.append("%s%s%s%s%s%s%s%s%s%s%s" % (str(scan).rjust(8), first_row.field('OBJECT').rjust(14), str(velocity.to(u.km/u.s)).rjust(15), (first_row.field('OBSMODE').split(':')[0]).rjust(6), str(first_row.field('PROCSEQN')).rjust(4), str(rest_freq.to(u.GHz)).rjust(15), str(len(ifnum)).rjust(6), str(nint).rjust(6), str(len(fdnum)).rjust(6), azimuth.rjust(8), elevation.rjust(8)))
        print "Load Time: ", time.time()-start_time

    def Summary(self):
        """Show fit file summary"""
        print "File: ", self.path
        header_summary = "      %s%s%s%s%s%s%s%s%s%s%s" % ("Scan".rjust(8), "Source".rjust(14), "Vel".rjust(15), 'Proc'.rjust(6), 'Seq'.rjust(4), 'RestF'.rjust(15), 'nIF'.rjust(6), 'nInt'.rjust(6), 'nFd'.rjust(6), 'Az'.rjust(8), 'El'.rjust(8))
        print header_summary
        print "=" * (len(header_summary) + 6)
        index = 0
        for line in self.summary:
            print "[%.4d]%s" % (index, line)
            index += 1
        print "Total scans: ", len(self.summary)

    def T_cal(self, scan):
        objects_in_scan = self.data[np.where(self.data.field('SCAN') == scan)]
        return objects_in_scan.field('TCAL')

    def T_sys(self, scan, force=False):
        if ((not(hasattr(self, 't_sys'))) or (force)):
            objects_in_scan = self.data[np.where(self.data.field('SCAN') == scan)]
            objects_in_cal_on = objects_in_scan[np.where(objects_in_scan.field('CAL') == 'T')]
            t_cal = objects_in_scan.field('TCAL')

            data_on = objects_in_cal_on.field('DATA')
            limit = int(len(data_on[0]) * 0.1)
            ref_80 = data_on[:,limit:len(data_on)-limit]
            ref_80_avg_on = np.mean(ref_80)

            objects_in_cal_off = objects_in_scan[np.where(objects_in_scan.field('CAL') == 'F')]
            data_off = objects_in_cal_off.field('DATA')
            limit = int(len(data_off[0]) * 0.1)
            ref_80 = data_off[:,limit:len(data_off)-limit]
            ref_80_avg_off = np.average(ref_80)
            self.t_sys = t_cal * (ref_80_avg_off) / (ref_80_avg_on - ref_80_avg_off) + t_cal*0.5
        return self.t_sys


############################################ MAIN ########################################
if __name__ == "__main__":
    ffits = fit_file("%s/%s.acs.fits" % (DATA_PATH, filename))
#    ffits.Summary()
#    print ffits.T_cal(11)



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
