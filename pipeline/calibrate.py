#! /usr/bin/env python
import numpy as np
import time
from astropy.io import fits
from scipy.optimize import curve_fit
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pexpect
import spectralines

def gauss_function(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2)) 

def lorentzian_function(x, I, x0, gamma):
    return I * gamma**2 / ((x - x0)**2 + gamma**2) 

class fit_file:
    def __init__(self, path, smoothing_window_size=0):
        print "Loading file %s..." % path
        start_time = time.time()
        self.path = path
        hdulist = fits.open(path)
        hdu = hdulist[1]
        self.data = hdu.data
        self.lines_summary = None
        hdulist.close()
        #set calibration constants
        self.BB = .0132 # Ruze equation parameter
        self.UNDER_2GHZ_TAU_0 = 0.008
        self.SMOOTHING_WINDOW = smoothing_window_size
        process = pexpect.spawn("gbtidl")
        self.process = process
        prompts = ["GBTIDL -> "]
        self.prompts = prompts
        i=process.expect(prompts)
        print "Access to gbtidl..."
        process.sendline("filein, '%s'" % path)
        print "File read"
#        i=process.expect(prompts, timeout=120)
#        process.sendline("getkanod, 91")
#        i=process.expect(prompts, timeout=120)
#        process.sendline("dcascii, file='aux.txt'")
#        i=process.expect(prompts, timeout=120)
#        process.sendline("exit")


        print "Load Time: ", time.time()-start_time

    def reload_gbtidl(self):
        try:
            self.process.sendline("exit")
        except:
            pass
        process = pexpect.spawn("gbtidl")
        self.process = process
        prompts = ["GBTIDL -> "]
        self.prompts = prompts
        i=process.expect(prompts)
        print "Access to gbtidl..."
        process.sendline("filein, '%s'" % self.path)
        print "File read"

    @property
    def summary(self, force=False):
        """Show fit file summary"""
        if (force):
            self.lines_summary = None
        print "File: ", self.path
        header_summary = "      %s%s%s%s%s%s%s%s%s%s%s" % ("Scan".rjust(8), "Source".rjust(14), "Vel".rjust(15), 'Proc'.rjust(6), 'Seq'.rjust(4), 'RestF'.rjust(15), 'nIF'.rjust(6), 'nInt'.rjust(6), 'nFd'.rjust(6), 'Az'.rjust(8), 'El'.rjust(8))
        print header_summary
        print "=" * (len(header_summary) + 6)
        if (self.lines_summary is None):
            scans = self.data.field('SCAN')
            scans_keys,counts=np.unique(scans, return_counts=True)
            self.lines_summary = []
            index = 0
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
                line = "%s%s%s%s%s%s%s%s%s%s%s" % (str(scan).rjust(8), first_row.field('OBJECT').rjust(14), str(velocity.to(u.km/u.s)).rjust(15), (first_row.field('OBSMODE').split(':')[0]).rjust(6), str(first_row.field('PROCSEQN')).rjust(4), str(rest_freq.to(u.GHz)).rjust(15), str(len(ifnum)).rjust(6), str(nint).rjust(6), str(len(fdnum)).rjust(6), azimuth.rjust(8), elevation.rjust(8))
                print "[%.4d]%s" % (index, line)
                index += 1
                self.lines_summary.append(line)
        else:
            index = 0
            for line in self.lines_summary:
                print "[%.4d]%s" % (index, line)
                index += 1
            print "Total scans: ", len(self.lines_summary)

    def T_cal(self, scan):
        objects_in_scan = self.data[np.where(self.data.field('SCAN') == scan)]
        return objects_in_scan.field('TCAL')

    def record(self, index):
        return self.data[index]

    def source(self, name, scan = None):
        if scan is not None:
            objects_in_scan = self.data[np.where(self.data.field('SCAN') == scan)]
            return objects_in_scan[np.where(self.data.field('OBJECT') == name)]
        return self.data[np.where(self.data.field('OBJECT') == name)]

    def list_sources(self, scan = None):
        if (scan == None):
            objects = self.data.field('OBJECT')
            objects_name,counts=np.unique(objects, return_counts=True)
        else:
            objects = self.data[np.where(self.data.field('SCAN') == scan)]
            objects_name,counts=np.unique(objects.field('OBJECT'), return_counts=True)
        return (objects_name, counts)

    def scan_duration(self, scan):
        objects_in_scan = self.data[np.where(self.data.field('SCAN') == scan)]
        return (sum(objects_in_scan.field('DURATION')))

    def obsmode(self, scan):
        objects = self.data[np.where(self.data.field('SCAN') == scan)]
        obsmode = np.unique(objects.field('OBSMODE'))
        obsmode = obsmode[0].split(':')[0]
        return obsmode

    def source_duration(self, source):
        objects_in = self.data[np.where(self.data.field('OBJECT') == source)]
        scan_keys =np.unique(objects.field('SCAN'))
        result = {}
        for scan in scan_keys:
            result[scan] = sum(objects_in.field('DURATION'))
        return (result)

    def scan(self, scan):
        objects_in_scan = self.data[np.where(self.data.field('SCAN') == scan)]
        return(objects_in_scan)

    def getrec(self, index, unit=u.GHz):
        """Detail important: the data is inverted"""
        record = self.data[index]
        central_freq = record.field('RESTFREQ') * u.Hz
        sky_freq = record.field('OBSFREQ') * u.Hz
        title = "Scan:%s    Vel:%s    Date: %s" % (record.field('SCAN'), str(record.field('VELOCITY')*(u.km / u.s)).rjust(14), record.field('DATE-OBS')[0:10])
        title = "%s    FO: %s    F$_{sky}$: %s" % (title, central_freq.to(unit), sky_freq.to(unit).round(2))
        title = "%s\n%s   T$_{sys}$: %s" % (title, record.field('OBJECT'), str(record.field('TSYS')*u.k))
        data = record.field('DATA')
        data = data[::-1]           # The data are inverted; we have to reverse them
        resolution = record.field('BANDWID') * u.Hz / len(data)
        low_freq = central_freq - len(data)/2 * resolution
        high_freq = central_freq + len(data)/2 * resolution
        freq_axis = np.linspace(low_freq.to(unit).value, high_freq.to(unit).value, len(data))
        fig, ax = plt.subplots(1,1)
        plt.title(r'%s' % title, size=10)
        ax.plot(freq_axis, data)
        ax.set_xlabel(unit)
        ax.set_ylabel('Counts')
        ax.grid(color='r', which='both')
        plt.show(block=False)

    def spectrum(self, scan):
        if (self.procseqn(scan) == 1):
            process = self.process
            prompts = self.prompts
            filename = self.path.split('/')[-1].split('.')[0]
            print "Getting spectrum..."
            if (self.obsmode(scan) == 'Track'):
                process.sendline("getfs, %d" % scan)
            else:
                process.sendline("getkanod, %d" % scan)
            i=process.expect(prompts, timeout=120)
            print "Generating file %s_scan_%d.txt..." % (filename,scan)
            process.sendline("dcascii, file='%s_scan_%d.txt'" % (filename,scan))
            i=process.expect(prompts, timeout=120)
        else:
            print "This scan (%d) is not the first in the sequence" % scan

    def read_spectrum(self, scan, median_filter=7, unit=u.GHz, snr = 7., toscreen=True, tofile = False, showfit = True, path_save="./images"):
        detections = []
        result_detections = []
        if (self.procseqn(scan) == 1):
            from scipy.signal import find_peaks_cwt, medfilt
            import os.path
            import time
            filename = self.path.split('/')[-1].split('.')[0]
            if (not os.path.exists("%s_scan_%d.txt" % (filename, scan))):
                self.spectrum(scan)
                time.sleep(5)
                print "File saved"
            print "Reading spectrum..."
            f = np.loadtxt("%s_scan_%d.txt" % (filename, scan), skiprows=3)
            x = f[:,0]
            y = f[:,1]
            x = np.flip(x, axis=0)
            y = np.flip(y, axis=0)
            #y= medfilt(y,median_filter)

            print "Detecting in scan %s..." % scan

            x80 = x[int(len(x)*0.1):int(len(x)*0.9)]
            y80 = y[int(len(x)*0.1):int(len(x)*0.9)]
            y80 = y80 - np.median(y80)
            widths = np.arange(5,30,2)
            indexes = find_peaks_cwt(y80, widths ,min_snr=snr)
            for i in indexes:
                maxStd = 0
                for width in widths:
                    y_signal = y80[i-(width / 2) : i + (width / 2)]
                    stdSig = np.std(y_signal**2)
                    if (stdSig > maxStd):
                        maxStd = stdSig
                        widthSel = width 
                detections.append([i,widthSel])

            if (showfit):
                ax=plt.subplot(len(detections)+1,1,1)
            else:
                ax=plt.subplot(1,1,1)

            data = self.data[np.where(self.data.field('SCAN') == scan)]
            record = data[0]
            central_freq = record.field('RESTFREQ') * u.Hz
            sky_freq = record.field('OBSFREQ') * u.Hz
            title = "Scan:%s    Vel:%s    Date: %s" % (record.field('SCAN'), str(record.field('VELOCITY')*(u.km / u.s)).rjust(14), record.field('DATE-OBS')[0:10])
            title = "%s    FO: %s    F$_{sky}$: %s" % (title, central_freq.to(unit), sky_freq.to(unit).round(2))
            title = "%s\n%s   T$_{sys}$: %s" % (title, record.field('OBJECT'), str(record.field('TSYS')*u.k))
            print (title)
            resolution = abs(x80[1]-x80[0]) * u.GHz 
            low_freq = x80[0] * resolution
            high_freq = x80[-1] * resolution
            freq_axis = np.linspace(low_freq.to(unit).value, high_freq.to(unit).value, len(x80))
            ax.set_title(title)

            bandwidth = (x80[-1]-x80[0]) * u.GHz
            halpha = spectralines.get_freq(central_freq, bandwidth).to(u.GHz)
            print "Halpha [%s, %s]: %s" % (central_freq, bandwidth, halpha)
            plt.plot(x80 * unit,y80,'b')
            try:
                plt.plot(halpha, 0.9*max(y80), 'go')
            except:
                pass
            if (len(indexes) > 0):
                plt.plot(x80[indexes], y80[indexes],'ro')
                ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%0.2lf"))
                ax.set_xlabel('GHz')
                for index,width in detections:
                    x_signal = x80[index-(width / 2) : index + (width / 2)]
                    y_signal = y80[index-(width / 2) : index + (width / 2)]
                    plt.fill_between(x_signal, y_signal, np.min(y_signal) ,color='y')
                    result_detections.append([x80[index] * unit, resolution])

                img=0
                if (showfit):
                    for index, width in detections:
                        print "Detection in %s" % (x80[i] * unit)
                        width = int(1.5*width)
                        left_limit = np.max([0, index-(width/2)])
                        right_limit = np.min([len(y80)-1, index + (width / 2)])
                        x_signal = x80[left_limit : right_limit]
                        y_signal = y80[left_limit : right_limit]
                        img = img +1
                        params=[np.max(y_signal), np.mean(x_signal),0.5*(x_signal[-1] - x_signal[0])]

                        try:
                            ax=plt.subplot(len(detections)+1,2,2*img +1)
                            popt,pcov = curve_fit(gauss_function, x_signal, y_signal, p0=params)
                            #print "Gauss ", params, popt
                            perfil_gauss = gauss_function(x_signal, *popt)
                            plt.plot(x_signal, y_signal, 'b')
                            plt.plot(x_signal, perfil_gauss, 'g')
                            ax.set_title('Gauss [%s]' % (x80[index] * unit))
                            ax.set_xlabel('GHz')
                            ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%0.5lf"))
                            ax.set_ylabel('Counts')
                            ax.plot(x80[index], y80[index],'rx', label=x80[index] * u.GHz)
                        except Exception as e:
                            print e
                            print "Gauss ", params
                            pass

                        try:
                            ax=plt.subplot(len(detections)+1,2,2*img+2)
                            popt,pcov = curve_fit(lorentzian_function, x_signal, y_signal, p0=params)
                            #print "Lorentz ", params, popt
                            perfil_lorentz = lorentzian_function(x_signal, *popt)
                            plt.plot(x_signal, y_signal, 'b')
                            plt.plot(x_signal, perfil_lorentz, 'r')
                            ax.set_title('Lorentz [%s]' % (x80[index] * unit))
                            ax.set_xlabel('GHz')
                            ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%0.5lf"))
                            ax.set_ylabel('Counts')
                            ax.plot([x80[index]], [y80[index]],'rx', label=x80[index] * u.GHz)
                        except Exception as e:
                            print e
                            print "Lorentz ", params
                            pass
            else:
                print "No detections in scan %d..." % scan
            fig = plt.gcf()
            if (toscreen):
                plt.show()
            if (tofile):
                if (not os.path.exists(path_save)):
                    os.makedirs(path_save)
                fig.savefig("%s/%s_scan_%d.png" % (path_save, filename, scan))
            if(not toscreen):
                plt.close()
        else:
            print "This scan (%d) is not the first in the sequence" % scan

        return (scan, result_detections)

    def procseqn(self, scan):
        data = self.data[np.where(self.data.field('SCAN') == scan)]
        seqn = np.unique(data.field('PROCSEQN'))
        return seqn[0]

    def tsys(self, scan, plnum=0):
        data = self.data[np.where((self.data.field('SCAN') == scan) & (self.data.field('SIG') == 'T') & (self.data.field('PLNUM')==plnum))]
        calon = data[np.where((data.field('CAL') == 'T'))].field('DATA')
        caloff = data[np.where((data.field('CAL') == 'F'))].field('DATA')
        tcal = data.field('TCAL')[::2]
        tsys = []
        for i in range(len(calon)):
            tsys_value = caloff[i] / (calon[i] - caloff[i]) * tcal[i]
            tsys.append(np.mean(tsys_value))

        return tsys



############################################ MAIN ########################################
if __name__ == "__main__":
    DATA_PATH = "../data/AGBT10B_045_01_SDFITS/"
    filename = "AGBT10B_045_01.raw"
    ffits = fit_file("%s/%s.acs.fits" % (DATA_PATH, filename))



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
