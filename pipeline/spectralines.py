import astropy.units as u
import numpy as np

def read_catalogue(path="./catalogues/splatalogue.tsv"):
    cat = np.genfromtxt(path, delimiter='\t',skip_header=1, dtype=None)
    return cat

def find_nearest(lines, value):
    idx = (np.abs(lines-value)).argmin()
    return lines[idx]

def find_in_range(lines, central_freq, bandwidth):
    lines_in_range = []
    try:
        lines_in_range = lines[np.where(lines <= central_freq + 0.5*bandwidth)]
        lines_in_range = lines_in_range[np.where(lines_in_range >= central_freq - 0.5*bandwidth)]
    except Exception as e:
        print "ERROR in find_in_range: %s" % e
        pass
    return (lines_in_range)

def get_freq(freq, bandwidth= 0, units=u.GHz):
    if (type(freq) != type(1*units)):
        freq = freq * units
    else:
        units = freq.unit
    cat = read_catalogue()
    lines = []
    for row in cat:
        lines.append(row[2])
    lines = lines * u.MHz
    if (bandwidth == 0):
        return (find_nearest(lines, freq).to(units))
    else:
        if (type(bandwidth) != type(1*units)):
           bandwidth = bandwidth * units
        return find_in_range(lines, freq, bandwidth)

