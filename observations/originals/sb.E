# GBT observing script for 10B-45: Cyanopolyynes + Methyl Cyanide
# observers  V. Hernandez  &  S. Kurtz
# January 2011 observations
#
# Standard Freq Switch observing

#  This is for the 2 high-Dec sources
#  SB should take about 1.25 hr to run

Src1 = "I20126"
Src2 = "DR21-OH"

# Scan duration in seconds
ScanDur=900   

# Load configuration file
execfile("/users/skurtz/10B_45/cyano_config_file")

# Load catalog files
Catalog("/users/skurtz/10B_45/cyano_catalog_file")
Catalog(pointing)
Catalog(fluxcal)

# set-up to 26/29 GHz before pointing/focus
Configure(FSConfigMMF1)

# Pointing/Focus check near the 1st source
AutoPeakFocus(Src1,frequency=28000)

# Reconfigure to undo the AutoPeakFocus config changes
Configure(FSConfigMMF1)

# Balance IF system
Balance()

# Observe for 0.5 hr

Track(Src1,None,ScanDur)
Track(Src2,None,ScanDur)

# set-up to the 36 GHz configuration
# Pointing still OK, but re-focus for 36 GHz
AutoFocus(Src1,frequency=36650)

# Reconfigure to undo the AutoFocus config changes
Configure(FSConfigMMF2)

# Balance IF system
Balance()

# Observe for 0.5 hour

Track(Src1,None,ScanDur)
Track(Src2,None,ScanDur)

# -----------------------------------

