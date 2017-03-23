# GBT observing script for 10B-45: Cyanopolyynes + Methyl Cyanide
# observers  V. Hernandez  &  S. Kurtz
# January 2011 observations
#
# Standard Freq Switch observing

#  This is for the final 3 molecular cores 
#  Run this AFTER you've run sb.A - sb.E
#  SB should take about 2 hr to run

Src1 = "G30_97MM1"
Src2 = "G31_97MM1"
Src3 = "G34_43MM4"

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

# Observe for 0.75 hr

Track(Src1,None,ScanDur)
Track(Src2,None,ScanDur)
Track(Src3,None,ScanDur)

# set-up to the 36 GHz configuration
# re-focus for 36 GHz  -  and re-point too; it's been 45 min.
AutoPeakFocus(Src1,frequency=36650)

# Reconfigure to undo the AutoFocus config changes
Configure(FSConfigMMF2)

# Balance IF system
Balance()

# Observe for 0.75 hour

Track(Src1,None,ScanDur)
Track(Src2,None,ScanDur)
Track(Src3,None,ScanDur)

# -----------------------------------

