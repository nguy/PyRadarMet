"""
radarmet.conversion
=========================

A grouping of functions that converts common radar units.

Created by Nick Guy.

"""
# HISTORY::
#   3 Feb 2014 - Nick Guy. NOAA/NSSL, NRC (nick.guy@noaa.gov) 
#
# NOTES::
#   Arrays seem to be able to be passed, but make sure they are float arrays
#    (e.g. created with numpy) and not lists
#
# FUNCTIONS::
# dBZ2Z - Convert from dBZ to linear Z units
# Z2dBZ - Convert from linear Z units to dBZ
#-------------------------------------------------------------------
# Load the needed packages
import numpy as np
#===============================================================
# DEFINE CONTSTANTS
#===============================================================
#
#===============================================================
# BEGIN FUNCTIONS
#===============================================================
def dBZ2Z(dBZ):
    """Conversion from dBZ (log) units to linear Z units
    
 INPUT::
  dBZ           = logarithmic reflectivity value
 OUTPUT::
  Zlin          = linear reflectivity units
 USAGE::
  Zlin = dBZ2Z(dBZ)
    """
# HISTORY::
#   3 Feb 2014 - Nick Guy NOAA/NSSL/WRDD, NRC Postdoc
#---------------------------------------
    Zlin = 10.**(dBZ/10.)

    return Z
#====================================================
def Z2dBZ(Zlin):
    """Conversion from linear Z units to dBZ (log) units
    
 INPUT::
  Zlin          = linear reflectivity value
 OUTPUT::
  dBZ           = logarithmic reflectivity units
 USAGE::
  dBZ = Z2dBZ(Zlin)
    """
# HISTORY::
#   3 Feb 2014 - Nick Guy NOAA/NSSL/WRDD, NRC Postdoc
#---------------------------------------
    dBZ = 10. * np.log10(Zlin)

    return dBZ
#====================================================
