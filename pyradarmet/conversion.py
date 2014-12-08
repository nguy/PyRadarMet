"""
pyradarmet.conversion
=========================

A grouping of functions that converts common radar units.

Author:
Nick Guy  NOAA/NSSL, NRC (nick.guy@noaa.gov)

3 Feb 2014 - Created

"""
# NOTES::
#   Arrays seem to be able to be passed, but make sure they are float arrays
#    (e.g. created with numpy) and not lists
#
# FUNCTIONS::
# dBZ2Z - Convert from dBZ to linear Z units
# Z2dBZ - Convert from linear Z units to dBZ
# si2kmh - Convert from SI wind units to km/hr
# si2mph - Convert from SI wind units to miles per hour
# si2kts - Convert from SI wind units to knots
#-------------------------------------------------------------------
# Load the needed packages
import numpy as np
###################
# DEFINE CONTSTANTS
###################

###################
# BEGIN FUNCTIONS
###################

def dBZ2Z(dBZ):
    """Conversion from dBZ (log) units to linear Z units
    
    INPUT::
    -----
    dBZ : float
        logarithmic reflectivity value
 
    OUTPUT::
    ------
    Zlin : float
        linear reflectivity units
    
    USAGE::
    -----
    Zlin = dBZ2Z(dBZ)
    """

    Zlin = 10.**(dBZ/10.)

    return Zlin
    
#============

def Z2dBZ(Zlin):
    """Conversion from linear Z units to dBZ (log) units
    
    INPUT::
    -----
    Zlin : float
        linear reflectivity units
 
    OUTPUT::
    ------
    dBZ : float
        logarithmic reflectivity value
    
    USAGE::
    -----
    dBZ = Z2dBZ(Zlin)
    """

    dBZ = 10. * np.log10(Zlin)

    return dBZ
    
#=============

def si2kmh(SI):
    """Conversion from SI wind units to km/hr
    
    INPUT::
    -----
    SI : float
        Wind in SI units (m/s)
 
    OUTPUT::
    ------
    kmh: float
        Wind in km/hr
    
    USAGE::
    -----
    Ukmh = si2kmh(Usi)
    """

    kmh = SI * 3600. / 1000.
    
    return kmh
    
#=============

def si2mph(SI):
    """Conversion from SI wind units to miles/hr
    
    INPUT::
    -----
    SI: float
        Wind in SI units (m/s)
 
    OUTPUT::
    ------
    mph: float
        Wind in miles per hour
    
    USAGE::
    -----
    Umph = si2mph(Usi)
    """

    mph = SI * 0.62137 / 1000. * 3600.
    
    return mph
    
#=============

def si2kts(SI):
    """Conversion from SI wind units to knots
    
    INPUT::
    -----
    SI: float
        Wind in SI units (m/s)
 
    OUTPUT::
    ------
    kts: float
        Wind in knots
    
    USAGE::
    -----
    Ukts = si2kts(Usi)
    """

    kts = SI * 0.51
    
    return kts
    
#=============

def kmh2si(kmh):
    """Conversion from km/hr to SI wind units
    
    INPUT::
    -----
    kmh: float
        Wind in km/hr
 
    OUTPUT::
    ------
    SI: float
        Wind in SI units (m/s)
    
    USAGE::
    -----
    Ukmh = si2mph(Usi)
    """

    SI = kmh * 1000. / 3600.
    
    return SI
    
#============

def mph2si(mph):
    """Conversion from miles/hr to SI wind units 
    
    INPUT::
    -----
    mph: float
        Wind in miles per hour
 
    OUTPUT::
    ------
    SI: float
        Wind in SI units (m/s)
    
    USAGE::
    -----
    Umph = mph2si(Usi)
    """

    SI = mph * 1000. / (0.62137 * 3600.)
    
    return SI
    
#============

def kts2si(kts):
    """Conversion from knots to SI wind units
    
    INPUT::
    -----
    kts: float
        Wind in knots
 
    OUTPUT::
    ------
    SI: float
        Wind in SI units (m/s)
    
    USAGE::
    -----
    Ukts = si2mph(Usi)
    """

    SI = kts / 0.51
    
    return SI
    
#============


