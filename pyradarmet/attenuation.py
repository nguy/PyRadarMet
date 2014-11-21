"""
pyradarmet.attenuation
=========================

A grouping of functions that calculate coefficients that can be used with attenuation calculations.

Author:
Nick Guy  NOAA/NSSL, NRC (nick.guy@noaa.gov)

3 Feb 2014 - Created
"""
# NOTES::
#   Arrays seem to be able to be passed, but make sure they are float arrays
#    (e.g. created with numpy) and not lists
#
#FUNCTIONS::
# abs_coeff - Absorption coefficient
# scat_coeff - Scattering coefficient
# ext_coeff - Extinction coefficient
#-------------------------------------------------------------------
# Load the needed packages
import numpy as np
###################
# DEFINE CONTSTANTS
###################

###################
# BEGIN FUNCTIONS
###################
def abs_coeff(D, lam, m):
    """Absorption coefficient of a spherical particle.

    From Doviak and Zrnic (1993), Eqn 3.14a or Battan (1973), Eqn 6.6
    
    INPUT::
    -----
    D : float
        Particle diameter [m]
    lam : float
        Radar wavelength [m]
    m : float
        Complex refractive index [unitless]
 
    OUTPUT::
    ------
    Qa :
        Absorption coefficient [unitless]
 
    USAGE::
    -----
    Qa = abs_coeff(D,lam,m)
 
    NOTES::
    -----
    An example from Battan (1973) is for water at 0C m=7.14-2.89j for a
       wavelength of 3.21 cm and for ice m=1.78-0.0024j for 
       wavelength range from 1-10 cm.
    See Battan (1973) Ch.4 , Tables 4.1 and 4.2 for values from 
       Gunn and East (1954).
    Also see Doviak and Zrnic (1993), Fig. 3.3 caption.
    """

    Km = (m**2 - 1) / (m**2 + 2)
    Qa = (np.pi**2 * D**3 / lam) * np.imag(-1 * Km)

    return Qa

#=============

def scat_coeff(D, lam, m):
    """Scattering coefficient of a spherical particle.
    
    From Doviak and Zrnic (1993), Eqn 3.14b or Battan (1973), Eqn 6.5
    
    INPUT::
    -----
    D : float
        Particle diameter [m]
    lam : float
        Radar wavelength [m]
    m : float
        Complex refractive index [unitless]
    
    OUTPUT::
    ------
    Qs : float
        Scattering coefficient [unitless]
    
    USAGE::
    -----
    Qs = scat_coeff(D,lam,m)
 
    NOTES::
    -----
    An example from Battan (1973) is for water at 0C m=7.14-2.89j for a
       wavelength of 3.21 cm and for ice m=1.78-0.0024j for 
       wavelength range from 1-10 cm.
    See Battan (1973) Ch.4 , Tables 4.1 and 4.2 for values from 
       Gunn and East (1954).
    Also see Doviak and Zrnic (1993), Fig. 3.3 caption.
    """

    Km = (m**2 - 1) / (m**2 + 2)
    Qs = (2 * np.pi**5 * D**6 / (3 * lam**4) * (np.absolute(Km))**2)
    
    return Qs
    
#=============

def ext_coeff(D, lam, m):
    """Extinction coefficient of a spherical particle.
    
    From Doviak and Zrnic (1993), Eqn 3.14b or Battan (1973), Eqn 6.5
    
    INPUT::
    D : float
        Particle diameter [m]
    lam : float
        Radar wavelength [m]
    m : float
        Complex refractive index [unitless]
    
    OUTPUT::
    ------
    Qe : float
        Scattering coefficient [unitless]
 
    USAGE::
    -----
    Qe = ext_coeff(D,lam,m)
 
    NOTES::
    -----
    The default is for a dielectric factor value for water.   
       This can be changed by the user, e.g. K=0.208, diameters 
       for particle sizes of equivalent melted or K=0.176 
       for particle sizes of equivalent ice spheres.
    """
    
    Qa = abs_coeff(D, lam, m)
    Qs = scat_coeff(D, lam, m)
    Qe = Qa + Qs

    return Qe
    
#=============
