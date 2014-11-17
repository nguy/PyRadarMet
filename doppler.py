"""
pyradarmet.doppler
=========================

A grouping of functions that calculate a number of radar characteristics for Doppler radar.

Author:
Nick Guy  NOAA/NSSL, NRC (nick.guy@noaa.gov)

3 Feb 2014 - Created

"""
# NOTES::
#   Arrays seem to be able to be passed, but make sure they are float arrays
#    (e.g. created with numpy) and not lists
#
# FUNCTIONS::
# freq - Frequency
# wavelength - Wavelength
# pulse_length - Pulse Length
# pulse_duration - Pulse Duration
# fmax - Maximum frequency
# Vmax - Nyquist or maximum unambiguous velocity
# Rmax - Maximum unambiguous range
# Dop_dilemma - The "Doppler dilemma" either Nyq vel or Rmax yields the other
# Vshift - Shift in Doppler velocity due to moving platform
# Vmax_dual - Maximum unambiguous velocity from dual PRF system
#-------------------------------------------------------------------
# Load the needed packages
import numpy as np
###################
# DEFINE CONTSTANTS
###################
SLP = 1013.25 # Sea-level Pressure [hPa]
P0 = 1000.  # Reference pressure [hPa]
c = 3e8 # Speed of light [m/s]
Re = 6374000 # Earth's radius [m]
R43 = Re*4./3. # 4/3 Approximation effective radius for standard atmosphere [m]
kBoltz = 1.381e-23 # Boltzmann's constant [ m^2 kg s^-2 K^-1]

###################
# BEGIN FUNCTIONS
###################
def freq(lam):
    """Frequency given wavelength.
    
    INPUT::
    -----
    lam : float
        Wavelength [m]
 
    OUTPUT::
    ------
    freq : float
        Frequency [Hz]
    
    USAGE::
    -----
    freq = freq(lam)
    """
    
    freq = c / lam

    return freq
    
#==============

def wavelength(freq):
    """Wavelength given frequency.
    
    INPUT::
    -----
    freq : float
        Frequency [Hz]
 
    OUTPUT::
    ------
    lam : float
        Wavelength [m]
    
    USAGE::
    -----
    lam = wavelength(freq)
    """

    lam = c / freq

    return lam
    
#==============

def pulse_duration(tau):
    """Pulse duration from pulse length.
    
    INPUT::
    -----
    tau : float
        Pulse length [m]
 
    OUTPUT::
    ------
    pDur : float
        Pulse duration [s]
    
    USAGE::
    -----
    pDur = pulse_duration(tau)
    """

    pDur = 2 * tau / c

    return pDur
    
#===============

def pulse_length(pDur):
    """Pulse length from pulse duration.
    
    INPUT::
    -----
    pDur : float
        Pulse duration [s]
 
    OUTPUT::
    ------
    tau : float
        Pulse length [m]
    
    USAGE::
    -----
    tau = pulse_length(pDur)
    """

    tau = c * pDur / 2

    return tau
    
#=============

def fmax(PRF):
    """Maximum frequency given PRF.
    
    From Rinehart (1997), Eqn 6.8
    
    INPUT::
    -----
    PRF : float
        Pulse repetition frequency [Hz]
 
    OUTPUT::
    ------
    fmax : float
        Maximum frequency [Hz]
    
    USAGE::
    -----
    fmax = fmax(PRF)
    """

    fmax = PRF/2.

    return fmax
    
#===============

def Vmax(PRF, lam):
    """Nyquist velocity, or maximum unambiguous Doppler velocity (+ or -).
    
    From Rinehart (1997), Eqn 6.7
    
    INPUT::
    -----
    PRF : float
        Radar pulse repetition frequency [Hz]
    lam : float
        Radar wavelength [m]
 
    OUTPUT::
    ------
    Vmax : float
        Nyquist velocity [m/s], +/-
    
    USAGE::
    -----
    Vmax = Vmax(f,lam)
    """

    Vmax = PRF * lam / 4.

    return Vmax
    
#===============

def Rmax(PRF):
    """Maximum unamiguous range.
    
    From Rinehart (1997), Eqn 6.11
    
    INPUT::
    -----
    PRF : float
        Pulse repetition frequency [Hz]
 
    OUTPUT::
    ------
    Rmax : float
        Maximum unambiguous range [m]
    
    USAGE::
    -----
    Rmax = Rmax(PRF)
    """

    Rmax = c / (2. * PRF)

    return Rmax
    
#==============

def Dop_dilemma(In, lam):
    """The "Doppler dilemma" is the fact that both the Nyquist velocity and 
      unambiguous maximum range of the radar are based upon the PRF of the system.
    However, they are inversely proportional, meaning that increasing one 
      requires a decrease in the other.  A trade-off inherent in Doppler radar
      systems.  This relationship allows a solution for one variable given the
      value of the other.

    From Rinehart (1997), Eqn 6.12
    
    INPUT::
    -----
    In : float
        Nyquist Velocity [m/s] or Maximum unambiguous range [m]
    lam : float
        Radar wavelength [m]
 
    OUTPUT::
    ------
    Out : float
        The In that is not used
    
    USAGE::
    -----
    Out = Dop_dilemma(In,lam)
    """

    Out = (c * lam / 8.) / In

    return Out
    
#==============

##################
# MOBILE PLATFORMS
##################

def Vshift(GS, psi):
    """Adjusted Doppler velocity from a mobile platform.
    
    From Jorgensen (1983), Eqn 2
    
    INPUT::
    -----
    GS : float
        Gound speed [m/s]
    psi : float
        Angle between actual azimuth and fore/aft angle [deg]
 
    OUTPUT::
    ------
    Vshift : float
        Shift in Doppler velocity from mobile aspect [m/s]
    
    USAGE::
    -----
    Vshift = Vshift(GS,psi)
    
    NOTES::
    -----
    In the case of a mobile platform (such as the NOAA P-3 aircraft, the
      Doppler velocity must be adjusted for movement of the scanning platform.

    The fore/aft angle is defined as the angle fore or aft from a plane 
      normal to the direction of motion  
    """

    Vshift = GS * np.cos(np.deg2rad(psi))

    return Vshift
    
#=================

def Vmax_dual(lam, PRF1, PRF2):
    """Doppler velocity from dual PRF scheme radar (+ or -).
    
    From Jorgensen (1983), Eqn 2
    
    INPUT::
    -----
    lam : float
        Radar wavelength [m]
    PRF1 : float
        First Pulse repetition frequency [Hz]
    PRF2 : float
        Second Pulse repetition frequency [Hz]
 
    OUTPUT::
    ------
    Vmax : float
        Doppler velocity [m/s]
    
    USAGE::
    -----
    Vmax = Vmax_dual(GS,psi)
 
    NOTES::
    -----
    In the case of a mobile platform (such as the NOAA P-3 aircraft, the
       Doppler velocity must be adjusted for movement of the scanning platform.

    The fore/aft angle is defined as the angle fore or aft from a plane 
       normal to the direction of motion  
    """

    Vmax = lam / (4 * ((1. / PRF1) - (1. / PRF2)))

    return Vmax
#====================================================

