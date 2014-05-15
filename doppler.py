"""
radarmet.doppler
=========================

A grouping of functions that calculate a number of radar characteristics for Doppler radar.

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
#===============================================================
# DEFINE CONTSTANTS
#===============================================================
SLP = 1013.25 # Sea-level Pressure [hPa]
P0 = 1000.  # Reference pressure [hPa]
c = 3e8 # Speed of light [m/s]
Re = 6374000 # Earth's radius [m]
R43 = Re*4./3. # 4/3 Approximation effective radius for standard atmosphere [m]
kBoltz = 1.381e-23 # Boltzmann's constant [ m^2 kg s^-2 K^-1]
#
#===============================================================
# BEGIN FUNCTIONS
#===============================================================
def freq(lam):
    """Frequency given wavelength.
    
 INPUT::
  lam            = Wavelength [m]
 OUTPUT::
  freq             = Frequency [Hz]
 USAGE::
  freq = freq(lam)
    """
# HISTORY::
#   3 Feb 2014 - Nick Guy NOAA/NSSL/WRDD, NRC Postdoc
#---------------------------------------
    freq = c/lam

    return freq
#====================================================
def wavelength(freq):
    """Wavelength given frequency.
    
 INPUT::
  freq            = Frequency [Hz]
 OUTPUT::
  lam             = Wavelength [m]
 USAGE::
  lam = wavelength(freq)
    """
# HISTORY::
#   3 Feb 2014 - Nick Guy NOAA/NSSL/WRDD, NRC Postdoc
#---------------------------------------
    lam = c/freq

    return lam
#====================================================
def pulse_duration(tau):
    """Pulse duration from pulse length.
    
 INPUT::
  tau             = Pulse length [m]
 OUTPUT::
  pDur            = Pulse duration [s]
 USAGE::
  pDur = pulse_duration(tau)
    """
# HISTORY::
#   3 Feb 2014 - Nick Guy NOAA/NSSL/WRDD, NRC Postdoc
#---------------------------------------
    tau = c * pDur/2
    pDur = 2 * tau/c

    return pDur
#====================================================
def pulse_length(pDur):
    """Pulse length from pulse duration.
    
 INPUT::
  pDur            = Pulse duration [s]
 OUTPUT::
  tau             = Pulse length [m]
 USAGE::
  tau = pulse_length(pDur)
    """
# HISTORY::
#   3 Feb 2014 - Nick Guy NOAA/NSSL/WRDD, NRC Postdoc
#---------------------------------------
    tau = c * pDur/2

    return tau
#====================================================
def fmax(PRF):
    """Maximum frequency given PRF.
    
  From Rinehart (1997), Eqn 6.8
 INPUT::
  PRF           = Pulse repetition frequency [Hz]
 OUTPUT::
  fmax          = Maximum frequency [Hz]
 USAGE::
  fmax = fmax(PRF)
    """
# HISTORY::
#   3 Feb 2014 - Nick Guy NOAA/NSSL/WRDD, NRC Postdoc
#---------------------------------------
    fmax = PRF/2.

    return fmax
#====================================================
def Vmax(PRF,lam):
    """Nyquist velocity, or maximum unambiguous Doppler velocity (+ or -).
    
  From Rinehart (1997), Eqn 6.7
 INPUT::
  PRF           = Radar pulse repetition frequency [Hz]
  lam           = Radar wavelength [m]
 OUTPUT::
  Vmax          = Nyquist velocity [m/s], +/-
 USAGE::
  Vmax = Vmax(f,lam)
    """
# HISTORY::
#   3 Feb 2014 - Nick Guy NOAA/NSSL/WRDD, NRC Postdoc
#---------------------------------------
    Vmax = PRF * lam/4.

    return Vmax
#====================================================
def Rmax(PRF):
    """Maximum unamiguous range.
    
  From Rinehart (1997), Eqn 6.11
 INPUT::
  PRF           = Pulse repetition frequency [Hz]
 OUTPUT::
  Rmax          = Maximum unambiguous range [m]
 USAGE::
  Rmax = Rmax(PRF)
    """
# HISTORY::
#   3 Feb 2014 - Nick Guy NOAA/NSSL/WRDD, NRC Postdoc
#---------------------------------------
    Rmax = c / (2. * PRF)

    return Rmax
#====================================================
def Dop_dilemma(In,lam):
    """The "Doppler dilemma" is the fact that both the Nyquist velocity and 
   unambiguous maximum range of the radar are based upon the PRF of the system.
   However, they are inversely proportional, meaning that increasing one 
   requires a decrease in the other.  A trade-off inherent in Doppler radar
   systems.  This relationship allows a solution for one variable given the
   value of the other.

  From Rinehart (1997), Eqn 6.12
 INPUT::
  In            = Nyquist Velocity [m/s] or Maximum unambiguous range [m]
  lam           = Radar wavelength [m]
 OUTPUT::
  Out           = The In that is not used
 USAGE::
  Out = Dop_dilemma(In,lam)
      """
# HISTORY::
#   3 Feb 2014 - Nick Guy NOAA/NSSL/WRDD, NRC Postdoc
#---------------------------------------
    Out = (c * lam/8.)/In

    return Out
#====================================================
##################
# MOBILE PLATFORMS
##################
#====================================================
def Vshift(GS,psi):
    """Adjusted Doppler velocity from a mobile platform.
    
  From Jorgensen (1983), Eqn 2
 INPUT::
  GS            = Gound speed [m/s]
  psi           = Angle between actual azimuth and fore/aft angle [deg]
 OUTPUT::
  Vshift        = Shift in Doppler velocity from mobile aspect [m/s]
 USAGE::
  Vshift = Vshift(GS,psi)
 NOTES::
  In the case of a mobile platform (such as the NOAA P-3 aircraft, the
   Doppler velocity must be adjusted for movement of the scanning platform.

  The fore/aft angle is defined as the angle fore or aft from a plane 
   normal to the direction of motion  
    """
# HISTORY::
#   3 Feb 2014 - Nick Guy NOAA/NSSL/WRDD, NRC Postdoc
#---------------------------------------
    Vshift = GS * np.cos(np.deg2rad(psi))

    return Vshift
#====================================================
def Vmax_dual(lam,PRF1,PRF2):
    """Doppler velocity from dual PRF scheme radar (+ or -).
    
  From Jorgensen (1983), Eqn 2
 INPUT::
  lam           = Radar wavelength [m]
  PRF1          = First Pulse repetition frequency [Hz]
  PRF2          = Second Pulse repetition frequency [Hz]
 OUTPUT::
  Vmax          = Doppler velocity [m/s]
 USAGE::
  Vmax = Vmax_dual(GS,psi)
 NOTES::
  In the case of a mobile platform (such as the NOAA P-3 aircraft, the
   Doppler velocity must be adjusted for movement of the scanning platform.

  The fore/aft angle is defined as the angle fore or aft from a plane 
   normal to the direction of motion  
    """
# HISTORY::
#   3 Feb 2014 - Nick Guy NOAA/NSSL/WRDD, NRC Postdoc
#---------------------------------------
    Vmax = lam / (4* ((1./PRF1) - (1./PRF2)))

    return Vmax
#====================================================

