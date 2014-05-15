"""
radarmet.system
=========================

A grouping of functions that calculate a number of radar system characteristics.

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
# gain_Pratio - Radar gain by power ratio
# freq - Frequency
# wavelength - Wavelength
# pulse_length - Pulse Length
# pulse_duration - Pulse duration
# radar_const - Radar constant
# ant_eff_area - Antenna effective area
# power_target - Power intecepted by a target
# xsec_bscatter_sphere - Backscattering cross-sectional area of sphere
# norm_xsec_bscatter_sphere - Normalized backscattering cross-sectional area of sphere
# size_param - Size parameter of scatterer
# power_return_target - Power returned by spherical target
# thermal_noise - Thermal noise power
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
def gain_Pratio(P1,P2):
    """Antenna gain via power ratio.  
    
  From Rinehart (1997), Eqn 2.1
 INPUT::
  P1            = Power on the beam axis [W]
  P2            = Power from an isotropic antenna [W]
 OUTPUT::
  G             = Gain [dB]
 USAGE::
  G = RadGain_Pratio(P1,P2)
 NOTES::
  Ensure that both powers have the same units!
    """
# HISTORY::
#   3 Feb 2014 - Nick Guy NOAA/NSSL/WRDD, NRC Postdoc
#---------------------------------------
    G = 10.*np.log10(P1/P2)

    return G
#====================================================
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
def pulse_length(pDur):
    """Pulse length given pulse duration.
    
 INPUT::
  pDur            = Pulse duration [s]
 OUTPUT::
  tau             = Pulse length [m]
 USAGE::
  tau = pulse_length(pDur)
 NOTES::
  This equation is only interested in pulses that return to radar, leading
   to the factor of 1/2.
    """
# HISTORY::
#   3 Feb 2014 - Nick Guy NOAA/NSSL/WRDD, NRC Postdoc
#---------------------------------------
    tau = c * pDur/2.

    return tau
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
def radar_const(Pt,G,Tau,lam,bwH,bwV,Lm,Lr):
    """Radar constant.
    
  From CSU Radar Meteorology notes, AT 741
 INPUT::
  Pt              = Transmitted power [W]
  G               = Antenna Gain [dB]
  Tau             = Pulse Width [s]
  lam             = Radar wavelength [m]
  bwH             = Horizontalntenna beamwidth [degrees]
  bwV             = Vertical antenna beamwidth [degrees]
  Lm              = Antenna/waveguide/coupler loss [dB]
  Lr              = Receiver loss [dB]
 OUTPUT::
  C               = Radar constant [unitless]
 USAGE::
  C = radar_constant(Pt,G,Tau,lam,bwidth,Lm,Lr)
    """
# HISTORY::
#   3 Feb 2014 - Nick Guy NOAA/NSSL/WRDD, NRC Postdoc
#---------------------------------------
    # Convert from dB to linear units
    Lmlin = 10**(Lm/10.)
    Lrlin = 10**(Lr/10.)
    Glin = 10**(G/10.)

    # Convert beamwidth to radians
    bwHr = np.deg2rad(bwH)
    bwVr = np.deg2rad(bwV)

    # Calculate the numerator
    Numer = np.pi**3 * c * Pt * Glin**2 * Tau * bwHr * bwVr * Lmlin * Lrlin

    # Calculate the denominator
    Denom = 1024. *np.log(2) * lam**2

    C = Numer/Denom

    return C
#====================================================
def ant_eff_area(G,lam):
    """Antenna effective area.
    
  From Rinehart (1997), Eqn 4.5
 INPUT::
  G               = Antenna Gain [dB]
  lam             = Radar wavelength [m]
 OUTPUT::
  Ae              = Antenna effective area [unitless]
 USAGE::
  Ae = ant_eff_area(G,lam)
    """
# HISTORY::
#   3 Feb 2014 - Nick Guy NOAA/NSSL/WRDD, NRC Postdoc
#---------------------------------------
    # Convert from dB to linear units
    Glin = 10**(G/10.)

    Ae = Glin * lam**2 / (4 * np.pi)

    return Ae
#====================================================
def power_target(Pt,G,Asig,r):
    """Power intercepted by target.
    
  From Rinehart (1997), Eqn 4.3
 INPUT::
  Pt              = Transmitted power [W]
  G               = Antenna gain [dB]
  Asig            = Area of target [m^2]
  r               = Distance to sample volume from radar [m]
 OUTPUT::
  Psig              = Power intecepted by target [m]
 USAGE::
  Psig = power_target(Pt,G,At,r)
    """
# HISTORY::
#   3 Feb 2014 - Nick Guy NOAA/NSSL/WRDD, NRC Postdoc
#---------------------------------------
    # Convert from dB to linear units
    Glin = 10**(G/10.)

    Psig = (Pt * Glin * Asig) / (4 * np.pi * r**2)

    return Psig
#====================================================
def xsec_bscatter_sphere(D,lam,K=0.93):
    """Backscatter cross-sectional area of a sphere using the Rayleigh approximation.
    
  From Rinehart (1997), Eqn 4.9 and 5.7
 INPUT::
  D               = Diamter of targer [m]
  lam             = Radar wavelength [m]
  K               = Dielectric factor [unitless]
 OUTPUT::
  sig             = Backscattering cross-section [m*2]
 USAGE::
  sig = xsec_bscatter_sphere(D,lam, [K=0.93])
 NOTES::
  The Rayleigh approximation is good when the diamter of a spherical particle
   is much smaller than the wavelength of the radar (D/wavelength= 1/16).  This 
   condition leads to the relationship that the area is proportional to the 
   sixth power of the diameter.

  The default is for a dielectric factor value for water.  This can be 
   changed by the user, e.g. K=0.208 for particle sizes of equivalent melted
   diameters or K=0.176 for particle sizes of equivalent ice spheres.
    """
# HISTORY::
#   3 Feb 2014 - Nick Guy NOAA/NSSL/WRDD, NRC Postdoc
#---------------------------------------
    sig = (np.pi**5 * K**2 * D**6) / lam**4

    return sig
#====================================================
def norm_xsec_bscatter_sphere(D,lam,K=0.93):
    """Normalized Backscatter cross-sectional area of a sphere using the Rayleigh approximation.
    
  From Rinehart (1997), Eqn 4.9 and 5.7 and Battan Ch. 4.5
 INPUT::
  D               = Diamter of targer [m]
  lam             = Radar wavelength [m]
  K               = Dielectric factor [unitless]
 OUTPUT::
  sigNorm         = Normalized backscatter cross-section [unitless]
 USAGE::
  sigNorm = norm_xsec_bscatter_sphere(D,lam, [K=0.93])
 NOTES::
  The Rayleigh approximation is good when the diamter of a spherical particle
   is much smaller than the wavelength of the radar (D/wavelength= 1/16).  This 
   condition leads to the relationship that the area is proportional to the 
   sixth power of the diameter.

  The default is for a dielectric factor value for water.  This can be 
   changed by the user, e.g. K=0.208 for particle sizes of equivalent melted
   diameters or K=0.176 for particle sizes of equivalent ice spheres.
    """
# HISTORY::
#  11 Feb 2014 - Nick Guy NOAA/NSSL/WRDD, NRC Postdoc
#---------------------------------------
    # Calculate the cross-sectional backscatter area
    sig = xsec_bscatter_sphere(D,lam,K)

    sigNorm = sig/ (np.pi * (D/2.)**2)
    return sigNorm
#====================================================
def size_param(D,lam):
    """Size parameter calculation.
    
  From Rinehart (1997), Eqn 4.9 and 5.7 and Battan Ch. 4.5
 INPUT::
  D               = Diamter of targer [m]
  lam             = Radar wavelength [m]
 OUTPUT::
  alpha           = Size parameter [unitless]
 USAGE::
  alpha = size_param(D,lam)
 NOTES::
  The size paramter can be used along with the backscattering cross-section to
   distinguish ice and water dielectric characteristics.  For example:
  Alpha < 2 the backscattering cross-section of ice is smaller than water,
   Alpha > 2 the opposite is true due to the fact that absorption in water
   exceeds that in ice.
    """
# HISTORY::
#  11 Feb 2014 - Nick Guy NOAA/NSSL/WRDD, NRC Postdoc
#---------------------------------------
    alpha = 2 * np.pi * D/2. / lam
    return alpha
#====================================================
def power_return_target(Pt,G,lam,sig,r):
    """Power returned y target located at the center of the antenna beam pattern.
    
  From Rinehart (1997), Eqn 4.7
 INPUT::
  Pt              = Transmitted power [W]
  G               = Antenna gain [dB]
  lam             = Radar wavelength [m]
  sig             = Backscattering cross-sectional area of target [m^2]
  r               = Distance to sample volume from radar [m]
 OUTPUT::
  Pr              = Power returned by target [m]
 USAGE::
  Pr = power_return_target(Pt,G,lam,sig,r)
    """
# HISTORY::
#   3 Feb 2014 - Nick Guy NOAA/NSSL/WRDD, NRC Postdoc
#---------------------------------------
    # Convert from dB to linear units
    Glin = 10**(G/10.)

    Pr = (Pt * Glin**2 * lam**2 * sig) / (64 * np.pi**3 * r**4)

    return Pr
#====================================================
def thermal_noise(Bn,Units,Ts=290.):
    """Thermal noise power.
    
  From CSU Radar Meteorology notes, AT741
 INPUT::
  Bn              = Receiver bandwidth [Hz]
  Units           = String of nits desired, can be 'W' or 'dBm'
  Ts              = Reciever noise temperature [K]
 OUTPUT::
  Nt              = Thermal noise power [W or 'dBm']
 USAGE::
  Nt = thermal_noise(Bn,Units,[Ts])
 NOTES::
  Reciever noise temp set to conventional 290K by default
    """
# HISTORY::
#   3 Feb 2014 - Nick Guy NOAA/NSSL/WRDD, NRC Postdoc
#---------------------------------------
    # Calculate the noise, convert if requested
    N = kBoltz * Ts * Bn
    print N
    if Units.upper()=='W':
        Nt = N
    elif Units.upper()=='DBM':
        Nt = 10. * np.log10(N/10**-3)
    else:
        print "Units must be in 'W' or 'dBm'"
        Nt = np.nan
    return Nt
#====================================================
