"""
pyradarmet.system
=========================

A grouping of functions that calculate a number of radar system characteristics.

Author:
Nick Guy  NOAA/NSSL, NRC (nick.guy@noaa.gov)

3 Feb 2014 - Created

"""
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

def gain_Pratio(P1, P2):
    """
    Antenna gain via power ratio.  
    
    From Rinehart (1997), Eqn 2.1
    
    INPUT::
    -----
    P1 : float
        Power on the beam axis [W]
    P2 : float
        Power from an isotropic antenna [W]
    
    OUPUT::
    -----
    G : float
        Gain [dB]
    
    USAGE::
    -----
    G = RadGain_Pratio(P1,P2)
    
    NOTES::
    -----
    Ensure that both powers have the same units!
    """

    G = 10. * np.log10(P1 / P2)

    return G
    
#==============

def freq(lam):
    """
    Frequency given wavelength.
    
    INPUT::
    -----
    lam : float
        Wavelength [m]
    
    OUPUT::
    -----
    freq : float
        Frequency [Hz]
    
    USAGE::
    -----
    freq = freq(lam)
    """

    freq = c / lam

    return freq
    
#===============

def wavelength(freq):
    """
    Wavelength given frequency.
    
    INPUT::
    -----
    freq : float
        Frequency [Hz]
    
    OUPUT::
    -----
    lam : float
        Wavelength [m]
    
    USAGE::
    -----
    lam = wavelength(freq)
    """

    lam = c / freq

    return lam
    
#==============

def pulse_length(pDur):
    """
    Pulse length given pulse duration.
    
    INPUT::
    -----
    pDur : float
        Pulse duration [s]
    
    OUPUT::
    -----
    tau : float
        Pulse length [m]
    
    USAGE::
    -----
    tau = pulse_length(pDur)
    
    NOTES::
    -----
    This equation is only interested in pulses that return to radar, leading
       to the factor of 1/2.
    """

    tau = c * pDur/2.

    return tau
    
#==============

def pulse_duration(tau):
    """
    Pulse duration from pulse length.
    
    INPUT::
    -----
    tau : float
        Pulse length [m]
    
    OUPUT::
    -----
    pDur : float
        Pulse duration [s]
    
    USAGE::
    -----
    pDur = pulse_duration(tau)
    """

    pDur = 2 * tau / c

    return pDur
    
#===============

def radar_const(Pt, G, Tau, lam, bwH, bwV, Lm, Lr):
    """
    Radar constant.
    
    From CSU Radar Meteorology notes, AT 741
    
    INPUT::
    -----
    Pt : float
        Transmitted power [W]
    G : float
        Antenna Gain [dB]
    Tau : float
        Pulse Width [s]
    lam : float
        Radar wavelength [m]
    bwH : float
        Horizontalntenna beamwidth [degrees]
    bwV : float
        Vertical antenna beamwidth [degrees]
    Lm : float
        Antenna/waveguide/coupler loss [dB]
    Lr : float
        Receiver loss [dB]
    
    OUPUT::
    -----
    C : float
        Radar constant [unitless]
    
    USAGE::
    -----
    C = radar_constant(Pt,G,Tau,lam,bwidth,Lm,Lr)
    """

    # Convert from dB to linear units
    Lmlin = 10**(Lm / 10.)
    Lrlin = 10**(Lr / 10.)
    Glin = 10**(G / 10.)

    # Convert beamwidth to radians
    bwHr = np.deg2rad(bwH)
    bwVr = np.deg2rad(bwV)

    # Calculate the numerator
    Numer = np.pi**3 * c * Pt * Glin**2 * Tau * bwHr * bwVr * Lmlin * Lrlin

    # Calculate the denominator
    Denom = 1024. *np.log(2) * lam**2

    C = Numer/Denom

    return C
    
#============

def ant_eff_area(G, lam):
    """
    Antenna effective area.
    
    From Rinehart (1997), Eqn 4.5
    
    INPUT::
    -----
    G : float
        Antenna Gain [dB]
    lam : float
        Radar wavelength [m]
    
    OUPUT::
    -----
    Ae : float
        Antenna effective area [unitless]
    
    USAGE::
    -----
    Ae = ant_eff_area(G,lam)
    """

    # Convert from dB to linear units
    Glin = 10**(G / 10.)

    Ae = Glin * lam**2 / (4 * np.pi)

    return Ae
    
#============

def power_target(Pt, G, Asig, r):
    """
    Power intercepted by target.
    
    From Rinehart (1997), Eqn 4.3
    
    INPUT::
    -----
    Pt : float
        Transmitted power [W]
    G : float
        Antenna gain [dB]
    Asig : float
        Area of target [m^2]
    r : float
        Distance to sample volume from radar [m]
    
    OUPUT::
    -----
    Psig : float
        Power intecepted by target [m]
    
    USAGE::
    -----
    Psig = power_target(Pt,G,At,r)
    """

    # Convert from dB to linear units
    Glin = 10**(G / 10.)

    Psig = (Pt * Glin * Asig) / (4 * np.pi * r**2)

    return Psig
    
#==============

def xsec_bscatter_sphere(D, lam, K=0.93):
    """
    Backscatter cross-sectional area of a sphere using the Rayleigh approximation.
    
    From Rinehart (1997), Eqn 4.9 and 5.7
    
    INPUT::
    -----
    D : float
        Diamter of targer [m]
    lam : float
        Radar wavelength [m]
    K : float
        Dielectric factor [unitless]
    
    OUPUT::
    -----
    sig : float
        Backscattering cross-section [m*2]
    
    USAGE::
    -----
    sig = xsec_bscatter_sphere(D,lam, [K=0.93])
    
    NOTES::
    -----
    The Rayleigh approximation is good when the diamter of a spherical particle
       is much smaller than the wavelength of the radar (D/wavelength= 1/16).  This 
       condition leads to the relationship that the area is proportional to the 
       sixth power of the diameter.

    The default is for a dielectric factor value for water.  This can be 
       changed by the user, e.g. K=0.208 for particle sizes of equivalent melted
       diameters or K=0.176 for particle sizes of equivalent ice spheres.
    """

    sig = (np.pi**5 * K**2 * D**6) / lam**4

    return sig
    
#==============

def norm_xsec_bscatter_sphere(D, lam, K=0.93):
    """
    Normalized Backscatter cross-sectional area of a sphere using the Rayleigh approximation.
    
    From Rinehart (1997), Eqn 4.9 and 5.7 and Battan Ch. 4.5
    
    INPUT::
    -----
    D : float
        Diamter of targer [m]
    lam : float
        Radar wavelength [m]
    K : float
        Dielectric factor [unitless]
    
    OUPUT::
    -----
    sigNorm : float
        Normalized backscatter cross-section [unitless]
    
    USAGE::
    -----
    sigNorm = norm_xsec_bscatter_sphere(D,lam, [K=0.93])
    
    NOTES::
    -----
    The Rayleigh approximation is good when the diamter of a spherical particle
       is much smaller than the wavelength of the radar (D/wavelength= 1/16).  This 
       condition leads to the relationship that the area is proportional to the 
       sixth power of the diameter.

    The default is for a dielectric factor value for water.  This can be 
       changed by the user, e.g. K=0.208 for particle sizes of equivalent melted
       diameters or K=0.176 for particle sizes of equivalent ice spheres.
    """

    # Calculate the cross-sectional backscatter area
    sig = xsec_bscatter_sphere(D, lam, K)

    sigNorm = sig/ (np.pi * (D/2.)**2)
    return sigNorm
    
#==================

def size_param(D, lam):
    """
    Size parameter calculation.
    
    From Rinehart (1997), Eqn 4.9 and 5.7 and Battan Ch. 4.5
    
    INPUT::
    -----
    D : float
        Diamter of targer [m]
    lam : float
        Radar wavelength [m]
    
    OUPUT::
    -----
    alpha : float
        Size parameter [unitless]
    
    USAGE::
    -----
    alpha = size_param(D,lam)
    
    NOTES::
    -----
    The size paramter can be used along with the backscattering cross-section to
       distinguish ice and water dielectric characteristics.  For example:
       Alpha < 2 the backscattering cross-section of ice is smaller than water,
       Alpha > 2 the opposite is true due to the fact that absorption in water
       exceeds that in ice.
    """

    alpha = 2 * np.pi * D/2. / lam
    return alpha
    
#================

def power_return_target(Pt, G, lam, sig, r):
    """
    Power returned y target located at the center of the antenna beam pattern.
    
    From Rinehart (1997), Eqn 4.7
    
    INPUT::
    -----
    Pt : float
        Transmitted power [W]
    G : float
        Antenna gain [dB]
    lam : float
        Radar wavelength [m]
    sig : float
        Backscattering cross-sectional area of target [m^2]
    r : float
        Distance to sample volume from radar [m]
    
    OUPUT::
    -----
    Pr : float
        Power returned by target [m]
    
    USAGE::
    -----
    Pr = power_return_target(Pt,G,lam,sig,r)
    """

    # Convert from dB to linear units
    Glin = 10**(G/10.)

    Pr = (Pt * Glin**2 * lam**2 * sig) / (64 * np.pi**3 * r**4)

    return Pr
    
#=============

def thermal_noise(Bn, Units, Ts=290.):
    """
    Thermal noise power.
    
    From CSU Radar Meteorology notes, AT741
    
    INPUT::
    -----
    Bn : float
        Receiver bandwidth [Hz]
    Units : float
        String of nits desired, can be 'W' or 'dBm'
    Ts : float
        Reciever noise temperature [K]
    
    OUPUT::
    -----
    Nt : float
        Thermal noise power [W or 'dBm']
    
    USAGE::
    -----
    Nt = thermal_noise(Bn,Units,[Ts])
    
    NOTES::
    -----
    Reciever noise temp set to conventional 290K by default
    """

    # Calculate the noise, convert if requested
    N = kBoltz * Ts * Bn
    
    if Units.upper()=='W':
        Nt = N
    elif Units.upper()=='DBM':
        Nt = 10. * np.log10(N/10**-3)
    else:
        print "Units must be in 'W' or 'dBm'"
        Nt = np.nan
    return Nt
    
#=============

