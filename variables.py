"""
radarmet.variables
=========================

A grouping of functions that calculate a number of radar-derived variables.

Created by Nick Guy.

"""
# HISTORY::
#   5 Feb 2014 - Nick Guy. NOAA/NSSL, NRC (nick.guy@noaa.gov) 
#
# NOTES::
#   Arrays seem to be able to be passed, but make sure they are float arrays
#    (e.g. created with numpy) and not lists
#
# FUNCTIONS::
# reflectivity - Radar reflectivity
# dop_vel - Doppler velocity
# CDR - Circular depolarization ratio
# LDR - Linear depolarization ratio
# ZDR - Differential reflectivity
# ZDP - Reflectivity difference
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
def reflectivity(Pt,G,Tau,lam,bwH,bwV,Lm,Lr,Pr,r,K=0.93):
    """Radar reflectivity.
    
  From Rinehart (1993), Eqn 5.17 (See Eqn 5.14-5.16 also)
 INPUT::
  Pt              = Transmitted power [W]
  G               = Antenna Gain [dB]
  Tau             = Pulse Width [s]
  lam             = Radar wavelength [m]
  bwidth          = Antenna beamwidth [degrees]
  Lm              = Antenna/waveguide/coupler loss [dB]
  Lr              = Receiver loss [dB]
  K               = Dielectric factor [unitless]
  Pr              = Returned power [W]
  r               = Range to target [m]
 OUTPUT::
  Ze               = Radar reflectivity [mm^6/m^3]
 USAGE::
  Ze = reflectivity(Pt,G,Tau,lam,bwH,bwV,Lm,Lr,Pr,r,[K=.93])
 NOTES::
  This routine calls the radar_constant function.
  The default is for a dielectric factor value for water.  This can be 
   changed by the user, e.g. K=0.208 for particle sizes of equivalent melted
   diameters or K=0.176 for particle sizes of equivalent ice spheres.
    """
# HISTORY::
#   3 Feb 2014 - Nick Guy NOAA/NSSL/WRDD, NRC Postdoc
#---------------------------------------
    # Call the radar constant function
    C1 = radar_constant(Pt,G,Tau,lam,bwH,bwV,Lm,Lr)

    # 
    Ze = Pr * r**2 / (C1 * K**2)

    return Ze
#====================================================
def rad_vel(f,lam):
    """Radial velocity.
    
  From Rinehart (1993), Eqn 6.6
 INPUT::
  f               = Frequency shift [Hz]
  lam             = Radar wavelength [m]
 OUTPUT::
  Vr              = Radial velocity [m/s]
 USAGE::
  Vr = rad_vel(f,lam)
    """
# HISTORY::
#   3 Feb 2014 - Nick Guy NOAA/NSSL/WRDD, NRC Postdoc
#---------------------------------------
    Vr = f * lam / 2.

    return Vr
#====================================================
def CDR(Zpar,Zorth):
    """Circular depolarization ratio.  
    
  From Rinehart (1997), Eqn 10.2
 INPUT::
  Zpar          = Reflectivity in the parallel channel [mm^6/m^3]
  Zorth         = Reflectivity in the orthogonal channel [mm^6/m^3]
 OUTPUT::
  CDR           = Circular depolarization ratio [dB]
 USAGE::
  CDR = CDR(Zpar,Zorth)
 NOTES::
  Ensure that both powers have the same units!

  Radars that transmit right-hand circular polarization and receive and
   receive both left- and right-hand circular polarization (using two
   antennas) and acquiring the same pulse.  
  The parallel (orthogonal) component refers to the same 
   (opposite) polarization as transmitted. Non-spericity of hydrometeors
   may be detected (inf long, thin scatterers have CDR = 0 dB, while perfect
   spheres have CDR = -infinity

  Can also use power measurements instead of reflectivity
    """
# HISTORY::
#   3 Feb 2014 - Nick Guy NOAA/NSSL/WRDD, NRC Postdoc
#---------------------------------------
    CDR = 10.* np.log10(Zpar/Zorth)

    return CDR
#====================================================
def LDR(Zh,Zv):
    """Linear depolarization ratio.
    
  From Rinehart (1997), Eqn 10.3
 INPUT::
  Zh            = Horizontal reflectivity [mm^6/m^3]
  Zv            = Vertical reflectivity [mm^6/m^3]
 OUTPUT::
  LDR           = Linear depolarization ratio [dB]
 USAGE::
  LDR = LDR(Zh,Zv)
 NOTES::
  Ensure that both powers have the same units!

  Uses both polarizations in a dual-pol radar from a single pulse  
   Perfect spheres yield LDR => -infinity (though antenna limitations limit
    LDR values to -40 dB for small spheres.
   Long, thin targets, LDR => 0
  Typically -15 > LDR > -35 dB
    """
# HISTORY::
#   3 Feb 2014 - Nick Guy NOAA/NSSL/WRDD, NRC Postdoc
#---------------------------------------
    LDR = 10.* np.log10(Zh/Zv)

    return LDR
#====================================================
def ZDR(Zh,Zv):
    """Differential reflectivity.
    
  From Rinehart (1997), Eqn 10.3 and Seliga and Bringi (1976)
 INPUT::
  Zh            = Horizontal reflectivity [mm^6/m^3]
  Zv            = Vertical reflectivity [mm^6/m^3]
 OUTPUT::
  ZDR           = Differential reflectivity [dB]
 USAGE::
  ZDR = ZDR(Zh,Zv)
 NOTES::
  Ensure that both powers have the same units!

  Alternating horizontally and linearly polarized pulses are averaged.
    """
# HISTORY::
#   3 Feb 2014 - Nick Guy NOAA/NSSL/WRDD, NRC Postdoc 
#---------------------------------------
    ZDR = 10.* np.log10(Zh/Zv)

    return ZDR
#====================================================
def ZDP(Zh,Zv):
    """Reflectivity difference.
    
  From Rinehart (1997), Eqn 10.3
 INPUT::
  Zh            = Horizontal reflectivity [mm^6/m^3]
  Zv            = Vertical reflectivity [mm^6/m^3]
 OUTPUT::
  ZDP           = Reflectivity difference [dB]
 USAGE::
  ZDP = ZDP(Zh,Zv)
 NOTES::
  Ensure that both powers have the same units!

  Alternating horizontally and linearly polarized pulses are averaged. 
    """
# HISTORY::
#   3 Feb 2014 - Nick Guy NOAA/NSSL/WRDD, NRC Postdoc
#---------------------------------------
    if Zh > Zv:
        ZDP = 10.* np.log10(Zh - Zv)
    else:
        print "Zh < Zv !"
        ZDP = np.nan
    return ZDP
#====================================================
def HDR(dBZh,ZDR):
    """Differential reflectivity hail signature.
    
  From Aydin et al. (1986), Eqns 4-5
 INPUT::
  Zh            = Horizontal reflectivity [dBZ]
  ZDR            = Differential reflectivity [dBZ]
 OUTPUT::
  ZDP           = Reflectivity difference [dB]
 USAGE::
  HDR = HDR(dBZh,ZDR)
 NOTES::
  Ensure that both powers have the same units!

  Positive HDR and strong gradients at edges signify ice. The larger HDR,
   the greater likelihood that ice is present.

  Considerations for this equation (see paper for more details): 
   1) Standar error of disdrometer data allowed for
   2) Drop oscillation accounted for based on 50% model of Seliga et al (1984)
   3) Lower (27) bound chose to provide constant Zh ref level
   4) Upper cutoff of 60 (may be too low)

  Picca and Ryzhkof (2012) mention that this does not take into account 
   the hail melting process.  So use at your own risk!
    """
# HISTORY::
#   3 Feb 2014 - Nick Guy NOAA/NSSL/WRDD, NRC Postdoc
#---------------------------------------
    # Set the f(Zdr) based upon observations
    if ZDR <= 0:
        f = 27.
    elif ZDR > 0 and ZDR <= 1.74:
        f = 19. * ZDR + 27.
    elif ZDR > 1.74:
        f = 60.

    # Calculate HDR
    HDR = dBZh - f

    return HDR
#====================================================
