"""
pyradarmet.geometry
=========================

A grouping of functions that calculate a number of radar geometry characteristics.

Author:
Nick Guy  NOAA/NSSL, NRC (nick.guy@noaa.gov)

3 Feb 2014 - Created

"""
# NOTES::
#   Arrays seem to be able to be passed, but make sure they are float arrays
#    (e.g. created with numpy) and not lists
#
# FUNCTIONS::
# r_effective - Effective radius
# half_power_radius - Half-power beam radius
# ray_height - Height of a ray
# sample_vol_ideal - Ideal Sample volume
# sample_vol_gauss - Sample volume from Gaussian-distributed beam
# range_correct - Range distance corrected for elevation angle "loss" along ground
# beam_block_frac - Partial Beam blockage fraction 
#-------------------------------------------------------------------
# Load the needed packages
import numpy as np

###################
# DEFINE CONTSTANTS
###################
SLP = 1013.25 # Sea-level Pressure [hPa]
P0 = 1000.  # Reference pressure [hPa]
c = 3e8 # Speed of light [m/s]
Re = 6371000 # Earth's average radius [m] assuming sphericity 
R43 = Re*4./3. # 4/3 Approximation effective radius for standard atmosphere [m]
kBoltz = 1.381e-23 # Boltzmann's constant [ m^2 kg s^-2 K^-1]

# Earth radius taken according to International Union of Geodesy and Geophysics
###################
# BEGIN FUNCTIONS
###################

def r_effective(dNdH=-39e-6):
    """
    Effective radius calculation.
    
    Rinehart (1997), Eqn 3.9, solved for R'
    
    INPUT::
    -----
    dNdH : float
        Refraction [N x10^-6/km]
    
    OUTPUT::
    -----
    R1 : float
        Effective radius [m]
    
    USAGE::
    -----
    R1 = r_effective([dNdH=-39.])
    
    NOTES::
    -----
    Effective radius of earth given a refraction.  If no refraction is given
       a "standard atmosphere" is assumed, the valued needed to have straight
       radar rays.
    """

    # Convert earth's radius to km for common dN/dH values and then
    # multiply by 1000 to return radius in meters
    R1 = (1. / ((1/(Re/1000.)) + (dNdH))) * 1000.

    return R1
    
#=============

def half_power_radius(r, bwhalf):
    """
    Half-power radius.
    
    Battan (1973), 
    
    INPUT::
    -----
    r : float
        Range [m]
    bwhalf : float
        Half-power beam width [degrees]
    
    OUTPUT::
    -----
    Rhalf : float
        Half-power radius [m]
    
    USAGE::
    -----
    Rhalf = half_power_radius(r,bwhalf)
    """

    # Convert earth's radius to km for common dN/dH values and then
    # multiply by 1000 to return radius in meters
    Rhalf = (r * np.deg2rad(bwhalf)) / 2.

    return Rhalf
    
#===============

def ray_height(r, elev, H0, R1=R43):
    """
    Center of radar beam height calculation.
    
    Rinehart (1997), Eqn 3.12, Bech et al. (2003) Eqn 3
    
    INPUT::
    -----
    r : float
        Range from radar to point of interest [m]
    elev : float
        Elevation angle of radar beam [deg]
    H0 : float
        Height of radar antenna [m]
    R1 : float
        Effective radius
    
    OUTPUT::
    -----
    H : float
        Radar beam height [m]
    
    USAGE::
    -----
    H = ray_height(r,elev,H0,[R1=6374000.*4/3])
    
    NOTES::
    -----
    If no Effective radius is given a "standard atmosphere" is assumed, 
       the 4/3 approximation.
    Bech et al. (2003) use a factor ke that is the ratio of earth's radius
       to the effective radius (see r_effective function) and Eqn 4 in B03
    """

    # Convert earth's radius to km for common dN/dH values and then
    # multiply by 1000 to return radius in meters
    Term1 = np.sqrt(r**2 +R1**2 + 2*r*R1*np.sin(np.deg2rad(elev)))
    H = Term1 - R1 + H0

    return H
    
#============

def sample_vol_ideal(r, bwH, bwV, pLength):
    """
    Sample volume (idealized) assuming all power in half-power beamwidths.
    
    From Rinehart (1997), Eqn 5.2
    
    INPUT::
    -----
    r : float
        Distance to sample volume from radar [m]
    bwH : float
        Horizontal beamwidth [deg]
    bwV : float
        Vertical beamwidth deg]
    pLength : float
        Pulse length [m]
    
    OUTPUT::
    -----
    SVol : float
        Sample Volume [m^3]
    
    USAGE::
    -----
    SVol = sample_vol_ideal(r,bwH,bwV,pLength)
    
    NOTES::
    -----
    This form assumes all transmitted energy is in the half-power beamwidths.
     A more realistic solution is found in the sample_vol_gauss function
    """

    SVol = np.pi * (r * np.deg2rad(bwH)/2.) * (r * np.deg2rad(bwV)/2.) * (pLength/2.)
    return SVol
    
#===============

def sample_vol_gauss(r, bwH, bwV, pLength):
    """
    Sample volume assuming transmitted energy in Gaussian beam shape.
    
    From Rinehart (1997), Eqn 5.4
    
    INPUT::
    -----
    r  : float
        Distance to sample volume from radar [m]
    bwH : float
        Horizontal beamwidth [deg]
    bwV : float
        Vertical beamwidth deg]
    pLength : float
        Pulse length [m]
    
    OUTPUT::
    -----
    SVol : float
        Sample Volume [m]
    
    USAGE::
    -----
    SVol = sample_vol_gauss(r,bwH,bwV,pLength)
    
    NOTES::
    -----
    This form assumes a Gaussian beam shape for transmitted energy and is more 
      realistic than the sample_vol_ideal.  Derived by Probert-Jones (1962).
    """

    Numer = np.pi * r**2 * np.deg2rad(bwH) * np.deg2rad(bwV) * pLength
    Denom = 16. * np.log(2)
    
    SVol = Numer / Denom
    return SVol
    
#===============

def range_correct(r, h, E):
    """
    A corrected range from radar that takes into account the "loss" of 
      ground distance because of the radar elevation angle.  This is a 
      cumulative effect at each gate along the ray.

    From CSU Radar Meteorology AT 741 Notes
    
    INPUT::
    -----
    r  : float
        Distance to sample volume from radar [m]
    h : float
        Height of the center of radar volume [m]
    E : float
        Elevation angle [deg]
    
    OUTPUT::
    -----
    rnew : float
        Adjusted range to sample volume [m]
    
    USAGE::
    -----
  rnew = range_correct(r,h,E)
    
    NOTES::
    -----
    This function requires that an array be passed!  If you need just one 
       point create a 2 element array with a begin point.
    This is now set up to only accept a 1D array I believe.  May need to 
       fix this in the future.
    """

    # Calculate the change in height along the ray
    dh1 = h[1:] - h[:-1]
    # Add the 0th place in the ray at the beginning
    dh2 = np.insert(dh1,0,h[0])

    # Calculate the change in distance at each gate
    a90r = np.pi/2.  # 90 degrees in radians
    dr = dh2 / (np.tan(a90r - np.deg2rad(E)))

    # Now calculate the corrected range at each gate
    rnew = r - np.cumsum(dr)

    return rnew
    
#===============

def beam_block_frac(Th, Bh, a):
    """Partial beam blockage fraction.
    
    From Bech et al. (2003), Eqn 2 and Appendix
    
    INPUT::
    -----
    Th : float
        Terrain height [m]
    Bh : float
        Beam height [m]
    a : float
        Half power beam radius [m]
    
    OUTPUT::
    -----
    PBB : float
        Partial beam blockage fraction [unitless]
    
    USAGE::
    -----
    PBB = beam_block_frac(Th,Bh,a)
    
    NOTES::
    -----
    This procedure uses a simplified interception function where no vertical
      gradient of refractivity is considered.  Other algorithms treat this
      more thoroughly.  However, this is accurate in most cases other than
      the super-refractive case.

    See the the half_power_radius function to calculate variable a

    The heights must be the same units!
    """

    # First find the difference between the terrain and height of
    # radar beam (Bech et al. (2003), Fig.3)
    y = Th - Bh

    Numer = (y * np.sqrt(a**2 - y**2)) + (a**2 * np.arcsin(y/a)) + (np.pi * a**2 /2.)

    Denom = np.pi * a**2
    
    PBB = Numer / Denom

    return PBB
    
#==============

