# -*- coding: utf-8 -*-
"""
pyradarmet.geometry
===================

Functions to calculate radar geometry characteristics.

References
----------
Rinehart (1997), Radar for Meteorologists.
Battan (1973), Radar Observations of the Atmosphere.
Bech et al. (2003; JAOT), The Sensitivity of Single Polarization Weather Radar
Beam Blockage Correction to Variability in the Vertical Refractivity Gradient.
"""
import numpy as np


earth_radius = 6371000 # Earth's average radius [m] assuming sphericity
r43 = earth_radius * 4./3. # 4/3 Approximation effective radius for standard atmosphere [m]

# Earth radius taken according to International Union of Geodesy and Geophysics

def r_effective(dndh=-39e-6):
    """
    Effective radius [m] calculation.

    Rinehart (1997), Eqn 3.9, solved for R'

    Parameters
    ----------
    dndh : float
        Refraction [N x10^-6/km]

    Notes
    -----
    Effective radius of earth given a refraction.  If no refraction is given
       a "standard atmosphere" is assumed, the valued needed to have straight
       radar rays.
    """
    # Convert earth's radius to km for common dN/dH values and then
    # multiply by 1000 to return radius in meters
    return (1. / ((1/(earth_radius/1000.)) + (dndh))) * 1000.


def half_power_radius(r, bwhalf):
    """
    Half-power radius [m].

    Battan (1973),

    Parameters
    ----------
    r : float or array
        Range [m]
    bwhalf : float
        Half-power beam width [degrees]
    """
    # Convert earth's radius to km for common dN/dH values and then
    # multiply by 1000 to return radius in meters
    return (np.asarray(r) * np.deg2rad(bwhalf)) / 2.


def ray_height(r, elev, h0, reff=r43):
    """
    Center of radar beam height [m] calculation.

    Rinehart (1997), Eqn 3.12, Bech et al. (2003) Eqn 3

    Parameters
    ----------
    r : float or array
        Range from radar to point of interest [m]
    elev : float
        Elevation angle of radar beam [deg]
    h0 : float
        Height of radar antenna [m]
    reff : float
        Effective radius

    Notes
    -----
    If no Effective radius is given a "standard atmosphere" is assumed,
    the 4/3 approximation.

    Bech et al. (2003) use a factor ke that is the ratio of earth's radius
    to the effective radius (see r_effective function) and Eqn 4 in B03
    """
    # Convert earth's radius to km for common dN/dH values and then
    # multiply by 1000 to return radius in meters
    term1 = (np.sqrt(np.asarray(r)**2 +reff**2 +
             2 * np.asarray(r) * reff * np.sin(np.deg2rad(elev))))
    h = term1 - reff + h0
    return h


def sample_vol_ideal(r, bw_h, bw_v, pulse_length):
    """
    Idealized Sample volume [m^3] assuming all power in half-power beamwidths.

    From Rinehart (1997), Eqn 5.2

    Parameters
    ----------
    r : float or array
        Distance to sample volume from radar [m]
    bw_h : float
        Horizontal beamwidth [deg]
    bw_v : float
        Vertical beamwidth deg]
    pulse_length : float
        Pulse length [m]

    Notes
    -----
    This form assumes all transmitted energy is in the half-power beamwidths.
    A more realistic solution is found in the sample_vol_gauss function
    """
    return (np.pi * (np.asarray(r) * np.deg2rad(bw_h)/2.) * (np.asarray(r) *
            np.deg2rad(bw_v)/2.) * (pulse_length/2.))


def sample_vol_gauss(r, bw_h, bw_v, pulse_length):
    """
    Sample volume [m^3] assuming transmitted energy in Gaussian beam shape.

    From Rinehart (1997), Eqn 5.4

    Parameters
    ----------
    r  : float or array
        Distance to sample volume from radar [m]
    bw_h : float
        Horizontal beamwidth [deg]
    bw_v : float
        Vertical beamwidth deg]
    pulse_length : float
        Pulse length [m]

    Notes
    -----
    This form assumes a Gaussian beam shape for transmitted energy and is more
    realistic than the sample_vol_ideal.  Derived by Probert-Jones (1962).
    """
    Numer = np.pi * np.asarray(r)**2 * np.deg2rad(bw_h) * np.deg2rad(bw_v) * pulse_length
    Denom = 16. * np.log(2)

    SVol = Numer / Denom
    return SVol


def range_correct(r, h, elev):
    """
    A corrected range [m] from radar that takes into account the "loss" of
    ground distance because of the radar elevation angle.  This is a
    cumulative effect at each gate along the ray.

    From CSU Radar Meteorology AT 741 Notes

    Parameters
    ----------
    r  : float or array
        Distance to sample volume from radar [m]
    h : float
        Height of the center of radar volume [m]
    elev : float
        Elevation angle [deg]

    Notes
    -----
    This function requires that an array be passed!  If you need just one
       point create a 2 element array with a begin point.
    This is now set up to only accept a 1D array I believe.  May need to
       fix this in the future.
    """
    # Calculate the change in height along the ray
    dh1 = h[1:] - h[:-1]
    # Add the 0th place in the ray at the beginning
    dh2 = np.insert(dh1, 0, h[0])

    # Calculate the change in distance at each gate
    a90r = np.pi/2.  # 90 degrees in radians
    dr = dh2 / (np.tan(a90r - np.deg2rad(elev)))

    # Now calculate the corrected range at each gate
    rnew = np.asarray(r) - np.cumsum(dr)
    return rnew

#===============

def beam_block_frac(Th, Bh, a):
    """Partial beam blockage fraction. Unitless.

    From Bech et al. (2003), Eqn 2 and Appendix

    Parameters
    ----------
    Th : float
        Terrain height [m]
    Bh : float
        Beam height [m]
    a : float
        Half power beam radius [m]

    Notes
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
