# -*- coding: utf-8 -*-
"""
pyradarmet.variables
====================

Functions to calculate radar-derived variables.

References
----------
Rinehart (1997), Radar for Meteorologists.
Aydin et al. (1986; JCAM), Remote Sensing of Hail with Dual Linear
Polarization Radar
"""
import numpy as np


def reflectivity(power_t, gain, pulse_width, wavelength, bw_h, bw_v,
                 aloss, rloss, power_return, r, dielectric=0.93):
    """
    Radar reflectivity [mm^6/m^3].

    From Rinehart (1993), Eqn 5.17 (See Eqn 5.14-5.16 also)

    Parameters
    ----------
    power_t : float
        Transmitted power [W]
    gain : float
        Antenna Gain [dB]
    pulse_width : float
        Pulse Width [s]
    wavelength : float
        Radar wavelength [m]
    bwidth : float
        Antenna beamwidth [degrees]
    aloss : float
        Antenna/waveguide/coupler loss [dB]
    rloss : float
        Receiver loss [dB]
    dielectric : float
        Dielectric factor [unitless]
    power_return : float
        Returned power [W]
    r : float or array
        Range to target [m]

    Notes
    -----
    This routine calls the radar_constant function.
    The default is for a dielectric factor value for water.  This can be
        changed by the user, e.g. K=0.208 for particle sizes of equivalent melted
        diameters or K=0.176 for particle sizes of equivalent ice spheres.
    """
    # Call the radar constant function
    C1 = radar_constant(power_t, gain, pulse_width, wavelength, bw_h, bw_v, aloss, rloss)
    return power_return * np.asarray(r)**2 / (C1 * dielectric**2)


def radial_velocity(frequency, wavelength):
    """
    Radial velocity [m/s].

    From Rinehart (1997), Eqn 6.6

    Parameters
    ----------
    frequency : float or array
        Frequency shift [Hz]
    wavelength : float or array
        Radar wavelength [m]

    Notes
    -----
    If arrays are used, either for one okay, or both must be the same length.
    """
    return np.asarray(frequency) * np.asarray(wavelength) / 2.


def cdr(refl_parallel, refl_orthogonal):
    """
    Circular depolarization ratio [dB].

    From Rinehart (1997), Eqn 10.2

    Parameters
    ----------
    refl_parallel : float or array
        Reflectivity in the parallel channel [mm^6/m^3]
    refl_orthogonal : float or array
        Reflectivity in the orthogonal channel [mm^6/m^3]

    Notes
    -----
    Ensure that both powers have the same units!

    Radars that transmit right-hand circular polarization and receive and
       receive both left- and right-hand circular polarization (using two
       antennas) and acquiring the same pulse.
    The parallel (orthogonal) component refers to the same
      (opposite) polarization as transmitted. Non-spericity of hydrometeors
      may be detected (inf long, thin scatterers have CDR = 0 dB, while perfect
     spheres have CDR = -infinity

    Can also use power measurements instead of reflectivity.

    If arrays are used, either for one okay, or both must be the same length.
    """
    return 10. * np.log10(np.asarray(refl_parallel)/np.asarray(refl_orthogonal))


def ldr(z_h, z_v):
    """
    Linear depolarization ratio [dB].

    From Rinehart (1997), Eqn 10.3

    Parameters
    ----------
    z_h : float or array
        Horizontal reflectivity [mm^6/m^3]
    z_v : float or array
        Vertical reflectivity [mm^6/m^3]

    Notes
    -----
    Ensure that both powers have the same units!

    Uses both polarizations in a dual-pol radar from a single pulse.

    Perfect spheres yield LDR => -infinity (though antenna limitations limit
    LDR values to -40 dB for small spheres.
    Long, thin targets, LDR => 0

    Typical values in the range -15 > LDR > -35 dB

    If arrays are used, either for one okay, or both must be the same length.
    """
    return 10. * np.log10(np.asarray(z_h) / np.asarray(z_v))


def zdr(z_h, z_v):
    """
    Differential reflectivity [dB].

    From Rinehart (1997), Eqn 10.3 and Seliga and Bringi (1976)

    Parameters
    ----------
    z_h : float or array
        Horizontal reflectivity [mm^6/m^3]
    z_v : float or array
        Vertical reflectivity [mm^6/m^3]

    Notes
    -----
    Ensure that both powers have the same units!

    Alternating horizontally and linearly polarized pulses are averaged.

    Notes
    -----
    If arrays are used, either for one okay, or both must be the same length.
    """
    return 10. * np.log10(np.asarray(z_h) / np.asarray(z_v))


def zdp(z_h, z_v):
    """
    Reflectivity difference [dB].

    From Rinehart (1997), Eqn 10.3

    Parameters
    ----------
    z_h : float
        Horizontal reflectivity [mm^6/m^3]
    z_v : float
        Horizontal reflectivity [mm^6/m^3]

    Notes
    -----
    Ensure that both powers have the same units!

    Alternating horizontally and linearly polarized pulses are averaged.
    """
    zh = np.atleast_1d(z_h)
    zv = np.atleast_1d(z_v)
    if len(zh) != len(zv):
        raise ValueError('Input variables must be same length')
        return

    zdp = np.full_like(zh, np.nan)
    good = np.where(zh > zv)
    zdp[good] = 10.* np.log10(zh[good] - zv[good])
    return zdp


def hdr(dbz_h, zdr):
    """
    Differential reflectivity [dB] hail signature.

    From Aydin et al. (1986), Eqns 4-5

    Parameters
    ----------
    dbz_h : float or array
        Horizontal reflectivity [dBZ]
    zdr : float or array
        Differential reflectivity [dBZ]

    Notes
    -----
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
    zdr = np.atleast_1d(zdr)
    # Set the f(zdr) based upon observations
    f = np.full_like(zdr, np.nan)
    negind = np.where(zdr <= 0)
    lowind = np.where((zdr > 0) & (zdr <= 1.74))
    highind = np.where(zdr > 1.74)
    f[negind] = 27.
    f[lowind] = 19. * zdr[lowind] + 27.
    f[highind] = 60.
    # Calculate HDR
    return np.asarray(dbz_h) - f
