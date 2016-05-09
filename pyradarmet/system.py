# -*- coding: utf-8 -*-
"""
pyradarmet.system
=================

Functions to calculate radar system characteristics.

References
----------
Rinehart (1997), Radar for Meteorologists.
"""
import numpy as np


speed_of_light = 3e8 # Speed of light [m/s]
kBoltz = 1.381e-23 # Boltzmann's constant [ m^2 kg s^-2 K^-1]

def gain_Pratio(p1, p2):
    """
    Antenna gain [dB] via power ratio.

    From Rinehart (1997), Eqn 2.1

    Parameters
    ----------
    p1 : float or array
        Power on the beam axis [W]
    p2 : float or array
        Power from an isotropic antenna [W]

    Notes
    -----
    Ensure that both powers have the same units!
    If arrays are used, either for one okay, or both must be the same length.
    """
    return 10. * np.log10(np.asarray(p1) / np.asarray(p2))


def freq(lam):
    """
    Frequency [Hz] given wavelength.

    Parameters
    ----------
    lam : float or array
        Wavelength [m]
    """
    return speed_of_light / np.asarray(lam)


def wavelength(freq):
    """
    Wavelength [m] given frequency.

    Parameters
    ----------
    freq : float or array
        Frequency [Hz]
    """
    return speed_of_light / np.asarray(freq)


def pulse_length(pdur):
    """
    Pulse length [m] given pulse duration.

    Parameters
    ----------
    pDur : float or array
        Pulse duration [s]

    Notes
    -----
    This equation is only interested in return pulses to radar, leading
    to the factor of 1/2.
    """
    return speed_of_light * np.asarray(pdur) / 2.


def pulse_duration(tau):
    """
    Pulse duration [s] from pulse length.

    Parameters
    ----------
    tau : float or array
        Pulse length [m]
    """
    return 2 * np.asarray(tau) / speed_of_light


def radar_const(power_t, gain, tau, lam, bw_h, bw_v, aloss, rloss):
    """
    Radar constant. Unitless.

    From CSU Radar Meteorology notes, AT 741

    Parameters
    ----------
    power_t : float
        Transmitted power [W]
    gain : float
        Antenna Gain [dB]
    tau : float
        Pulse Width [s]
    lam : float
        Radar wavelength [m]
    bw_h : float
        Horizontalntenna beamwidth [degrees]
    bw_v : float
        Vertical antenna beamwidth [degrees]
    aloss : float
        Antenna/waveguide/coupler loss [dB]
    rloss : float
        Receiver loss [dB]
    """
    # Convert from dB to linear units
    alosslin = 10**(aloss / 10.)
    rlosslin = 10**(rloss / 10.)
    gainlin = 10**(gain / 10.)

    # Convert beamwidth to radians
    bw_hr = np.deg2rad(bw_h)
    bw_vr = np.deg2rad(bw_v)

    # Calculate the numerator
    Numer = (np.pi**3 * speed_of_light * power_t * gainlin**2 * tau *
             bw_hr * bw_vr * alosslin * rlosslin)

    # Calculate the denominator
    Denom = 1024. * np.log(2) * lam**2
    return Numer/Denom


def ant_eff_area(gain, lam):
    """
    Antenna effective area. [m^-2]

    From Rinehart (1997), Eqn 4.5

    Parameters
    ----------
    gain : float or array
        Antenna Gain [dB]
    lam : float
        Radar wavelength [m]
    """
    # Convert from dB to linear units
    gainlin = 10**(np.asarray(gain) / 10.)
    return gainlin * lam**2 / (4 * np.pi)


def power_target(power_t, gain, areat, r):
    """
    Power [W] intercepted by target.

    From Rinehart (1997), Eqn 4.3

    Parameters
    ----------
    power_t : float
        Transmitted power [W]
    gain : float
        Antenna gain [dB]
    areat : float
        Area of target [m^2]
    r : float or array
        Distance to sample volume from radar [m]
    """
    # Convert from dB to linear units
    gainlin = 10**(gain / 10.)
    return (power_t * gainlin * areat) / (4 * np.pi * np.asarray(r)**2)


def xsec_bscatter_sphere(diam, lam, dielectric=0.93):
    """
    Backscatter cross-sectional area [m^-2] of a sphere using the Rayleigh approximation.

    From Rinehart (1997), Eqn 4.9 and 5.7

    Parameters
    ----------
    diam : float or array
        Diamter of target [m]
    lam : float
        Radar wavelength [m]
    dielectric : float
        Dielectric factor [unitless]

    Notes
    -----
    The Rayleigh approximation is good when the diamter of a spherical particle
    is much smaller than the wavelength of the radar (D/wavelength= 1/16).  This
    condition leads to the relationship that the area is proportional to the
    sixth power of the diameter.

    The default is for a dielectric factor value for water.  This can be
    changed by the user, e.g. K=0.208 for particle sizes of equivalent melted
    diameters or K=0.176 for particle sizes of equivalent ice spheres.
    """
    return (np.pi**5 * dielectric**2 * np.asarray(diam)**6) / lam**4


def norm_xsec_bscatter_sphere(diam, lam, dielectric=0.93):
    """
    Normalized Backscatter cross-sectional area [m^2] of a sphere using the Rayleigh approximation.

    From Rinehart (1997), Eqn 4.9 and 5.7 and Battan Ch. 4.5

    Parameters
    ----------
    diam : float or array
        Diamter of targer [m]
    lam : float
        Radar wavelength [m]
    dielectric : float
        Dielectric factor [unitless]

    Notes
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
    sig = xsec_bscatter_sphere(np.asarray(diam), lam, dielectric)
    return sig/ (np.pi * (np.asarray(diam)/2.)**2)


def size_param(diam, lam):
    """
    Size parameter calculation. Unitless.

    From Rinehart (1997), Eqn 4.9 and 5.7 and Battan Ch. 4.5

    Parameters
    ----------
    diam : float or float
        Diamter of target [m]
    lam : float
        Radar wavelength [m]

    Notes
    -----
    The size paramter can be used along with the backscattering cross-section to
    distinguish ice and water dielectric characteristics.  For example:
    Alpha < 2 the backscattering cross-section of ice is smaller than water,
    Alpha > 2 the opposite is true due to the fact that absorption in water
    exceeds that in ice.
    """
    return 2 * np.pi * np.asarray(diam)/2. / lam


def power_return_target(power_t, gain, lam, sig, r):
    """
    Power [W] returned y target located at the center of the antenna beam pattern.

    From Rinehart (1997), Eqn 4.7

    Parameters
    ----------
    power_t : float
        Transmitted power [W]
    gain : float
        Antenna gain [dB]
    lam : float
        Radar wavelength [m]
    sig : float
        Backscattering cross-sectional area of target [m^2]
    r : float or array
        Distance to sample volume from radar [m]
    """
    # Convert from dB to linear units
    gainlin = 10**(gain/10.)
    return ((power_t * gainlin**2 * lam**2 * sig) /
            (64 * np.pi**3 * np.asarray(r)**4))


def thermal_noise(bandwidth, Units, noise_temp=290.):
    """
    Thermal noise power [W or 'dBm'].

    From CSU Radar Meteorology notes, AT741

    Parameters
    ----------
    bandwidth : float or array
        Receiver bandwidth [Hz]
    Units : float
        String of nits desired, can be 'W' or 'dBm'
    Ts : float
        Reciever noise temperature [K]

    Notes
    -----
    Reciever noise temp set to conventional 290K by default
    """
    # Calculate the noise, convert if requested
    noise = kBoltz * noise_temp * np.asarray(bandwidth)

    if Units.upper()=='W':
        noiset = noise
    elif Units.upper()=='DBM':
        noiset = 10. * np.log10(noise/10**-3)
    else:
        print("Units must be in 'W' or 'dBm'")
        noiset = np.nan
    return noiset
