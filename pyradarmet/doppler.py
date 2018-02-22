# -*- coding: utf-8 -*-
"""
pyradarmet.doppler
==================

Functions to calculate radar characteristics for Doppler radar.

References
----------
Rinehart (1997), Radar for Meteorologists.
Jorgensen (1983; JCAM), Feasibility Test of an Airborne Pulse-Doppler
Meteorological Radar.
"""
import numpy as np


speed_of_light = 3e8  # Speed of light [m/s]


def freq(lam):
    """Frequency [Hz] given wavelength.

    Parameters
    ----------
    lam : float or array
        Wavelength [m]
    """
    return speed_of_light / np.asarray(lam)


def wavelength(freq):
    """Wavelength [m] given frequency.

    Parameters
    ----------
    freq : float or array
        Frequency [Hz]
    """
    return speed_of_light / np.asarray(freq)


def pulse_duration(tau):
    """Pulse duration [s] from pulse length.

    Parameters
    ----------
    tau : float or array
        Pulse length [m]
    """
    return 2 * np.asarray(tau) / speed_of_light


def pulse_length(pdur):
    """Pulse length [m] from pulse duration.

    Parameters
    ----------
    pDur : float or array
        Pulse duration [s]
    """
    return speed_of_light * np.asarray(pdur) / 2


def fmax(prf):
    """Maximum frequency [Hz] given PRF.

    From Rinehart (1997), Eqn 6.8

    Parameters
    ----------
    PRF : float or array
        Pulse repetition frequency [Hz]
    """
    return np.asarray(prf) / 2.


def Vmax(PRF, lam):
    """Nyquist velocity, or maximum unambiguous Doppler velocity (+ or -) [m/s].

    From Rinehart (1997), Eqn 6.7

    Parameters
    ----------
    PRF : float or array
        Radar pulse repetition frequency [Hz]
    lam : float or array
        Radar wavelength [m]
    """
    return np.asarray(PRF) * lam / 4.


def Rmax(PRF):
    """Maximum unamiguous range [m].

    From Rinehart (1997), Eqn 6.11

    Parameters
    ----------
    PRF : float or array
        Pulse repetition frequency [Hz]
    """
    return speed_of_light / (2. * np.asarray(PRF))


def doppler_dilemma(varin, lam):
    """
    The "Doppler dilemma" is the fact that both the Nyquist velocity and
    unambiguous maximum range of the radar are based upon
    the PRF of the system.

    However, they are inversely proportional, meaning that increasing one
    requires a decrease in the other.  A trade-off inherent in Doppler radar
    systems.  This relationship allows a solution for one variable given the
    value of the other.

    From Rinehart (1997), Eqn 6.12

    Parameters
    ----------
    varin : float or array
        Nyquist Velocity [m/s] or Maximum unambiguous range [m]
    lam : float
        Radar wavelength [m]
    """
    return (speed_of_light * lam / 8.) / np.asarray(varin)

######################
#  MOBILE PLATFORMS  #
######################


def Vshift(ground_speed, psi):
    """
    Adjusted Doppler velocity [m/s] from a mobile platform.
    Shift in Doppler velocity from mobile perspective.

    Jorgensen (1983), Eqn 2

    Parameters
    ----------
    ground_speed : float or array
        Gound speed [m/s]
    psi : float or array
        Angle between actual azimuth and fore/aft angle [deg]

    Notes
    -----
    In the case of a mobile platform (such as the NOAA P-3 aircraft, the
      Doppler velocity must be adjusted for movement of the scanning platform.

    The fore/aft angle is defined as the angle fore or aft from a plane
      normal to the direction of motion
    """
    len_gs = len(np.asarray(ground_speed))
    len_psi = len(np.asarray(psi))
#    if  len_gs != len_psi and len_gs != 1 and len_psi != 1:
    return np.asarray(ground_speed) * np.cos(np.deg2rad(psi))


def Vmax_dual(lam, prf1, prf2):
    """Doppler velocity [m/s] from dual PRF scheme radar (+ or -).

    From Jorgensen (1983), Eqn 2

    Parameters
    ----------
    lam : float
        Radar wavelength [m]
    prf1 : float
        First Pulse repetition frequency [Hz]
    prf2 : float
        Second Pulse repetition frequency [Hz]

    Notes
    -----
    In the case of a mobile platform (such as the NOAA P-3 aircraft, the
    Doppler velocity must be adjusted for movement of the scanning platform.

    The fore/aft angle is defined as the angle fore or aft from a plane
    normal to the direction of motion
    """

    Vmax = lam / (4 * ((1. / prf1) - (1. / prf2)))

    return Vmax
