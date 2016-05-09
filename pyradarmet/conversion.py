# -*- coding: utf-8 -*-
"""
pyradarmet.conversion
=====================

Functions for converting common radar units.
"""
import numpy as np


def dBZ2Z(dBZ):
    """
    Convert from log [dBZ] to linear Z [mm^6 m^−3] units.

    Parameters
    ----------
    dBZ : float or array
        logarithmic reflectivity value
    """
    return 10.**(np.asarray(dBZ)/10.)


def Z2dBZ(Zlin):
    """
    Convert from linear Z [mm^6 m^−3] to log [dBZ] units.

    Parameters
    ----------
    Zlin : float or array
        linear reflectivity units
    """
    return 10. * np.log10(np.asarray(Zlin))


def si2kmh(SI):
    """
    Convert from SI [m/s] wind units to km/h.

    Parameters
    ----------
    SI : float or array
        Wind in SI units (m/s)
    """
    return np.asarray(SI) * 3600. / 1000.


def si2mph(SI):
    """
    Convert from SI wind units to miles/h [mph].

    Parameters
    ----------
    SI: float or array
        Wind in SI units (m/s)
    """
    return np.asarray(SI) * 0.62137 / 1000. * 3600.


def si2kts(SI):
    """
    Convert from SI wind units to knots [kt].

    Parameters
    ----------
    SI: float or array
        Wind in SI units (m/s)
    """
    return np.asarray(SI) * 0.51


def kmh2si(kmh):
    """
    Convert from km/h to SI wind units [m/s].

    Parameters
    ----------
    kmh: float or array
        Wind in km/hr
    """
    return np.asarray(kmh) * 1000. / 3600.


def mph2si(mph):
    """
    Convert from miles/h to SI wind units [m/s].

    Parameters
    ----------
    mph: float or array
        Wind in miles per hour
    """
    return np.asarray(mph) * 1000. / (0.62137 * 3600.)


def kts2si(kts):
    """
    Convert from knots to SI wind units [m/s].

    Parameters
    ----------
    kts: float or array
        Wind in knots
    """
    return np.asarray(kts) / 0.51
