# -*- coding: utf-8 -*-
r"""
===============================================================================
Submodule -- vapor_pressure
===============================================================================

Methods for predicing the vapor pressure of pure species

"""
import scipy as sp

def antoine(phase,A,B,C,**kwargs):
    r"""
    Uses Antoine equation [1]_ to estimate vapor pressure of a pure component

    Parameters
    ----------
    A, B, C :  float, array_like
            Antoine vapor pressure constants for pure compounds.  Note that
            these constants are traditionally reported such that they give
            vapor pressure in mmHg. This function converts pressure to Pascals.

    [1] Antoine, C. (1888), Vapor Pressure: a new relationship between pressure 
        and temperature, Comptes Rendus des Séances de l'Académie des Sciences 
        (in French) 107: 681–684, 778–780, 836–837
    
    """
    T = phase['pore.temperature']
    value = (10**(A-B/(C+T-273.15)))*133.3
    return value
