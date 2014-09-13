r"""
===============================================================================
Submodule -- mixture_props
===============================================================================

Models for adopting properties from the mixture into each component phase

"""
import scipy as sp

def temperature(phase,mixture,**kwargs):
    r"""
    Adopts the temperature of the parent mixture phase
    """
    vals = mixture['pore.temperature']
    return vals

def pressure(phase,mixture,**kwargs):
    r"""
    Adopts the pressure of the parent mixture phase
    """
    vals = mixture['pore.pressure']
    return vals
