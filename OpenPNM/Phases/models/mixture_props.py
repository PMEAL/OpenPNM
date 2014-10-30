r"""
===============================================================================
Submodule -- mixture_props
===============================================================================

Models for adopting properties from the mixture into each component phase

"""
import scipy as sp

def temperature(phase,**kwargs):
    r"""
    Adopts the temperature of the parent mixture phase
    """
    mixture = phase._phases[0]
    vals = mixture['pore.temperature']
    return vals

def pressure(phase,**kwargs):
    r"""
    Adopts the pressure of the parent mixture phase
    """
    mixture = phase._phases[0]
    vals = mixture['pore.pressure']
    return vals
