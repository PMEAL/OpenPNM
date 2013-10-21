
"""
module Diffusivity
===============================================================================

"""
import OpenPNM
import scipy as sp

def constant(DAB,**params):
    return DAB

def na(**params):
    return 'n/a'

def Fuller(network,MA=31.99,MB=28.01,vA=16.6,vB=17.9,**params):
    r"""
    Uses the Fuller model estimate diffusion coefficient at conditions of interest
    """
    T = network.pore_conditions['temperature']
    P = network.pore_conditions['pressure']
    MAB = 2*(1/MA+1/MB)**(-1)
    DAB = 0.00143*T**1.75/(P*MAB**0.5*(vA**(1/3)+vB**(1/3))**2)*1e-4
    return DAB

