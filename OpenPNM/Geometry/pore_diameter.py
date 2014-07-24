r"""
===============================================================================
Submodule -- pore_diameter
===============================================================================

"""
import scipy as sp
import scipy.stats as spst
from scipy.special import cbrt

def constant(geometry,
             network,
             propname,
             value,
             **params):
    r"""
    Assign specified constant value
    """
    network.set_data(prop=propname,pores=geometry.pores(),data=value)

def sphere(geometry,
           network,
           propname,
           seed='seed',
           **params):
    r"""
    Calculate pore diameter from seed values for a spherical pore body
    """
    prob_fn = getattr(spst,params['name'])
    P = prob_fn(params['shape'],loc=params['loc'],scale=params['scale'])
    value = P.ppf(network.get_data(prop=seed,pores=geometry.pores()))
    network.set_data(prop=propname,pores=geometry.pores(),data=value)

def voronoi(geometry,
            network,
            propname,
            **params):
    r"""
    Calculate pore diameter from equivalent sphere - volumes must be calculated first
    """
    pore_vols = network.get_data(prop='volume',pores=geometry.pores())
    value = cbrt(6*pore_vols/sp.pi)
    network.set_data(prop=propname,pores=geometry.pores(),data=value)