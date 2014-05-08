r"""
===============================================================================
Submodule -- pore_diameter
===============================================================================

"""
import scipy as sp
import scipy.stats as spst


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

    