
"""
module PoreDiameter
===============================================================================

"""
import scipy as sp
import scipy.stats as spst

def constant(network,value,**params):
    r"""
    Assign specified constant value
    """
    network.pore_properties['diameter'] = value

def sphere(network,**params):
    r"""
    Calculate pore diameter from seed values for a spherical pore body
    """
    prob_fn = getattr(spst,params['name'])
    P = prob_fn(params['shape'],loc=params['loc'],scale=params['scale'])
    network.pore_properties['diameter'] = P.ppf(network.pore_properties['seed'])
    