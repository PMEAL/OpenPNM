
"""
module pore_diameter
===============================================================================

"""
import scipy as sp
import scipy.stats as spst
import os
propname = os.path.splitext(os.path.basename(__file__))[0]
propname = propname.split('_')[1]

def constant(geometry,network,value,**params):
    r"""
    Assign specified constant value
    """
    network.set_pore_data(labels=geometry,prop=propname,data=value)

def sphere(geometry,network,**params):
    r"""
    Calculate pore diameter from seed values for a spherical pore body
    """
    prob_fn = getattr(spst,params['name'])
    P = prob_fn(params['shape'],loc=params['loc'],scale=params['scale'])
    value = P.ppf(network.get_pore_data(prop='seed'))
    network.set_pore_data(labels=geometry,prop=propname,data=value)

    