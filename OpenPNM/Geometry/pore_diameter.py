
"""
module pore_diameter
===============================================================================

"""
import scipy as sp
import scipy.stats as spst

def constant(geometry,network,value,**params):
    r"""
    Assign specified constant value
    """
    network.set_pore_data(prop='diameter',data=value)

def sphere(geometry,network,**params):
    r"""
    Calculate pore diameter from seed values for a spherical pore body
    """
    try:
        prob_fn = getattr(spst,params['name'])
        P = prob_fn(params['shape'],loc=params['loc'],scale=params['scale'])
        network.set_pore_data(prop='diameter',data=P.ppf(network.get_pore_data(prop='seed')))
    except: 
        print('given network does not contain seed values')

    