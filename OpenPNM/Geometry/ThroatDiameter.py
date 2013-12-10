
"""
module ThroatDiameter
===============================================================================

"""
import scipy as sp
import scipy.stats as spst

def constant(network, value,**params):
    r"""
    Assigns specified constant value
    """
    network.throat_properties['diameter'] = value

def cylinder(network,**params):
    r"""
    Calculate throat diameter from seeds for a cylindrical throat
    """
    prob_fn = getattr(spst,params['name'])
    P = prob_fn(params['shape'],loc=params['loc'],scale=params['scale'])
    network.throat_properties['diameter'] = P.ppf(network.throat_properties['seed'])

def cuboid(**params):
    r"""
    Calculate throat diameter from seeds for a cuboidal throat
    """
    print('cuboid: nothing yet')
    
def voronoi(**params):
    r"""
    Calculate throat diameter from analysis of Voronoi facets
    """
    print('voronoi: nothing yet')