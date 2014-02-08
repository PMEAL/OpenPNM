
"""
module throat_diameter
===============================================================================

"""
import scipy as sp
import scipy.stats as spst
import os
propname = os.path.splitext(os.path.basename(__file__))[0]
propname = propname.split('_')[1]

def constant(geometry,network,value,**params):
    r"""
    Assigns specified constant value
    """
    network.set_throat_data(subdomain=geometry,prop=propname,data=value)

def cylinder(geometry,network,**params):
    r"""
    Calculate throat diameter from seeds for a cylindrical throat
    """
    prob_fn = getattr(spst,params['name'])
    P = prob_fn(params['shape'],loc=params['loc'],scale=params['scale'])
    value=P.ppf(network.get_throat_data(prop='seed'))
    network.set_throat_data(subdomain=geometry,prop=propname,data=value)

def cuboid(geometry,network,**params):
    r"""
    Calculate throat diameter from seeds for a cuboidal throat
    """
    print('cuboid: nothing yet')
    
def voronoi(geometry,network,**params):
    r"""
    Calculate throat diameter from analysis of Voronoi facets
    """
    print('voronoi: nothing yet')