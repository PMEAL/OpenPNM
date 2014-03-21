r"""
===============================================================================
Submodule -- throat_diameter
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
    Assigns specified constant value
    """
    network.set_throat_data(locations=geometry,prop=propname,data=value)

def cylinder(geometry,
             network,
             propname,
             seed='seed',
             **params):
    r"""
    Calculate throat diameter from seeds for a cylindrical throat
    """
    prob_fn = getattr(spst,params['name'])
    P = prob_fn(params['shape'],loc=params['loc'],scale=params['scale'])
    value=P.ppf(network.get_throat_data(prop=seed,locations=geometry))
    network.set_throat_data(locations=geometry,prop=propname,data=value)

def cuboid(geometry,
           network,
           propname,
           **params):
    r"""
    Calculate throat diameter from seeds for a cuboidal throat
    """
    print('cuboid: nothing yet')
    
def voronoi(geometry,
            network,
            propname,
            **params):
    r"""
    Calculate throat diameter from analysis of Voronoi facets
    """
    print('voronoi: nothing yet')