r"""
===============================================================================
Submodule -- throat_diameter
===============================================================================

"""
import scipy as sp
import scipy.stats as spst
import numpy as np

def constant(geometry,
             network,
             propname,
             value,
             **params):
    r"""
    Assigns specified constant value
    """
    network.set_data(prop=propname,throats=geometry.throats(),data=value)

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
    value=P.ppf(network.get_data(prop=seed,throats=geometry.throats()))
    network.set_data(prop=propname,throats=geometry.throats(),data=value)

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
    Equivalent circular diameter from voronoi area
    Could do better here and work out minimum diameter from verts
    """
    areas = network.get_throat_data(prop='area')   
    value = 2*np.sqrt(areas/np.pi)#64 bit sqrt doesn't work!
    network.set_data(prop=propname,throats=geometry.throats(),data=value)
    #print('voronoi: nothing yet')
