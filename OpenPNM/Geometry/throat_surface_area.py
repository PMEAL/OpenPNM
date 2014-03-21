r"""
===============================================================================
Submodule -- throat_surface_area
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
             diameter='diameter',
             length='length',
             **params):
    r"""
    Calculate throat area for a cylindrical throat
    """
    D = network.get_throat_data(prop=diameter,locations=geometry)
    L = network.get_throat_data(prop=length,locations=geometry)
    value = sp.constants.pi/(D)*L
    network.set_throat_data(locations=geometry,prop=propname,data=value)

def cuboid(geometry,
           network,
           propname,
           diameter='diameter',
           length='length',
           **params):
    r"""
    Calculate throat area for a cuboid throat
    """
    D = network.get_throat_data(prop=diameter,locations=geometry)
    L = network.get_throat_data(prop=length,locations=geometry)
    value = 4*D*L
    network.set_throat_data(locations=geometry,prop=propname,data=value)
    
