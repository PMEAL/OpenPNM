
"""
module throat_area
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
    network.set_throat_data(locations=geometry.get_throat_locations(),prop=propname,data=value)

def cylinder(geometry,
             network,
             propname,
             diameter='diameter',
             **params):
    r"""
    Calculate throat area for a cylindrical throat
    """
    D = network.get_throat_data(prop=diameter,locations=geometry)
    value = sp.constants.pi/4*(D)**2
    network.set_throat_data(locations=geometry.get_throat_locations(),prop=propname,data=value)

def cuboid(geometry,
           network,
           propname,
           diameter='diameter',
           **params):
    r"""
    Calculate throat area for a cuboid throat
    """
    D = network.get_throat_data(prop=diameter,locations=geometry)
    value = (D)**2
    network.set_throat_data(locations=geometry.get_throat_locations(),prop=propname,data=value)
    
