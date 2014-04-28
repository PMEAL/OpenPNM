
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
    network.set_data(prop=propname,throats=geometry.throats,data=value)

def cylinder(geometry,
             network,
             propname,
             diameter='diameter',
             length='length',
             **params):
    r"""
    Calculate throat area for a cylindrical throat
    """
    D = network.get_data(prop=diameter,throats=geometry.throats)
    L = network.get_data(prop=length,throats=geometry.throats)
    value = sp.constants.pi/(D)*L
    network.set_data(prop=propname,throats=geometry.throats,data=value)

def cuboid(geometry,
           network,
           propname,
           diameter='diameter',
           length='length',
           **params):
    r"""
    Calculate throat area for a cuboid throat
    """
    D = network.get_data(prop=diameter,throats=geometry.throats)
    L = network.get_data(prop=length,throats=geometry.throats)
    value = 4*D*L
    network.set_data(prop=propname,throats=geometry.throats,data=value)
    
