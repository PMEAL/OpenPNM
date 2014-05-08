r"""
===============================================================================
Submodule -- throat_area
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
    network.set_data(prop=propname,throats=geometry.throats(),data=value)

def cylinder(geometry,
             network,
             propname,
             diameter='diameter',
             **params):
    r"""
    Calculate throat area for a cylindrical throat
    """
    D = network.get_data(prop=diameter,throats=geometry.throats())
    value = sp.constants.pi/4*(D)**2
    network.set_data(prop=propname,throats=geometry.throats(),data=value)

def cuboid(geometry,
           network,
           propname,
           diameter='diameter',
           **params):
    r"""
    Calculate throat area for a cuboid throat
    """
    D = network.get_data(prop=diameter,throats=geometry.throats())
    value = (D)**2
    network.set_data(prop=propname,throats=geometry.throats(),data=value)
    
