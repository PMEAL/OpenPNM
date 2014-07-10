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
    network.set_data(prop=propname,throats=geometry.throats(),data=value)

def cylinder(geometry,
             network,
             propname,
             diameter='diameter',
             length='length',
             **params):
    r"""
    Calculate throat area for a cylindrical throat
    """
    D = network.get_data(prop=diameter,throats=geometry.throats())
    L = network.get_data(prop=length,throats=geometry.throats())
    value = sp.constants.pi/(D)*L
    network.set_data(prop=propname,throats=geometry.throats(),data=value)

def cuboid(geometry,
           network,
           propname,
           diameter='diameter',
           length='length',
           **params):
    r"""
    Calculate throat area for a cuboid throat
    """
    D = network.get_data(prop=diameter,throats=geometry.throats())
    L = network.get_data(prop=length,throats=geometry.throats())
    value = 4*D*L
    network.set_data(prop=propname,throats=geometry.throats(),data=value)
    
def voronoi(geometry,
            network,
            propname,
            perimeter='perimeter',
            length='length',
            **params):
    r"""
    Calculate surface area from perimeter and lenght - 
    perimeter calculated when throat area is calculated so must be run in correct order
    """
    P = network.get_data(prop=perimeter,throats=geometry.throats())
    L = network.get_data(prop=length,throats=geometry.throats())
    value = P*L
    network.set_data(prop=propname,throats=geometry.throats(),data=value)