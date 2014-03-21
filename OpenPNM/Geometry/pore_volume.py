r"""
===============================================================================
Submodule -- pore_volume
===============================================================================

"""
import scipy as sp

def constant(geometry,
             network,
             propname,
             value,
             **params):
    r"""
    Assigns specified constant value
    """
    network.set_pore_data(locations=geometry,prop=propname,data=value)

def sphere(geometry,
           network,
           propname,
           pore_diameter='diameter',
           **params):
    r"""
    Calculate pore volume from diameter for a spherical pore body
    """
    value=sp.pi/6*network.get_pore_data(prop=pore_diameter,locations=geometry)**3
    network.set_pore_data(locations=geometry,prop=propname,data=value)
    
def cube(geometry,
         network,
         propname,
         pore_diameter='diameter',
         **params):
    r"""
    Calculate pore volume from diameter for a cubic pore body
    """
    value=network.get_pore_data(prop=pore_diameter,locations=geometry)**3
    network.set_pore_data(locations=geometry,prop=propname,data=value)