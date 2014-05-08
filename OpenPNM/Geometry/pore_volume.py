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
    network.set_data(prop=propname,pores=geometry.pores(),data=value)

def sphere(geometry,
           network,
           propname,
           diameter='diameter',
           **params):
    r"""
    Calculate pore volume from diameter for a spherical pore body
    """
    value=sp.pi/6*network.get_data(prop=diameter,pores=geometry.pores())**3
    network.set_data(prop=propname,pores=geometry.pores(),data=value)
    
def cube(geometry,
         network,
         propname,
         diameter='diameter',
         **params):
    r"""
    Calculate pore volume from diameter for a cubic pore body
    """
    value=network.get_data(prop=diameter,pores=geometry.pores())**3
    network.set_data(prop=propname,pores=geometry.pores(),data=value)
