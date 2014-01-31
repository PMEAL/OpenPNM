
"""
module pore_volume
===============================================================================

"""
import scipy as sp

def constant(geometry,network,value,**params):
    r"""
    Assigns specified constant value
    """
    network.set_pore_data(prop='volume',data=value)

def sphere(geometry,network,**params):
    r"""
    Calculate pore volume from diameter for a spherical pore body
    """
    network.set_pore_data(prop='volume',data=sp.pi/6*network.get_pore_data(prop='diameter')**3)
    
def cube(geometry,network,**params):
    r"""
    Calculate pore volume from diameter for a cubic pore body
    """
    network.set_pore_data(prop='volume',data=network.get_pore_data(prop='diameter')**3)