
"""
module pore_volume
===============================================================================

"""
import scipy as sp

def constant(geometry,network,value,**params):
    r"""
    Assigns specified constant value
    """
    network.pore_properties['volume'] = value

def sphere(geometry,network,**params):
    r"""
    Calculate pore volume from diameter for a spherical pore body
    """
    network.pore_properties['volume'] = sp.pi/6*network.pore_properties['diameter']**3
    
def cube(geometry,network,**params):
    r"""
    Calculate pore volume from diameter for a cubic pore body
    """
    network.pore_properties['volume'] = network.pore_properties['diameter']**3