
"""
module pore_volume
===============================================================================

"""
import scipy as sp
import os
propname = os.path.splitext(os.path.basename(__file__))[0]
propname = propname.split('_')[1]

def constant(geometry,network,value,**params):
    r"""
    Assigns specified constant value
    """
    network.set_pore_data(subdomain=geometry,prop=propname,data=value)

def sphere(geometry,network,**params):
    r"""
    Calculate pore volume from diameter for a spherical pore body
    """
    value=sp.pi/6*network.get_pore_data(prop='diameter')**3
    network.set_pore_data(subdomain=geometry,prop=propname,data=value)
    
def cube(geometry,network,**params):
    r"""
    Calculate pore volume from diameter for a cubic pore body
    """
    value=network.get_pore_data(prop='diameter')**3
    network.set_pore_data(subdomain=geometry,prop=propname,data=value)