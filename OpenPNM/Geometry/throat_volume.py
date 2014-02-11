
"""
module throat_volume
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
    network.set_throat_data(labels=geometry,prop=propname,data=value)

def cylinder(geometry,network,**params):
    r"""
    Calculate throat diameter from seeds for a cylindrical throat
    - note: this will need to account for volume taken up by spherical pore bodies
    """
    value=sp.pi/4*network.get_throat_data(prop='length')*network.get_throat_data(prop='diameter')**2
    network.set_throat_data(labels=geometry,prop=propname,data=value)

def cuboid(geometry,network,**params):
    r"""
    Calculate throat volume of cuboidal throat
    - note: this will need to account for volume taken up by spherical pore bodies
    """
    value=network.get_throat_data(prop='length')*network.get_throat_data(prop='diameter')**2
    network.set_throat_data(labels=geometry,prop=propname,data=value)