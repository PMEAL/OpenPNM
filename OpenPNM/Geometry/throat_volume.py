
"""
module throat_volume
===============================================================================

"""
import scipy as sp

def constant(geometry,network,value,**params):
    r"""
    Assigns specified constant value
    """
    network.set_throat_data(prop='volume',data=value)

def cylinder(geometry,network,**params):
    r"""
    Calculate throat diameter from seeds for a cylindrical throat
    - note: this will need to account for volume taken up by spherical pore bodies
    """
    try:
        network.set_throat_data(prop='volume',data=sp.pi/4*network.get_throat_data(prop='length')*network.get_throat_data(prop='diameter')**2)
    except:
        print('Cannot calculate volume, some required information is missing')

def cuboid(geometry,network,**params):
    r"""
    Calculate throat volume of cuboidal throat
    - note: this will need to account for volume taken up by spherical pore bodies
    """
    network.set_throat_data(prop='volume',data=network.get_throat_data(prop='length')*network.get_throat_data(prop='diameter')**2)