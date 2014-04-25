
"""
module throat_volume
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
    network.set_throat_data(locations=geometry.get_throat_locations(),prop=propname,data=value)

def cylinder(geometry,
             network,
             propname,
             throat_length='length',
             throat_diameter='diameter',
             **params):
    r"""
    Calculate throat diameter from seeds for a cylindrical throat
    - note: this will need to account for volume taken up by spherical pore bodies
    """
    value=sp.pi/4*network.get_throat_data(prop=throat_length,locations=geometry)*network.get_throat_data(prop=throat_diameter,locations=geometry)**2
    network.set_throat_data(locations=geometry.get_throat_locations(),prop=propname,data=value)

def cuboid(geometry,
           network,
           propname,
           throat_length='length',
           throat_diameter='diameter',
           **params):
    r"""
    Calculate throat volume of cuboidal throat
    - note: this will need to account for volume taken up by spherical pore bodies
    """
    value=network.get_throat_data(prop=throat_length,locations=geometry)*network.get_throat_data(prop=throat_diameter,locations=geometry)**2
    network.set_throat_data(locations=geometry.get_throat_locations(),prop=propname,data=value)