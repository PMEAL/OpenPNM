
"""
module ThroatVolume
===============================================================================

"""
import scipy as sp

def constant(network, value,**params):
    r"""
    Assigns specified constant value
    """
    network.throat_properties['volume'] = value

def cylinder(network,**params):
    r"""
    Calculate throat diameter from seeds for a cylindrical throat
    - note: this will need to account for volume taken up by spherical pore bodies
    """
    network.throat_properties['volume'] = sp.pi/4*network.throat_properties['length']*network.throat_properties['diameter']**2

def cuboid(network, **params):
    r"""
    Calculate throat volume of cuboidal throat
    - note: this will need to account for volume taken up by spherical pore bodies
    """
    network.throat_properties['volume'] = network.throat_properties['length']*network.throat_properties['diameter']**2