
"""
module throat_volume
===============================================================================

"""
import scipy as sp

def constant(geometry,network,value,**params):
    r"""
    Assigns specified constant value
    """
    network.throat_properties['volume'] = value

def cylinder(geometry,network,**params):
    r"""
    Calculate throat diameter from seeds for a cylindrical throat
    - note: this will need to account for volume taken up by spherical pore bodies
    """
    try:
        network.throat_properties['volume'] = sp.pi/4*network.throat_properties['length']*network.throat_properties['diameter']**2
    except:
        print('Cannot calculate volume, some required information is missing')

def cuboid(geometry,network,**params):
    r"""
    Calculate throat volume of cuboidal throat
    - note: this will need to account for volume taken up by spherical pore bodies
    """
    network.throat_properties['volume'] = network.throat_properties['length']*network.throat_properties['diameter']**2