r"""
===============================================================================
Submodule -- throat_volume
===============================================================================

"""
import scipy as _sp

def cylinder(network,
             throats,
             throat_length='throat.length',
             throat_diameter='throat.diameter',
             **kwargs):
    r"""
    Calculate throat diameter from seeds for a cylindrical throat
    - note: this will need to account for volume taken up by spherical pore bodies
    """
    leng = network[throat_length][throats]
    diam = network[throat_diameter][throats]
    value = _sp.pi/4*leng*diam**2
    return value

def cuboid(network,
           throats,
           throat_length='throat.length',
           throat_diameter='throat.diameter',
           **kwargs):
    r"""
    Calculate throat volume of cuboidal throat
    - note: this will need to account for volume taken up by spherical pore bodies
    """
    leng = network[throat_length][throats]
    diam = network[throat_diameter][throats]
    value = leng*diam**2
    return value

def extrusion(network,
              throats,
              throat_length='throat.length',
              throat_area='throat.area',
              **kwargs):
    r"""
    Calculate volume from the throat area and the throat length
    """
    leng = network[throat_length][throats]
    area = network[throat_area][throats]
    value = leng*area
    return value
    
def voronoi(network,
            throats,
            throat_length='throat.length',
            throat_area='throat.area',
            **kwargs):
    r"""
    Calculate volume from the throat area and the throat length
    """
    leng = network[throat_length][throats]
    area = network[throat_area][throats]
    value = leng*area
    return value