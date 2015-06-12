r"""
===============================================================================
Submodule -- throat_volume
===============================================================================

"""
import scipy as _sp


def cylinder(geometry, throat_length='throat.length',
             throat_diameter='throat.diameter', **kwargs):
    r"""
    Calculate throat diameter from seeds for a cylindrical throat
    - note: this will need to account for volume taken up by spherical pore
    bodies
    """
    leng = geometry[throat_length]
    diam = geometry[throat_diameter]
    value = _sp.pi/4*leng*diam**2
    return value


def cuboid(geometry, throat_length='throat.length',
           throat_diameter='throat.diameter', **kwargs):
    r"""
    Calculate throat volume of cuboidal throat
    - note: this will need to account for volume taken up by spherical pore bodies
    """
    leng = geometry[throat_length]
    diam = geometry[throat_diameter]
    value = leng*diam**2
    return value


def extrusion(geometry, throat_length='throat.length',
              throat_area='throat.area', **kwargs):
    r"""
    Calculate volume from the throat area and the throat length
    """
    leng = geometry[throat_length]
    area = geometry[throat_area]
    value = leng*area
    return value
