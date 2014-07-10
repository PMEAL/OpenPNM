r"""
===============================================================================
Submodule -- throat_volume
===============================================================================

"""
import scipy as _sp

def cylinder(network,throats,**kwargs):
    r"""
    Calculate throat diameter from seeds for a cylindrical throat
    - note: this will need to account for volume taken up by spherical pore bodies
    """
    leng = network['throat.length'][throats]
    diam = network['throat.diameter'][throats]
    value = _sp.pi/4*leng*diam**2
    return value

def cuboid(network,throats,**kwargs):
    r"""
    Calculate throat volume of cuboidal throat
    - note: this will need to account for volume taken up by spherical pore bodies
    """
    leng = network['throat.length'][throats]
    diam = network['throat.diameter'][throats]
    value = leng*diam**2
    return value

def extruded(network,throats,**kwargs):
    r"""
    Calculate volume from the throat area and the throat length
    """
    leng = network['throat.length'][throats]
    area = network['throat.area'][throats]
    value = leng*area
    return value