r"""
===============================================================================
pore_volume --
===============================================================================

"""
import scipy as _sp


def sphere(target, pore_diameter='pore.diameter'):
    r"""
    Calculate pore volume from diameter assuming a spherical pore body

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    pore_diameter : string
        The dictionary key to the pore diameter values

    """
    diams = target[pore_diameter]
    value = _sp.pi/6*diams**3
    return value


def cube(target, pore_diameter='pore.diameter'):
    r"""
    Calculate pore volume from diameter assuming a cubic pore body

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    pore_diameter : string
        The dictionary key to the pore diameter values

    """
    diams = target[pore_diameter]
    value = diams**3
    return value
