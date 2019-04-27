r"""
These models calculate pore volumes depending on the specified shape

.. autofunction:: openpnm.models.geometry.pore_volume.sphere
.. autofunction:: openpnm.models.geometry.pore_volume.cube

"""

from numpy import pi as _pi


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
        The dictionary key of the pore diameter values

    Returns
    -------
    value : NumPy ndarray
        Array containing pore volume values.

    """
    return _pi/6*target[pore_diameter]**3


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
        The dictionary key of the pore diameter values

    Returns
    -------
    value : NumPy ndarray
        Array containing pore volume values.

    """
    return target[pore_diameter]**3


def circle(target, pore_diameter='pore.diameter'):
    r"""
    Calculate pore volume from diameter assuming a spherical pore body

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    pore_diameter : string
        The dictionary key of the pore diameter values

    """
    return _pi/4 * target[pore_diameter]**2


def square(target, pore_diameter='pore.diameter'):
    r"""
    Calculate pore volume from diameter assuming a cubic pore body

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    pore_diameter : string
        The dictionary key of the pore diameter values

    """
    return target[pore_diameter]**2


def cylinder(target, pore_diameter='pore.diameter'):
    r"""
    Calculate pore volume from diameter assuming a cylindrical pore body
    with a height equal to its diameter.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    pore_diameter : string
        The dictionary key of the pore diameter values

    Returns
    -------
    value : NumPy ndarray
        Array containing pore volume values.

    """
    diams = target[pore_diameter]
    value = _pi/4*diams**2*diams
    return value
