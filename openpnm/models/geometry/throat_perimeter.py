r"""
"""
import numpy as _np


def cylinder(target, throat_diameter='throat.diameter'):
    r"""
    Calcuate the throat perimeter assuming a circular cross-section

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    throat_diameter : string
        The dictionary key of the array containing the throat diameter values

    Returns
    -------
    value : NumPy ndarray
        Array containing throat perimeter values.

    """
    return target[throat_diameter]*_np.pi


def cuboid(target, throat_diameter='throat.diameter'):
    r"""
    Calcuate the throat perimeter assuming a square cross-section

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    throat_diameter : string
        The dictionary key of the array containing the throat diameter values.

    Returns
    -------
    value : NumPy ndarray
        Array containing throat perimeter values.

    """
    return target[throat_diameter]*4


def rectangle(target, throat_diameter='throat.diameter'):
    r"""
    Calcuate the throat perimeter assuming a rectangular cross-section (2D)

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    throat_diameter : string
        The dictionary key of the array containing the throat diameter values

    Returns
    -------
    value : NumPy ndarray
        Array containing throat perimeter values.

    """
    return 1.0
