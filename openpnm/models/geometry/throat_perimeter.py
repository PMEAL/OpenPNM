import numpy as _np
from openpnm.utils import Docorator


docstr = Docorator()


@docstr.get_sections(base='models.geometry.throat_perimeter',
                     sections=['Parameters', 'Returns'])
def cylinder(target, throat_diameter='throat.diameter'):
    r"""
    Calcuate the throat perimeter assuming a circular cross-section

    Parameters
    ----------
    target : OpenPNM Base object
        Object with which this model is associated. This controls
        the length of the calculated array, and also provides access to
        other necessary properties.
    throat_diameter : str
        Name of the dictionary key on ``target`` where the array containing
        throat diameter values is stored

    Returns
    -------
    perimeters : ndarray
        A numpy ndarray containing throat perimeter values

    """
    return target[throat_diameter]*_np.pi


@docstr.dedent
def cuboid(target, throat_diameter='throat.diameter'):
    r"""
    Calcuate the throat perimeter assuming a square cross-section

    Parameters
    ----------
    %(models.geometry.throat_perimeter.parameters)s

    Returns
    -------
    %(models.geometry.throat_perimeter.results)s

    """
    return target[throat_diameter]*4


def rectangle(target, throat_diameter='throat.diameter'):
    r"""
    Calcuate the throat perimeter assuming a rectangular cross-section (2D)

    Parameters
    ----------
    %(models.geometry.throat_perimeter.parameters)s

    Returns
    -------
    %(models.geometry.throat_perimeter.results)s

    """
    return 1.0
