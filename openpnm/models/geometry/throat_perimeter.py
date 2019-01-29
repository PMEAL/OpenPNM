r"""

.. autofunction:: openpnm.models.geometry.throat_perimeter.cylinder
.. autofunction:: openpnm.models.geometry.throat_perimeter.cuboid

"""
import scipy as _sp


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
    """
    return target[throat_diameter]*_sp.pi


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
        The dictionary key of the array containing the throat diameter values
    """
    return target[throat_diameter]*4
