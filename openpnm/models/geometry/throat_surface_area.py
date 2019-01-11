r"""

.. autofunction:: openpnm.models.geometry.throat_surface_area.cylinder
.. autofunction:: openpnm.models.geometry.throat_surface_area.cuboid
.. autofunction:: openpnm.models.geometry.throat_surface_area.extrusion

"""
import scipy as _sp


def cylinder(target, throat_diameter='throat.diameter',
             throat_length='throat.length'):
    r"""
    Calculate surface area for a cylindrical throat

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    throat_diameter : string
        Dictionary key to the throat diameter array.  Default is
        'throat.diameter'.

    throat_length : string
        Dictionary key to the throat length array.  Default is 'throat.length'.
    """
    D = target[throat_diameter]
    L = target[throat_length]
    value = _sp.pi*D*L
    return value


def cuboid(target, throat_diameter='throat.diameter',
           throat_length='throat.length'):
    r"""
    Calculate surface area for a cuboid throat

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    throat_diameter : string
        Dictionary key to the throat diameter array.  Default is
        'throat.diameter'.

    throat_length : string
        Dictionary key to the throat length array.  Default is 'throat.length'.
    """
    D = target[throat_diameter]
    L = target[throat_length]
    value = 4*D*L
    return value


def extrusion(target, throat_perimeter='throat.perimeter',
              throat_length='throat.length'):
    r"""
    Calculate surface area for an arbitrary shaped throat give the perimeter
    and length.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    throat_perimeter : string
        Dictionary key to the throat perimeter array.  Default is
        'throat.perimeter'.

    throat_length : string
        Dictionary key to the throat length array.  Default is 'throat.length'.

    """
    P = target[throat_perimeter]
    L = target[throat_length]
    value = P*L
    return value
