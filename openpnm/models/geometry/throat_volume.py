r"""

.. autofunction:: openpnm.models.geometry.throat_volume.cylinder
.. autofunction:: openpnm.models.geometry.throat_volume.cuboid
.. autofunction:: openpnm.models.geometry.throat_volume.extrusion

"""
import scipy as _sp


def cylinder(target, throat_length='throat.length',
             throat_diameter='throat.diameter'):
    r"""
    Calculate throat volume assuing a cylindrical shape

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    throat_length and throat_diameter : strings
        The dictionary keys containing the arrays with the throat diameter and
        length values.

    Notes
    -----
    At present this models does NOT account for the volume reprsented by the
    intersection of the throat with a spherical pore body.
    """
    leng = target[throat_length]
    diam = target[throat_diameter]
    value = _sp.pi/4*leng*diam**2
    return value


def cuboid(target, throat_length='throat.length',
           throat_diameter='throat.diameter'):
    r"""
    Calculate throat volume assuing a square cross-section

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    throat_length and throat_diameter : strings
        The dictionary keys containing the arrays with the throat diameter and
        length values.

    Notes
    -----
    At present this models does NOT account for the volume reprsented by the
    intersection of the throat with a spherical pore body.
    """
    leng = target[throat_length]
    diam = target[throat_diameter]
    value = leng*diam**2
    return value


def extrusion(target, throat_length='throat.length',
              throat_area='throat.area'):
    r"""
    Calculate throat volume from the throat area and the throat length. This
    method is useful for abnormal shaped throats.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    throat_length and throat_area : strings
        The dictionary keys containing the arrays with the throat area and
        length values.

    Notes
    -----
    At present this models does NOT account for the volume reprsented by the
    intersection of the throat with a spherical pore body.
        """
    leng = target[throat_length]
    area = target[throat_area]
    value = leng*area
    return value
