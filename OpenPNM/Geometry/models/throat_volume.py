r"""
===============================================================================
Submodule -- throat_volume
===============================================================================

"""
import scipy as _sp


def cylinder(geometry, throat_length='throat.length',
             throat_diameter='throat.diameter', **kwargs):
    r"""
    Calculate throat volume assuing a cylindrical shape

    Parameters
    ----------
    geometry : OpenPNM Geometry object
        The Geometry object which this model is associated.  This is needed
        to access the other relevant geometrical properties for this
        calculation.

    throat_length and throat_diameter : strings
        The dictionary keys containing the arrays with the throat diameter and
        length values.

    Notes
    -----
    At present this models does NOT account for the volume reprsented by the
    intersection of the throat with a spherical pore body.
    """
    leng = geometry[throat_length]
    diam = geometry[throat_diameter]
    value = _sp.pi/4*leng*diam**2
    return value


def cuboid(geometry, throat_length='throat.length',
           throat_diameter='throat.diameter', **kwargs):
    r"""
    Calculate throat volume assuing a square cross-section

    Parameters
    ----------
    geometry : OpenPNM Geometry object
        The Geometry object which this model is associated.  This is needed
        to access the other relevant geometrical properties for this
        calculation.

    throat_length and throat_diameter : strings
        The dictionary keys containing the arrays with the throat diameter and
        length values.

    Notes
    -----
    At present this models does NOT account for the volume reprsented by the
    intersection of the throat with a spherical pore body.
    """
    leng = geometry[throat_length]
    diam = geometry[throat_diameter]
    value = leng*diam**2
    return value


def extrusion(geometry, throat_length='throat.length',
              throat_area='throat.area', **kwargs):
    r"""
    Calculate throat volume from the throat area and the throat length. This
    method is useful for abnormal shaped throats.

    Parameters
    ----------
    geometry : OpenPNM Geometry object
        The Geometry object which this model is associated.  This is needed
        to access the other relevant geometrical properties for this
        calculation.

    throat_length and throat_area : strings
        The dictionary keys containing the arrays with the throat area and
        length values.

    Notes
    -----
    At present this models does NOT account for the volume reprsented by the
    intersection of the throat with a spherical pore body.
        """
    leng = geometry[throat_length]
    area = geometry[throat_area]
    value = leng*area
    return value
