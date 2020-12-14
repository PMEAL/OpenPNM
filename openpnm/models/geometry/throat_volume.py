r"""
"""
import numpy as _np


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

    Returns
    -------
    value : NumPy ndarray
        Array containing throat volume values.

    Notes
    -----
    At present this models does NOT account for the volume reprsented by the
    intersection of the throat with a spherical pore body.

    """
    leng = target[throat_length]
    diam = target[throat_diameter]
    value = _np.pi/4*leng*diam**2
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

    Returns
    -------
    value : NumPy ndarray
        Array containing throat volume values.

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

    Returns
    -------
    value : NumPy ndarray
        Array containing throat volume values.

    Notes
    -----
    At present this models does NOT account for the volume reprsented by the
    intersection of the throat with a spherical pore body.

    """
    leng = target[throat_length]
    area = target[throat_area]
    value = leng*area
    return value


def rectangle(target, throat_length='throat.length',
              throat_diameter='throat.diameter'):
    r"""
    Calculate throat volume assuing a rectangular shape

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
    return target[throat_length] * target[throat_diameter]


def lens(target, throat_diameter='throat.diameter',
         pore_diameter='pore.diameter'):
    r"""
    Calculates the volume residing the hemispherical caps formed by the
    intersection between cylindrical throats and spherical pores.

    This volume should be subtracted from throat volumes if the throat lengths
    were found using throat end points.

    Parameters
    ----------
    throat_diameter : string
        The dictionary keys containing the array with the throat diameter
        values.
    pore_diameter : string
        The dictionary keys containing the array with the pore diameter
        values.

    Returns
    -------
    volume : ND-array
        The volume that should be subtracted from each throat volume to prevent
        double counting the volume of overlapping area.

    Notes
    -----
    This model does NOT consider the possibility that multiple throats might
    overlap in the same location which could happen if throats are large and
    connectivity is random.

    See Also
    --------
    pendular_ring
    """
    network = target.network
    conns = network['throat.conns']
    Rp = target[pore_diameter]
    Rt = target[throat_diameter]
    a = _np.atleast_2d(Rt).T
    q = _np.arcsin(a/Rp[conns])
    b = Rp[conns]*_np.cos(q)
    h = Rp[conns] - b
    V = 1/6*_np.pi*h*(3*a**2 + h**2)
    return _np.sum(V, axis=1)


def pendular_ring(target, throat_diameter='throat.diameter',
                  pore_diameter='pore.diameter'):
    r"""
    Calculates the volume of the pendular rings residing between the end of
    a cylindrical throat and spherical pores that are in contact but not
    overlapping.

    This volume should be added to the throat volume if the throat length was
    found as the center-to-center distance less the pore radii.

    Parameters
    ----------
    throat_diameter : string
        The dictionary keys containing the array with the throat diameter
        values.
    pore_diameter : string
        The dictionary keys containing the array with the pore diameter
        values.

    Returns
    -------
    volume : ND-array
        The volume that should be added to each throat volume to account for
        under-represented void volume at the pore-throat junctions.

    Notes
    -----
    This model does NOT consider the possibility that multiple throats might
    overlap in the same location which could happen if throats are large and
    connectivity is random.

    See Also
    --------
    lens
    """
    network = target.network
    conns = network['throat.conns']
    Rp = target[pore_diameter]
    Rt = target[throat_diameter]
    a = _np.atleast_2d(Rt).T
    q = _np.arcsin(a/Rp[conns])
    b = Rp[conns]*_np.cos(q)
    h = Rp[conns] - b
    Vlens = 1/6*_np.pi*h*(3*a**2 + h**2)
    Vcyl = _np.pi*(a)**2*h
    V = Vcyl - Vlens
    return _np.sum(V, axis=1)
