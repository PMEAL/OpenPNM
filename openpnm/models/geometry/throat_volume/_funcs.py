import numpy as _np
from openpnm.utils import Docorator


__all__ = ["cylinder",
           "cuboid",
           "rectangle",
           "extrusion",
           "lens",
           "pendular_ring"]
docstr = Docorator()


@docstr.get_sections(base='models.geometry.throat_volume',
                     sections=['Parameters', 'Returns'])
@docstr.dedent
def cylinder(target, throat_length='throat.length',
             throat_diameter='throat.diameter'):
    r"""
    Calculate throat volume assuing a cylindrical shape

    Parameters
    ----------
    %(models.target.parameters)s
    %(models.geometry.tlen)s
    %(models.geometry.tdia)s

    Returns
    -------
    volumes : ndarray
        A numpy ndarray containing throat volume values

    Notes
    -----
    This models does not account for the volume reprsented by the
    intersection of the throat with a spherical pore body.  Use the ``lens``
    or ``pendular_ring`` models in addition to this one to account for this
    volume.

    """
    leng = target[throat_length]
    diam = target[throat_diameter]
    value = _np.pi/4*leng*diam**2
    return value


@docstr.dedent
def cuboid(target, throat_length='throat.length',
           throat_diameter='throat.diameter'):
    r"""
    Calculate throat volume assuing a square cross-section

    Parameters
    ----------
    %(models.geometry.throat_volume.parameters)s

    Returns
    -------
    %(models.geometry.throat_volume.returns)s

    Notes
    -----
    At present this models does NOT account for the volume reprsented by the
    intersection of the throat with a spherical pore body.

    """
    leng = target[throat_length]
    diam = target[throat_diameter]
    value = leng*diam**2
    return value


@docstr.dedent
def rectangle(target, throat_length='throat.length',
              throat_diameter='throat.diameter'):
    r"""
    Calculate throat volume assuing a rectangular shape

    Parameters
    ----------
    %(models.geometry.throat_volume.parameters)s

    Returns
    -------
    %(models.geometry.throat_volume.returns)s

    Notes
    -----
    At present this models does NOT account for the volume reprsented by the
    intersection of the throat with a spherical pore body.
    """
    return target[throat_length] * target[throat_diameter]


@docstr.dedent
def extrusion(target, throat_length='throat.length',
              throat_area='throat.cross_sectional_area'):
    r"""
    Calculate throat volume from the throat area and the throat length. This
    method is useful for abnormal shaped throats.

    Parameters
    ----------
    ----------
    %(models.target.parameters)s
    %(models.geometry.tlen)s
    %(models.geometry.tarea)s

    Returns
    -------
    %(models.geometry.throat_volume.returns)s

    Notes
    -----
    At present this models does NOT account for the volume reprsented by the
    intersection of the throat with a spherical pore body.

    """
    leng = target[throat_length]
    area = target[throat_area]
    value = leng*area
    return value


@docstr.dedent
def lens(target, throat_diameter='throat.diameter',
         pore_diameter='pore.diameter'):
    r"""
    Calculates the volume residing the hemispherical caps formed by the
    intersection between cylindrical throats and spherical pores.

    This volume should be subtracted from throat volumes if the throat lengths
    were found using throat end points.

    Parameters
    ----------
    %(models.target.parameters)s
    %(models.geometry.tdia)s
    %(models.geometry.pdia)s

    Returns
    -------
    %(models.geometry.throat_volume.returns)s

    Notes
    -----
    This model does not consider the possibility that multiple throats might
    overlap in the same location which could happen if throats are large and
    connectivity is random.

    See Also
    --------
    pendular_ring
    """
    network = target.network
    conns = network['throat.conns']
    Rp = target[pore_diameter]/2
    Rt = target[throat_diameter]/2
    a = _np.atleast_2d(Rt).T
    q = _np.arcsin(a/Rp[conns])
    b = Rp[conns]*_np.cos(q)
    h = Rp[conns] - b
    V = 1/6*_np.pi*h*(3*a**2 + h**2)
    return _np.sum(V, axis=1)


@docstr.dedent
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
    %(models.target.parameters)s
    %(models.geometry.tdia)s
    %(models.geometry.pdia)s

    Returns
    -------
    %(models.geometry.throat_volume.returns)s

    Notes
    -----
    This model does not consider the possibility that multiple throats might
    overlap in the same location which could happen if throats are large and
    connectivity is random.

    See Also
    --------
    lens
    """
    network = target.network
    conns = network['throat.conns']
    Rp = target[pore_diameter]/2
    Rt = target[throat_diameter]/2
    a = _np.atleast_2d(Rt).T
    q = _np.arcsin(a/Rp[conns])
    b = Rp[conns]*_np.cos(q)
    h = Rp[conns] - b
    Vlens = 1/6*_np.pi*h*(3*a**2 + h**2)
    Vcyl = _np.pi*(a)**2*h
    V = Vcyl - Vlens
    return _np.sum(V, axis=1)
