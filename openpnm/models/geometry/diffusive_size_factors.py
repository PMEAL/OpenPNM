import numpy as _np
import openpnm.models.geometry.conduit_lengths as _conduit_lengths
import openpnm.geometry.GenericGeometry as _GenericGeometry

__all__ = [
    "spheres_and_cylinders",
    "circles_and_rectangles",
    "cones_and_cylinders",
    "pyramids_and_cuboids",
    "cubes_and_cuboids",
    "intersecting_cones",
    "intersecting_pyramids",
    "ncylinders_in_series"
]


def spheres_and_cylinders(
    target: _GenericGeometry,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
):
    r"""
    Computes diffusive size factor for conduits of spheres and cylinders.

    Parameters
    ----------
    target : _GenericGeometry
        Geometry object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.
    pore_diameter : str
        Dictionary key of the pore diameter values.
    throat_diameter : str
        Dictionary key of the throat diameter values.

    Notes
    -----
    The diffusive size factor is the geometrical part of the pre-factor in
    Fick's law:

    .. math::

        n_A = S_{diffusive} D_{AB} \Delta C_A

    Thus :math:`S_{diffusive}` represents the combined effect of the area and
    length of the *conduit*, which consists of a throat and 1/2 of the pore
    on each end.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    cn = network["throat.conns"][throats]

    D1 = network[pore_diameter][cn[:, 0]]
    D2 = network[pore_diameter][cn[:, 1]]
    Dt = network[throat_diameter][throats]
    L1, Lt, L2 = _conduit_lengths.spheres_and_cylinders(
        target, pore_diameter=pore_diameter, throat_diameter=throat_diameter
    ).T

    # Fi is the integral of (1/A) dx, x = [0, Li]
    F1 = 2 / (D1 * _np.pi) * _np.arctanh(2 * L1 / D1)
    F2 = 2 / (D2 * _np.pi) * _np.arctanh(2 * L2 / D2)
    Ft = Lt / (_np.pi / 4 * Dt ** 2)
    g1, g2, gt = 1 / F1, 1 / F2, 1 / Ft

    return {"pore1": g1, "throat": gt, "pore2": g2}


def circles_and_rectangles(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
):
    r"""
    Compute diffusive shape coefficient for conduits of circles and rectangles

    Parameters
    ----------
    target : _GenericGeometry
        Geometry object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.
    pore_diameter : str
        Dictionary key of the pore diameter values.
    throat_diameter : str
        Dictionary key of the throat diameter values.

    Notes
    -----
    The diffusive shape coefficient is the geometrical part of the pre-factor
    in Fick's law:

    .. math::

        n_A = S_{diffusive} D_{AB} \Delta C_A

    Thus :math:`S_{diffusive}` represents the combined effect of the area and
    length of the *conduit*, which consists of a throat and 1/2 of the pore
    on each end.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    cn = network["throat.conns"][throats]

    D1 = network[pore_diameter][cn[:, 0]]
    D2 = network[pore_diameter][cn[:, 1]]
    Dt = network[throat_diameter][throats]
    L1, Lt, L2 = _conduit_lengths.circles_and_rectangles(
        target, pore_diameter=pore_diameter, throat_diameter=throat_diameter
    ).T

    # Fi is the integral of (1/A) dx, x = [0, Li]
    F1 = 0.5 * _np.arctanh(2 * L1 / _np.sqrt(D1 ** 2 - 4 * L1 ** 2))
    F2 = 0.5 * _np.arctanh(2 * L2 / _np.sqrt(D2 ** 2 - 4 * L2 ** 2))
    Ft = Lt / Dt
    g1, g2, gt = 1 / F1, 1 / F2, 1 / Ft

    return {"pore1": g1, "throat": gt, "pore2": g2}


def cones_and_cylinders(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
    pore_area="pore.area",
    throat_area="throat.area",
    conduit_lengths="throat.conduit_lengths",
):
    r"""
    Computes diffusive shape coefficient assuming pores are truncated pyramids
    and throats are cylinders.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.
    pore_diameter : string
        Dictionary key of the pore diameter values.
    throat_diameter : string
        Dictionary key of the throat diameter values.
    pore_area : string
        Dictionary key of the pore area values.
    throat_area : string
        Dictionary key of the throat area values.
    conduit_lengths : string
        Dictionary key of the conduit lengths' values.

    Returns
    -------
    A dictionary containing the diffusive_shape_coefficient, which can be
    accessed via the dict keys 'pore1', 'pore2', and 'throat'.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    cn = network["throat.conns"][throats]

    D1 = network[pore_diameter][cn[:, 0]]
    D2 = network[pore_diameter][cn[:, 1]]
    Dt = network[throat_diameter][throats]
    At = network[throat_area][throats]
    L1 = network[conduit_lengths + ".pore1"][throats]
    L2 = network[conduit_lengths + ".pore2"][throats]
    Lt = network[conduit_lengths + ".throat"][throats]

    # Fi is the integral of (1/A) dx, x = [0, Li]
    F1 = 4 * L1 / (D1 * Dt * _np.pi)
    F2 = 4 * L2 / (D2 * Dt * _np.pi)
    Ft = Lt / At
    g1, g2, gt = 1 / F1, 1 / F2, 1 / Ft

    return {"pore1": g1, "throat": gt, "pore2": g2}


def pyramids_and_cuboids(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
):
    r"""
    Computes diffusive size factor for conduits of pyramids and cuboids.

    Parameters
    ----------
    target : _GenericGeometry
        Geometry object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.
    pore_diameter : str
        Dictionary key of the pore diameter values.
    throat_diameter : str
        Dictionary key of the throat diameter values.

    Notes
    -----
    The diffusive size factor is the geometrical part of the pre-factor in
    Fick's law:

    .. math::

        n_A = S_{diffusive} D_{AB} \Delta C_A

    Thus :math:`S_{diffusive}` represents the combined effect of the area and
    length of the *conduit*, which consists of a throat and 1/2 of the pore
    on each end.

    """
    raise NotImplementedError


def cubes_and_cuboids(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
    pore_aspect=[1, 1, 1],
    throat_aspect=[1, 1, 1],
):
    r"""
    Computes diffusive size factor for conduits of cubes and cuboids.

    Parameters
    ----------
    target : _GenericGeometry
        Geometry object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.
    pore_diameter : str
        Dictionary key of the pore diameter values.
    throat_diameter : str
        Dictionary key of the throat diameter values.
    pore_aspect : list
        Aspect ratio of the pores
    throat_aspect : list
        Aspect ratio of the throats

    Notes
    -----
    The diffusive size factor is the geometrical part of the pre-factor in
    Fick's law:

    .. math::

        n_A = S_{diffusive} D_{AB} \Delta C_A

    Thus :math:`S_{diffusive}` represents the combined effect of the area and
    length of the *conduit*, which consists of a throat and 1/2 of the pore
    on each end.

    """
    raise NotImplementedError


def intersecting_cones(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
    midpoint=None,
):
    r"""
    Computes diffusive size factor for conduits of intersecting cones.

    Parameters
    ----------
    target : _GenericGeometry
        Geometry object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.
    pore_diameter : str
        Dictionary key of the pore diameter values.
    throat_diameter : str
        Dictionary key of the throat diameter values.
    midpoint : str, optional
        Dictionary key of the midpoint values.

    Notes
    -----
    The diffusive size factor is the geometrical part of the pre-factor in
    Fick's law:

    .. math::

        n_A = S_{diffusive} D_{AB} \Delta C_A

    Thus :math:`S_{diffusive}` represents the combined effect of the area and
    length of the *conduit*, which consists of a throat and 1/2 of the pore
    on each end.

    """
    network = target.network
    R1, R2 = (target[pore_diameter][network.conns] / 2).T
    Rt = target[throat_diameter] / 2
    Lt = 0.
    L1 = Lt - R1
    L2 = Lt - R2
    alpha1 = (R1 - Rt) / L1
    beta1 = 1 / (1 / (Rt ** 3) - 1 / (R1 ** 3))
    alpha2 = (R2 - Rt) / L2
    beta2 = 1 / (1 / (Rt ** 3) - 1 / (R2 ** 3))

    g1 = (3 * alpha1 * _np.pi / 8) * beta1
    g2 = (3 * alpha2 * _np.pi / 8) * beta2
    gt = _np.pi * Rt ** 4 / (8 * Lt)

    return {"pore1": g1, "throat": gt, "pore2": g2}


def intersecting_pyramids(
    target: _GenericGeometry,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
    midpoint=None,
):
    r"""
    Computes diffusive size factor for conduits of intersecting pyramids.

    Parameters
    ----------
    target : _GenericGeometry
        Geometry object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.
    pore_diameter : str
        Dictionary key of the pore diameter values.
    throat_diameter : str
        Dictionary key of the throat diameter values.
    midpoint : str, optional
        Dictionary key of the midpoint values.

    Notes
    -----
    The diffusive size factor is the geometrical part of the pre-factor in
    Fick's law:

    .. math::

        n_A = S_{diffusive} D_{AB} \Delta C_A

    Thus :math:`S_{diffusive}` represents the combined effect of the area and
    length of the *conduit*, which consists of a throat and 1/2 of the pore
    on each end.

    """
    raise NotImplementedError


def ncylinders_in_series(
    target: _GenericGeometry,
    pore_diameter="pore.equivalent_diameter",
    throat_diameter="throat.equivalent_diameter",
    throat_length="throat.length",
    n=5,
):
    r"""
    Computes diffusive size factors for conduits of spheres and cylinders
    with the spheres approximated as N cylinders in series.

    Parameters
    ----------
    target : _GenericGeometry
        Geometry object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.
    pore_diameter : str
        Dictionary key pointing to the pore diameter values.
    throat_diameter : str
        Dictionary key pointing to the throat diameter values.
    throat_length : str
        Dictionary key pointing to the throat_length values.
    n : int
        Number of cylindrical divisions for each pore

    Notes
    -----
    The diffusive size factor is the geometrical part of the pre-factor in
    Fick's law:

    .. math::

        n_A = S_{diffusive} D_{AB} \Delta C_A

    Thus :math:`S_{diffusive}` represents the combined effect of the area and
    length of the *conduit*, which consists of a throat and 1/2 of the pore
    on each end.

    """
    project = target.project
    network = project.network
    P12 = network["throat.conns"]
    D1, D2 = network[pore_diameter][P12].T
    Dt = network[throat_diameter]

    # Ensure throats are never bigger than connected pores
    Dt = _np.minimum(Dt, 0.99 * _np.minimum(D1, D2))
    L1 = D1 / 2 * (_np.cos(_np.arcsin(Dt / D1)))
    L2 = D2 / 2 * (_np.cos(_np.arcsin(Dt / D2)))
    dL1 = _np.linspace(0, L1, num=n, endpoint=False)
    dL2 = _np.linspace(0, L2, num=n, endpoint=False)
    r1 = D1 / 2 * _np.sin(_np.arccos(dL1 / (D1 / 2)))
    r2 = D2 / 2 * _np.sin(_np.arccos(dL2 / (D2 / 2)))

    gtemp = (_np.pi * r1**2 / (L1 / n)).T
    g1 = 1 / _np.sum(1 / gtemp, axis=1)
    gtemp = (_np.pi * r2 ** 2 / (L2 / n)).T
    g2 = 1 / _np.sum(1 / gtemp, axis=1)
    Lt = network[throat_length]
    gt = (_np.pi * (Dt / 2)**2 / (Lt)).T

    return {"pore1": g1, "throat": gt, "pore2": g2}
