import logging
import numpy as _np
from transforms3d import _gohlketransforms as tr
from openpnm.models import physics as pm
from openpnm.models import _doctxt
from openpnm.models.physics._utils import _get_key_props


logger = logging.getLogger(__name__)

__all__ = [
    "washburn",
    "purcell",
    "ransohoff_snap_off",
    "purcell_bidirectional",
    "sinusoidal_bidirectional"
]


@_doctxt
def washburn(target,
             surface_tension="throat.surface_tension",
             contact_angle="throat.contact_angle",
             diameter="throat.diameter"):
    r"""
    Computes the capillary entry pressure assuming the throat in a
    cylindrical tube.

    Parameters
    ----------
    %(target_blurb)s
    surface_tension : str
        %(dict_blurb)s surface tension. If a pore property is given, it is
        interpolated to a throat list.
    contact_angle : str
        %(dict_blurb)s contact angle. If a pore property is given, it is
        interpolated to a throat list.
    diameter : str
        %(dict_blurb)s throat diameter

    Returns
    -------
    %(return_arr)s capillary entry pressure

    Notes
    -----
    The Washburn equation is:

    .. math::
        P_c = -\frac{2\sigma(cos(\theta))}{r}

    This is the most basic approach to calculating entry pressure and is
    suitable for highly non-wetting invading phases in most materials.

    """
    network = target.network
    phase = target
    sigma = phase[surface_tension]
    theta = phase[contact_angle]
    r = network[diameter] / 2
    value = -2 * sigma * _np.cos(_np.radians(theta)) / r
    if diameter.split(".")[0] == "throat":
        value = value[phase.throats(target.name)]
    else:
        value = value[phase.pores(target.name)]
    value[_np.absolute(value) == _np.inf] = 0
    return value


@_doctxt
def purcell(target,
            r_toroid,
            surface_tension="pore.surface_tension",
            contact_angle="pore.contact_angle",
            diameter="throat.diameter"):
    r"""
    Computes the throat capillary entry pressure assuming the throat is a
    toroid.

    Parameters
    ----------
    %(target_blurb)s
    r_toroid : float or array_like
        The radius of the toroid surrounding the pore
    surface_tension : str
        %(dict_blurb)s surface tension. If a pore property is given, it is
        interpolated to a throat list.
    contact_angle : str
        %(dict_blurb)s contact angle. If a pore property is given, it is
        interpolated to a throat list.
    diameter : str
        %(dict_blurb)s throat diameter

    Returns
    -------
    %(return_arr)s capillary entry pressure

    Notes
    -----
    This approach accounts for the converging-diverging nature of many throat
    types.  Advancing the meniscus beyond the apex of the toroid requires an
    increase in capillary pressure beyond that for a cylindical tube of the
    same radius. The details of this equation are described by Mason and
    Morrow [1]_, and explored by Gostick [2]_ in the context of a pore network
    model.

    References
    ----------
    .. [1] G. Mason, N. R. Morrow, Effect of contact angle on capillary
           displacement curvatures in pore throats formed by spheres. J.
           Colloid Interface Sci. 168, 130 (1994).
    .. [2] J. Gostick, Random pore network modeling of fibrous PEMFC gas
           diffusion media using Voronoi and Delaunay tessellations. J.
           Electrochem. Soc. 160, F731 (2013).

    """
    network = target.project.network
    phase = target
    element, sigma, theta = _get_key_props(
        phase=phase,
        diameter=diameter,
        surface_tension=surface_tension,
        contact_angle=contact_angle,
    )
    r = network[diameter] / 2
    R = r_toroid
    alpha = (
        theta - 180 + _np.rad2deg(_np.arcsin(_np.sin(_np.radians(theta)) / (1 + r / R)))
    )
    value = (-2 * sigma / r) * (
        _np.cos(_np.radians(theta - alpha))
        / (1 + R / r * (1 - _np.cos(_np.radians(alpha))))
    )
    if diameter.split(".")[0] == "throat":
        value = value[phase.throats(target.name)]
    else:
        value = value[phase.pores(target.name)]
    return value


@_doctxt
def ransohoff_snap_off(target,
                       shape_factor=2.0,
                       wavelength=5e-6,
                       require_pair=False,
                       surface_tension="pore.surface_tension",
                       contact_angle="pore.contact_angle",
                       diameter="throat.diameter",
                       vertices="throat.offset_vertices",
                       **kwargs):
    r"""
    Computes the capillary snap-off pressure assuming the throat is
    cylindrical with converging-diverging change in diamater - like the
    Purcell model. The wavelength of the change in diamater is the fiber
    radius.

    Parameters
    ----------
    %(target_blurb)s
    shape_factor : float
        A constant dependent on the shape of throat cross-section
        1.75 - 2.0, see Ref [1]
    wavelength : float or array like
        The transverse interfacial radius of curvature at the neck
        (fiber radius in fibrous media)
    require_pair : bool
        Controls whether snap-off requires a pair of arc meniscii to occur.
    surface_tension : str
        %(dict_blurb)s surface tension. If a pore property is given, it is
        interpolated to a throat list.
    contact_angle : str
        %(dict_blurb)s contact angle. If a pore property is given, it is
        interpolated to a throat list.
    diameter : str
        %(dict_blurb)s throat diameter

    Returns
    -------
    %(return_arr)s capillary entry pressure

    References
    ----------
    [1]: Ransohoff, T.C., Gauglitz, P.A. and Radke, C.J., 1987. Snap‐off of gas
    bubbles in smoothly constricted noncircular capillaries. AIChE Journal,
    33(5), pp.753-765.

    """
    phase = target
    geometry = target.project.find_geometry(target)
    element, sigma, theta = _get_key_props(
        phase=phase,
        diameter=diameter,
        surface_tension=surface_tension,
        contact_angle=contact_angle,
    )
    try:
        all_verts = geometry[vertices]
        # Work out whether throat geometry can support at least one pair of
        # adjacent arc menisci that can grow and merge to form snap-off
        # Only works if throat vertices are in convex hull order
        angles_ok = _np.zeros(geometry.Nt, dtype=bool)
        for T in range(geometry.Nt):
            verts = all_verts[T]
            x = verts[:, 0]
            y = verts[:, 1]
            z = verts[:, 2]
            # PLus
            p = 1
            # Minus
            m = -1
            verts_p = _np.vstack((_np.roll(x, p), _np.roll(y, p), _np.roll(z, p))).T
            verts_m = _np.vstack((_np.roll(x, m), _np.roll(y, m), _np.roll(z, m))).T
            v1 = verts_p - verts
            v2 = verts_m - verts
            corner_angles = _np.rad2deg(tr.angle_between_vectors(v1, v2, axis=1))
            # Logical test for existence of arc menisci
            am = theta[T] <= 90 - corner_angles / 2
            if require_pair:
                # Logical test for two adjacent arc menisci
                pair_p = _np.logical_and(am, _np.roll(am, +p))
                pair_m = _np.logical_and(am, _np.roll(am, +m))
                am_pair = _np.any(_np.logical_or(pair_p, pair_m))
                angles_ok[T] = am_pair
            else:
                # Logical test for any arc menisci
                angles_ok[T] = _np.any(am)
    except Exception:
        logger.warning("Model is designed to work with property: " + vertices)
        angles_ok = _np.ones(geometry.Nt, dtype=bool)

    # Condition for arc menisci to form in corners
    rad_Ts = geometry[diameter] / 2
    # Ransohoff and Radke eq. 4
    C = 1 / rad_Ts - 1 / wavelength
    value = sigma[phase.throats(target.name)] * C
    # Only throats that can support arc menisci can snap-off
    value[~angles_ok] = _np.nan
    logger.info(
        "Snap off pressures calculated for "
        + str(_np.around(100 * _np.sum(angles_ok) / _np.size(angles_ok), 0))
        + "% of throats"
    )
    return value


@_doctxt
def purcell_bidirectional(target,
                          r_toroid=5e-6,
                          num_points=1000,
                          surface_tension="pore.surface_tension",
                          contact_angle="pore.contact_angle",
                          throat_diameter="throat.diameter",
                          pore_diameter="pore.diameter"):
    r"""
    Computes the throat capillary entry pressure assuming the throat is a
    toroid. Makes use of the toroidal meniscus model with mode touch.
    This model accounts for mensicus protrusion into adjacent pores and
    touching solid features.

    It is bidirectional becauase the connected pores generally have
    different sizes and this determines how far the meniscus can protrude.

    Parameters
    ----------
    %(target_blurb)s
    r_toroid : float or array_like
        The radius of the toroid surrounding the pore
    num_points : float, default 100
        The number of divisions to make along the profile length to assess the
        meniscus properties in order to find the touch length.
    surface_tension : str
        %(dict_blurb)s surface tension. If a pore property is given, it is
        interpolated to a throat list.
    contact_angle : str
        %(dict_blurb)s contact angle. If a pore property is given, it is
    interpolated to a throat list.
    throat_diameter : str
        %(dict_blurb)s throat diameter
    pore_diameter : str
        %(dict_blurb)s pore diameter

    Returns
    -------
    %(return_arr)s capillary entry pressure

    """
    network = target.project.network
    conns = network["throat.conns"]
    values = {}
    for p in range(2):
        network["throat.temp_diameter"] = network[pore_diameter][conns[:, p]]
        key = "throat.touch_pore_" + str(p)
        target.add_model(
            propname=key,
            model=pm.meniscus.purcell,
            mode="touch",
            r_toroid=r_toroid,
            num_points=num_points,
            throat_diameter=throat_diameter,
            surface_tension=surface_tension,
            contact_angle=contact_angle,
            touch_length="throat.temp_diameter",
        )
        values[p] = target[key]
        target.remove_model(key)
    del network["throat.temp_diameter"]
    return _np.vstack((values[0], values[1])).T


@_doctxt
def sinusoidal_bidirectional(target,
                             r_toroid=5e-6,
                             num_points=1e3,
                             surface_tension="pore.surface_tension",
                             contact_angle="pore.contact_angle",
                             throat_diameter="throat.diameter",
                             pore_diameter="pore.diameter"):
    r"""
    Computes the throat capillary entry pressure assuming the throat has a
    sinusoisal profile.

    Makes use of the toroidal meniscus model with mode touch. This model
    accounts for mensicus protrusion into adjacent pores and touching
    solid features. It is bidirectional becauase the connected pores
    generally have different sizes and this determines how far the
    meniscus can protrude.

    Parameters
    ----------
    %(target_blurb)s
    r_toroid : float or array_like
        The radius of the toroid surrounding the pore
    num_points : float, default 100
        The number of divisions to make along the profile length to assess the
        meniscus properties in order to find the touch length.
    surface_tension : str
        %(dict_blurb)s surface tension. If a pore property is given, it is
        interpolated to a throat list.
    contact_angle : str
        %(dict_blurb)s contact angle. If a pore property is given, it is
        interpolated to a throat list.
    throat_diameter : str
        %(dict_blurb)s throat diameter
    pore_diameter : str
        %(dict_blurb)s pore diameter

    Returns
    -------
    %(return_arr)s capillary entry pressure

    """
    network = target.project.network
    conns = network["throat.conns"]
    values = {}
    for p in range(2):
        network["throat.temp_diameter"] = network[pore_diameter][conns[:, p]]
        key = "throat.touch_pore_" + str(p)
        target.add_model(
            propname=key,
            model=pm.meniscus.sinusoidal,
            mode="touch",
            r_toroid=r_toroid,
            num_points=num_points,
            surface_tension=surface_tension,
            contact_angle=contact_angle,
            throat_diameter=throat_diameter,
            touch_length="throat.temp_diameter",
        )
        values[p] = target[key]
        target.remove_model(key)
    del network["throat.temp_diameter"]
    return _np.vstack((values[0], values[1])).T
