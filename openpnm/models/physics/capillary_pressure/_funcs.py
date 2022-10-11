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
]


@_doctxt
def washburn(phase,
             surface_tension="throat.surface_tension",
             contact_angle="throat.contact_angle",
             diameter="throat.diameter"):
    r"""
    Computes the capillary entry pressure assuming the throat in a
    cylindrical tube.

    Parameters
    ----------
    %(phase)s
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
    network = phase.network
    sigma = phase[surface_tension]
    theta = phase[contact_angle]
    r = network[diameter] / 2
    value = -2 * sigma * _np.cos(_np.radians(theta)) / r
    if diameter.split(".")[0] == "throat":
        pass
    else:
        value = value[phase.pores()]
    value[_np.absolute(value) == _np.inf] = 0
    return value


@_doctxt
def purcell(phase,
            r_toroid,
            surface_tension="throat.surface_tension",
            contact_angle="throat.contact_angle",
            diameter="throat.diameter"):
    r"""
    Computes the throat capillary entry pressure assuming the throat is a
    toroid.

    Parameters
    ----------
    %(phase)s
    r_toroid : float or array_like
        The radius of the toroid surrounding the pore
    surface_tension : str
        %(dict_blurb)s surface tension.
    contact_angle : str
        %(dict_blurb)s contact angle.
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
    network = phase.network
    sigma = phase[surface_tension]
    theta = phase[contact_angle]
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
        pass
    else:
        value = value[phase.pores()]
    return value
