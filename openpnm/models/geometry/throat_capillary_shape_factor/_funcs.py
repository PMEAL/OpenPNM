import numpy as _np
from openpnm.models.geometry import _geodocs


__all__ = [
    "mason_morrow",
    "jenkins_rao",
]


@_geodocs
def mason_morrow(
    network,
    throat_perimeter='throat.perimeter',
    throat_area='throat.cross_sectional_area',
):
    r"""
    Mason and Morrow relate the capillary pressure to the shape factor in a
    similar way to Mortensen but for triangles.

    Parameters
    ----------
    %(network)s
    %(Pt)s
    %(At)s

    Returns
    -------

    References
    ----------
    Mason, G. and Morrow, N.R.. Capillary behavior of a perfectly wetting
    liquid in irregular triangular tubes. Journal of Colloid and Interface
    Science, 141(1), pp.262-274 (1991).

    """
    # Only apply to throats with an area
    ts = network.throats()[network[throat_area] <= 0]
    P = network[throat_perimeter]
    A = network[throat_area]
    value = A/(P**2)
    value[ts] = 1/(4*_np.pi)
    return value


def jenkins_rao(
    network,
    throat_perimeter='throat.perimeter',
    throat_area='throat.cross_sectional_area',
    throat_diameter='throat.indiameter',
):
    r"""
    Jenkins and Rao relate the capillary pressure in an eliptical throat to
    the aspect ratio

    Parameters
    ----------
    %(network)s
    %(Pt)s
    %(At)s
    %(Dt)s

    Returns
    -------

    References
    ----------
    Jenkins, R.G. and Rao, M.B., The effect of elliptical pores on
    mercury porosimetry results. Powder technology, 38(2), pp.177-180. (1984)

    """
    P = network[throat_perimeter]
    A = network[throat_area]
    r = network[throat_diameter]/2
    # Normalized by value for perfect circle
    value = (P/A)/(2/r)
    return value
