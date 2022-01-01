import numpy as _np
from openpnm.utils import Docorator


__all__ = ["compactness",
           "mason_morrow",
           "jenkins_rao"]
docstr = Docorator()


@docstr.get_sections(base='models.geometry.throat_capillary_shape_factor',
                     sections=['Parameters', 'Returns'])
@docstr.dedent
def compactness(target, throat_perimeter='throat.perimeter',
                throat_area='throat.cross_sectional_area'):
    r"""
    Mortensen et al. have shown that the Hagen-Poiseuille hydraluic resistance
    is linearly dependent on the compactness. Defined as perimeter^2/area.
    The dependence is not universal as shapes with sharp corners provide more
    resistance than those that are more elliptical. Count the number of
    vertices and apply the right correction.

    Parameters
    ----------
    %(models.target.parameters)s
    throat_perimeter : string
        The dictionary key of the array containing the throat perimeter values
    throat_area : string
        The dictionary key of the array containing the throat area values

    Returns
    -------
    alpha : NumPy ndarray
        Array containing throat compactness values.

    References
    ----------
    Mortensen N.A, Okkels F., and Bruus H. Reexamination of Hagen-Poiseuille
    flow: Shape dependence of the hydraulic resistance in microchannels.
    Physical Review E, v.71, pp.057301 (2005).

    """
    # Only apply to throats with an area
    ts = target.throats()[target[throat_area] > 0]
    P = target[throat_perimeter]
    A = target[throat_area]
    C = _np.ones(target.num_throats())
    C[ts] = P[ts]**2/A[ts]
    alpha = _np.ones_like(C)*8*_np.pi
    if 'throat.offset_vertices' in target.props():
        verts = target['throat.offset_vertices']
        for i in ts:
            if ~_np.any(_np.isnan(verts[i])):
                if len(verts[i]) == 3:
                    # Triangular Correction
                    alpha[i] = C[i]*(25/17) + (40*_np.sqrt(3)/17)
                elif len(verts[i]) == 4:
                    # Rectangular Correction
                    alpha[i] = C[i]*(22/7) - (65/3)
                elif len(verts[i]) > 4:
                    # Approximate Elliptical Correction
                    alpha[i] = C[i]*(8/3) - (8*_np.pi/3)
    # For a perfect circle alpha = 8*pi so normalize by this
    alpha /= 8*_np.pi
    # Very small throats could have values less than one
    alpha[alpha < 1.0] = 1.0
    return alpha


def mason_morrow(target, throat_perimeter='throat.perimeter',
                 throat_area='throat.cross_sectional_area'):
    r"""
    Mason and Morrow relate the capillary pressure to the shape factor in a
    similar way to Mortensen but for triangles.

    Parameters
    ----------
    %(models.geometry.throat_capillary_shape_factor.parameters)s

    Returns
    -------
    %(models.geometry.throat_capillary_shape_factor.returns)s

    References
    ----------
    Mason, G. and Morrow, N.R.. Capillary behavior of a perfectly wetting
    liquid in irregular triangular tubes. Journal of Colloid and Interface
    Science, 141(1), pp.262-274 (1991).

    """
    # Only apply to throats with an area
    ts = target.throats()[target[throat_area] <= 0]
    P = target[throat_perimeter]
    A = target[throat_area]
    value = A/(P**2)
    value[ts] = 1/(4*_np.pi)
    return value


def jenkins_rao(target, throat_perimeter='throat.perimeter',
                throat_area='throat.cross_sectional_area',
                throat_diameter='throat.indiameter'):
    r"""
    Jenkins and Rao relate the capillary pressure in an eliptical throat to
    the aspect ratio

    Parameters
    ----------
    %(models.geometry.throat_capillary_shape_factor.returns)s
    %(models.geometry.tdia)s

    Returns
    -------
    %(models.geometry.throat_capillary_shape_factor.returns)s

    References
    ----------
    Jenkins, R.G. and Rao, M.B., The effect of elliptical pores on
    mercury porosimetry results. Powder technology, 38(2), pp.177-180. (1984)

    """
    P = target[throat_perimeter]
    A = target[throat_area]
    r = target[throat_diameter]/2
    # Normalized by value for perfect circle
    value = (P/A)/(2/r)
    return value
