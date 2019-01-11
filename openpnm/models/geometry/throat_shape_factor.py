r"""

.. autofunction:: openpnm.models.geometry.throat_shape_factor.compactness
.. autofunction:: openpnm.models.geometry.throat_shape_factor.mason_morrow
.. autofunction:: openpnm.models.geometry.throat_shape_factor.jenkins_rao

"""

import scipy as _sp


def compactness(target, throat_perimeter='throat.perimeter',
                throat_area='throat.area'):
    r"""
    Mortensen et al. have shown that the Hagen-Poiseuille hydraluic resistance
    is linearly dependent on the compactness. Defined as perimeter^2/area.
    The dependence is not universal as shapes with sharp corners provide more
    resistance than those that are more elliptical. Count the number of
    vertices and apply the right correction.

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
    C = _sp.ones(target.num_throats())
    C[ts] = P[ts]**2/A[ts]
    alpha = _sp.ones_like(C)*8*_sp.pi
    if 'throat.offset_vertices' in target.props():
        verts = target['throat.offset_vertices']
        for i in ts:
            if ~_sp.any(_sp.isnan(verts[i])):
                if len(verts[i]) == 3:
                    # Triangular Correction
                    alpha[i] = C[i]*(25/17) + (40*_sp.sqrt(3)/17)
                elif len(verts[i]) == 4:
                    # Rectangular Correction
                    alpha[i] = C[i]*(22/7) - (65/3)
                elif len(verts[i]) > 4:
                    # Approximate Elliptical Correction
                    alpha[i] = C[i]*(8/3) - (8*_sp.pi/3)
    # For a perfect circle alpha = 8*pi so normalize by this
    alpha /= 8*_sp.pi
    # Very small throats could have values less than one
    alpha[alpha < 1.0] = 1.0
    return alpha


def mason_morrow(target, throat_perimeter='throat.perimeter',
                 throat_area='throat.area'):
    r"""
    Mason and Morrow relate the capillary pressure to the shaped factor in a
    similar way to Mortensen but for triangles.

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
    value[ts] = 1/(4*_sp.pi)
    return value


def jenkins_rao(target, throat_perimeter='throat.perimeter',
                throat_area='throat.area',
                throat_diameter='throat.indiameter'):
    r"""
    Jenkins and Rao relate the capillary pressure in an eliptical throat to
    the aspect ratio

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
