r"""
===============================================================================
Submodule -- throat_shape_factor
===============================================================================

"""
import scipy as _sp


def compactness(geometry, throat_perimeter='throat.perimeter',
                throat_area='throat.area', **kwargs):
    r"""
    Mortensen et al. have shown that the Hagen-Poiseuille hydraluic resistance is
    linearly dependent on the compactness. Defined as perimeter^2/area.
    The dependence is not universal as shapes with sharp corners provide more
    resistance than those that are more elliptical. Count the number of vertices
    and apply the right correction.
    """
    # Only apply to throats with an area
    ts = geometry.throats()[geometry[throat_area] > 0]
    P = geometry[throat_perimeter]
    A = geometry[throat_area]
    C = _sp.ones(geometry.num_throats())
    C[ts] = P[ts]**2/A[ts]
    verts = geometry['throat.offset_vertices']
    alpha = _sp.ones_like(C)
    for i in ts:
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
