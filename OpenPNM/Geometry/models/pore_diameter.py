r"""
===============================================================================
pore_diameter
===============================================================================

"""
import scipy as _sp
import numpy as np

from scipy.spatial import ConvexHull


def sphere(geometry, psd_name, psd_shape, psd_loc, psd_scale,
           pore_seed='pore.seed', psd_offset=0, **kwargs):
    r"""
    Calculate pore diameter from given seed values.

    Parameters
    ----------
    geometry : OpenPNM Geometry Object
        The Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    psd_name : string
        The name of the statistical distribution to use. This model uses the
        Scipy.stats module, so any of the distributions available there are
        suitable options.

    psd_shape, loc and scale : scalars
        The parameters to send to the selected statistics model.  Most of the
        Scipy.stats models accept the same keyword arguments.  Note that the
        psd_ prefix is added by OpenPNM to indicate 'pore size distribution'.

    psd_offset : scalar
        Controls the minimum value in the pore size distribution by shifting
        the entire set of values by the given offset.  This is useful for
        avoiding pore sizes too close to zero.

    """
    import scipy.stats as spst
    prob_fn = getattr(spst, psd_name)
    P = prob_fn(psd_shape, loc=psd_loc, scale=psd_scale)
    value = P.ppf(geometry[pore_seed]) + psd_offset
    return value


def equivalent_sphere(geometry, pore_volume='pore.volume', **kwargs):
    r"""
    Calculate pore diameter as the diameter of a sphere with an equivalent
    volume.

    Parameters
    ----------
    geometry : OpenPNM Geometry Object
        The Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    pore_volume : string
        The dictionary key containing the pore volume values
    """
    from scipy.special import cbrt
    pore_vols = geometry[pore_volume]
    value = cbrt(6*pore_vols/_sp.pi)
    return value


def equivalent_cube(geometry, pore_volume='pore.volume', **kwargs):
    r"""
    Calculate pore diameter as the width of a cube with an equivalent volume.

    Parameters
    ----------
    geometry : OpenPNM Geometry Object
        The Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    pore_volume : string
        The dictionary key containing the pore volume values
    """
    from scipy.special import cbrt
    pore_vols = geometry[pore_volume]
    value = cbrt(pore_vols)
    return value


def voronoi(geometry, pore_volume='pore.volume', **kwargs):
    r"""
    Calculate pore diameter from equivalent sphere - volumes must be calculated
    first.
    """
    from scipy.special import cbrt
    pore_vols = geometry[pore_volume]
    value = cbrt(6*pore_vols/_sp.pi)
    return value


def centroids(network, geometry, **kwargs):
    r"""
    Calculate the diameter representing an inclosed sphere. The maximum is very
    difficult to caluclate for irregular polygons with more than 4 faces so an
    average distance from the pore centroid to the throat centroid is an
    approximation
    """
    Np = geometry.num_pores()
    value = _sp.zeros(Np)
    pore_map = geometry.map_pores(geometry.pores(), geometry._net)
    for geom_pore, net_pore in pore_map:
        net_throats = geometry._net.find_neighbor_throats(net_pore)
        geom_throats = geometry._net.map_throats(net_throats, geometry)[:, 1]
        tcs = geometry["throat.centroid"][geom_throats]
        pc = geometry["pore.centroid"][geom_pore]
        value[geom_pore] = _sp.mean(_sp.sqrt(((tcs - pc) * (tcs - pc))[:, 0] +
                                             ((tcs-pc)*(tcs-pc))[:, 1] +
                                             ((tcs-pc)*(tcs-pc))[:, 2]))*2
    return value


def insphere(network, geometry, **kwargs):
    r"""
    Calculate the diameter of the insphere using linear programming
    """
    import warnings
    try:
        import pulp as pu
        Np = geometry.num_pores()
        value = _sp.zeros(Np)
        pore_map = geometry.map_pores(geometry.pores(), geometry._net)
        for geom_pore, net_pore in pore_map:
            net_throats = geometry._net.find_neighbor_throats(net_pore)
            geom_throats = geometry._net.map_throats(net_throats, geometry)[:, 1]
            verts = geometry['throat.offset_vertices'][geom_throats]
            if len(verts) > 1:
                try:
                    pts = np.vstack((i for i in verts if len(i) > 0))
                except ValueError:
                    pts = []
                if len(pts) > 4:
                    "Work out central point to use as initial guess"
                    c0 = np.mean(pts, axis=0)
                    "Compute convex hull to find points lying on the hull in order"
                    hull = ConvexHull(pts, qhull_options='QJ Pp')
                    "For each simplex making up the hull collect the end points"
                    A = pts[hull.simplices[:, 0]]
                    B = pts[hull.simplices[:, 1]]
                    C = pts[hull.simplices[:, 2]]
                    "Normal of the simplices"
                    N = np.cross((B-A), (C-A), axis=1)
                    "Normalize the normal vector"
                    L = np.linalg.norm(N, axis=1)
                    F = np.vstack((L, L, L)).T
                    N /= F
                    "If normals point out of hull change sign to point in"
                    pointing_out = (np.sum((A-c0)*N, axis=1) > 0)
                    N[pointing_out] *= -1
                    "Define Linear Program Variables"
                    "The centre of the incircle adjustment"
                    cx = pu.LpVariable('cx', None, None, pu.LpContinuous)
                    cy = pu.LpVariable('cy', None, None, pu.LpContinuous)
                    cz = pu.LpVariable('cy', None, None, pu.LpContinuous)
                    "Radius of the incircle"
                    R = pu.LpVariable('R', 0, None, pu.LpContinuous)
                    # Slack variables for shortest distance between centre
                    # and simplices
                    S = pu.LpVariable.dict('SlackVariable', range(len(A)), 0,
                                           None, pu.LpContinuous)
                    "Set up LP problem"
                    prob = pu.LpProblem('FindInRadius', pu.LpMaximize)
                    "Objective Function"
                    prob += R
                    for i in range(len(A)):
                        "Ni.(C-Ai)-Si = 0"
                        prob += N[i][0]*(c0[0] + cx) + N[i][1]*(c0[1] + cy) + \
                            N[i][2]*(c0[2] + cz) - N[i][0]*A[i][0] - \
                            N[i][1]*A[i][1] - N[i][2]*A[i][2] - S[i] == 0
                        "Si >= R"
                        prob += S[i] >= R
                    "Solve the LP"
                    with warnings.catch_warnings():
                        warnings.simplefilter('ignore')
                        prob.solve()
                    # As the radius is the objective function we can get it from the
                    # objective or as R.value()
                    rad = prob.objective.value()
                    value[geom_pore] = rad*2

        return value

    except ImportError:
        print('Cannot use insphere method without installing pulp package')
