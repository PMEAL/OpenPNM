import scipy as sp
from scipy.spatial import ConvexHull
from transforms3d import _gohlketransforms as tr


def voronoi(target, **kwargs):
    r"""
    Update the throat normals from the voronoi vertices
    """
    verts = target["throat.vertices"]
    value = sp.ndarray([len(verts), 3])
    for i in range(len(verts)):
        if len(sp.unique(verts[i][:, 0])) == 1:
            verts_2d = sp.vstack((verts[i][:, 1], verts[i][:, 2])).T
        elif len(sp.unique(verts[i][:, 1])) == 1:
            verts_2d = sp.vstack((verts[i][:, 0], verts[i][:, 2])).T
        else:
            verts_2d = sp.vstack((verts[i][:, 0], verts[i][:, 1])).T
        hull = ConvexHull(verts_2d, qhull_options='QJ Pp')
        sorted_verts = verts[i][hull.vertices]
        v1 = sorted_verts[1]-sorted_verts[0]
        v2 = sorted_verts[-1]-sorted_verts[0]
        value[i] = sp.cross(v1, v2)

    return value


def pore_coords(target, **kwargs):
    r"""
    Unit vector from P1 to P2 as defined in throat.conns
    """
    network = target.project.network
    conns = network['throat.conns']
    P1 = conns[:, 0]
    P2 = conns[:, 1]
    coords = network['pore.coords']
    vec = coords[P2] - coords[P1]
    unit_vec = tr.unit_vector(vec, axis=1)
    return unit_vec
