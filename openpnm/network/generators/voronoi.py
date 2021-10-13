import scipy.spatial as sptl
import numpy as np
from openpnm.topotools import vor_to_am, isoutside
from openpnm.network.generators import tools


def voronoi(points, shape=[1, 1, 1], crop=False):
    r"""
    Generate a network based on a Voronoi tessellation of base points

    Parameters
    ----------
    points : array_like or int
        Can either be an N-by-3 array of point coordinates which will be used,
        or a scalar value indicating the number of points to generate.
    shape : array_like
        Indicates the size and shape of the domain.
    crop : boolean (default = ``False``)
        If ``True`` (default) any points in the tessellation that lie outside
        the domain (defined by ``shape``) are removed, along with any throats
        which were connected to them

    Returns
    -------
    network : dict
        A dictionary containing ``'pore.coords'`` and ``'throat.conns'``
    vor : Voronoi tessellation object
        The Voronoi tessellation object produced by ``scipy.spatial.Voronoi``

    """
    points = tools.parse_points(points=points, shape=shape)
    mask = ~np.all(points == 0, axis=0)
    # Perform tessellation
    vor = sptl.Voronoi(points=points[:, mask])
    # Convert to adjecency matrix
    coo = vor_to_am(vor)
    # Write values to dictionary
    d = {}
    conns = np.vstack((coo.row, coo.col)).T
    d['throat.conns'] = conns
    if np.any(mask == False):
        verts = np.zeros([np.shape(vor.vertices)[0], 3])
        for i, col in enumerate(np.where(mask)[0]):
            verts[:, col] = vor.vertices[:, i]
    else:
        verts = np.copy(vor.vertices)
    d['pore.coords'] = verts

    if crop:
        Ps = isoutside(verts, shape=shape)
        d = tools.trim(network=d, pores=np.where(Ps)[0])

    return d, vor


if __name__ == "__main__":
    vn, vor = voronoi(points=50, shape=[1, 0, 1])
    print(vn.keys())
    print(vn['pore.coords'].shape)
    print(vn['throat.conns'].shape)
    vn, vor = voronoi(points=50, shape=[1, 0, 1], crop=True)
    print(vn.keys())
    print(vn['pore.coords'].shape)
    print(vn['throat.conns'].shape)
