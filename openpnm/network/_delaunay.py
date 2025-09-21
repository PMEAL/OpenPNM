from openpnm.network import Network
from openpnm.utils import Docorator
from openpnm._skgraph.generators import delaunay


__all__ = ['Delaunay']
docstr = Docorator()


@docstr.dedent
class Delaunay(Network):
    r"""
    Random network formed by Delaunay tessellation of arbitrary base points

    Parameters
    ----------
    points : array_like or int
        Can either be an N-by-D array of point coordinates which will be used,
        or a scalar value indicating the number of points to generate. If 3D points
        are supplied, and a 2D shape is specified, then the z-coordinate is ignored.
    shape : array_like
        The size and shape of the domain:

        ========== ============================================================
        shape      result
        ========== ============================================================
        [x, y, z]  A 3D cubic domain of dimension x, y and z
        [x, y, 0]  A 2D square domain of size x by y
        ========== ============================================================

    reflect : bool, optional
        If ``True`` then the base points will be reflected across
        all the faces of the domain prior to performing the tessellation. This
        feature is best combined with ``trim=True`` to prevent unreasonably long
        connections between points on the surfaces.
    f : float
        The fraction of points which should be reflected.  The default is 1 which
        reflects all the points in the domain, but this can lead to a lot of
        unnecessary points, so setting to 0.1 or 0.2 helps speed, but risks that
        the tessellation may not have smooth faces if not enough points are
        reflected.
    trim : bool, optional
        If ``True`` (default) then all vertices laying outside the domain will
        be removed. This is only useful if ``reflect=True``.

    %(Network.parameters)s

    Notes
    -----
    It is also possible to generate circular ``[r, 0]``, cylindrical ``[r, z]``, and
    spherical domains ``[r]``, but this feature does not quite work as desired.
    It does not produce a truly clean outer surface since the tessellation are
    conducted in cartesian coordinates, so the circular and spherical surfaces
    have artifacts. Scipy recently added the ability to do tessellations on
    spherical surfaces, for geological applications, but this is not flexible
    enough, yet.

    """

    def __init__(self, shape, points, reflect=True, trim=True, f=1, **kwargs):
        super().__init__(**kwargs)
        net, tri = delaunay(
            points=points,
            shape=shape,
            reflect=reflect,
            f=f,
            trim=trim,
            node_prefix='pore',
            edge_prefix='throat',
        )
        self.update(net)
        self._post_init()
        self.tri = tri
