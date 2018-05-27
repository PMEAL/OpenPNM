from openpnm.network import GenericNetwork
from openpnm.topotools import generate_base_points, trim
import scipy as sp


class GenericTessellation(GenericNetwork):

    def _parse_points(self, shape, points, num_points):
        # Deal with input arguments
        if points is None:
            if num_points is None:
                raise Exception('Must specify either "points" or "num_points"')
            points = generate_base_points(num_points=num_points,
                                          domain_size=shape)

        # Deal with points that are only 2D...they break Delaunay
        if points.shape[1] == 3 and len(sp.unique(points[:, 2])) == 1:
            points = points[:, :2]

        return points

    def _find_external_pores(self, shape=None):
        r"""
        Identifies pores that lie outside the specified domain.

        Parameters
        ----------
        domain_size : array_like
            The size and shape of the domain beyond which points should be
            trimmed. The argument is treated as follows:

            **sphere** : If a scalar or single element list is received, it's
            treated as the radius [r] of a sphere centered on [0, 0, 0].

            **cylinder** : If a two-element list is received it's treated as
            the radius and height of a cylinder [r, z] whose central axis
            starts at [0, 0, 0] and extends in the positive z-direction.

            **rectangle** : If a three element list is received, it's treated
            as the outer corner of rectangle [x, y, z] whose opposite corner
            lies at [0, 0, 0].

        Returns
        -------
        An Np-long mask of with True values indicating pores that lie outside
        the domain.

        """
        # Label external pores for trimming below
        if len(shape) == 1:  # Spherical
            # Find external pores
            r = sp.sqrt(sp.sum(self['pore.coords']**2, axis=1))
            Ps = (r > shape)
        elif len(shape) == 2:  # Cylindrical
            # Find external pores outside radius
            r = sp.sqrt(sp.sum(self['pore.coords'][:, [0, 1]]**2, axis=1))
            Ps = (r > shape[0])
            # Find external pores above and below cylinder
            Ps1 = self['pore.coords'][:, 2] > shape[1]
            Ps2 = self['pore.coords'][:, 2] < 0
            Ps = Ps*(Ps1 + Ps2)
        elif len(shape) == 3:  # Rectilinear
            shape = sp.array(shape, dtype=float)
            try:
                lo_lim = shape[:, 0]
                hi_lim = shape[:, 1]
            except IndexError:
                lo_lim = sp.array([0, 0, 0])
                hi_lim = shape
            Ps1 = sp.any(self['pore.coords'] > hi_lim, axis=1)
            Ps2 = sp.any(self['pore.coords'] < lo_lim, axis=1)
            Ps = Ps1 + Ps2
        return Ps
