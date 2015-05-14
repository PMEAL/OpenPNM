"""
===============================================================================
DelaunayCubic: Generate semi-random networks based on Delaunay Tessellations and
perturbed cubic lattices
===============================================================================

"""
import OpenPNM
import scipy as sp
import sys
import numpy as np
from OpenPNM.Network.__Delaunay__ import Delaunay
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)


class DelaunayCubic(Delaunay):
    r"""
    This class contains the methods for creating a *Delaunay* network topology
    based connecting pores with a Delaunay tessellation.

    This Subclass of Delaunay generates points on a cubic lattice and then perturbs
    them to prevent degeneracy

    Parameters
    ----------
    name : string
        A unique name for the network

    Examples
    --------
    >>> import OpenPNM
    >>> pn = OpenPNM.Network.DelaunayCubic(shape=[5, 5, 5],
    ...                                    spacing=[4e-5, 4e-5, 4e-5],
    ...                                    jiggle_factor=0.01)
    >>> pn.num_pores()
    125

    """

    def __init__(self, shape=None, template=None, spacing=[1, 1, 1],
                 jiggle_factor=0.1, arrangement='SC', **kwargs):
        if shape is not None:
            self._arr = np.atleast_3d(np.empty(shape))
        elif template is not None:
            self._arr = sp.array(template, ndmin=3, dtype=bool)
        else:
            self._arr = np.atleast_3d(np.empty([3, 3, 3]))

        # Store original network shape
        self._shape = sp.shape(self._arr)
        # Store network spacing instead of calculating it
        self._spacing = sp.asarray(spacing)
        self._num_pores = np.prod(np.asarray(self._shape))
        self._domain_size = np.asarray(self._shape) * self._spacing
        self._jiggle_factor = jiggle_factor
        self._arrangement = arrangement
        super().__init__(num_pores=self._num_pores,
                         domain_size=self._domain_size,
                         **kwargs)

    def _generate_pores(self):
        r"""
        Generate the pores with numbering scheme.
        """

        points = np.array([i for i, v in np.ndenumerate(self._arr)], dtype=float)
        points += 0.5

        # 2D Orthorhombic adjustment - shift even rows back a bit and odd rows
        # forward a bit
        "   0   0   0   "
        " 0   0   0   0 "
        "   0   0   0   "
        if self._arrangement == 'O':
            shift_y = np.array([0, 0.25, 0])
            shift_x = np.array([0.25, 0, 0])
            points[(points[:, 0] % 2 == 0)] -= shift_y
            points[(points[:, 2] % 2 != 0)] -= shift_x
            points[(points[:, 0] % 2 != 0)] += shift_y
            points[(points[:, 2] % 2 == 0)] += shift_x
        # BCC = Body Centre Cubic
        if self._arrangement == 'BCC':
            body_points = []
            for i in range(1, self._shape[0]):
                for j in range(1, self._shape[1]):
                    for k in range(1, self._shape[2]):
                        body_points.append([i, j, k])
            body_points = np.asarray(body_points)
            points = np.concatenate((points, body_points))

        # FCC = Face Centre Cubic
        if self._arrangement == 'FCC':
            face_points = []
            for i in range(1, self._shape[0]):
                for j in range(1, self._shape[1]):
                    for k in range(1, self._shape[2]):
                        left = [i-0.5, j, k]
                        right = [i+0.5, j, k]
                        back = [i, j-0.5, k]
                        front = [i, j+0.5, k]
                        bottom = [i, j, k-0.5]
                        top = [i, j, k+0.5]
                        if left not in face_points:
                            face_points.append(left)
                        if right not in face_points:
                            face_points.append(right)
                        if back not in face_points:
                            face_points.append(back)
                        if front not in face_points:
                            face_points.append(front)
                        if bottom not in face_points:
                            face_points.append(bottom)
                        if top not in face_points:
                            face_points.append(top)
            face_points = np.asarray(face_points)
            points = np.concatenate((points, face_points))

        jiggle = (np.random.rand(len(points), 3)-0.5)*self._jiggle_factor
        points += jiggle
        points *= self._spacing

        self['pore.coords'] = points
        logger.debug(sys._getframe().f_code.co_name + ': End of method')
