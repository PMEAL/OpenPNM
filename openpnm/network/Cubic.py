# -*- coding: utf-8 -*-
"""
===============================================================================
Cubic: Generate lattice-like networks
===============================================================================

"""
import numpy as np
import scipy as sp
from openpnm.network import GenericNetwork
from openpnm import topotools


class Cubic(GenericNetwork):
    r"""
    This class generates a cubic network of the specified size and shape.

    Parameters
    ----------
    name : string
        A unique name for the network.  If none is given a random name will
        be assigned

    shape : tuple of ints
        The (i,j,k) size and shape of the network.

    connectivity : int
        The number of connections to neighboring pores.  Connections are made
        symmetrically to any combination of face, edge, or corners neighbors.

        Options are:

        - 6: Faces only
        - 8: Corners only
        - 12: Edges only
        - 14: Faces and Corners
        - 18: Faces and Edges
        - 20: Edges and Corners
        - 26: Faces, Edges and Corners

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[3,4,5])
    >>> pn.Np
    60

    """
    def __init__(self, shape, spacing=[1, 1, 1], connectivity=6, name=None,
                 project=None):

        arr = np.atleast_3d(np.empty(shape))

        # Store original network shape
        self._shape = sp.shape(arr)
        # Store network spacing
        self._spacing = sp.ones(3)*sp.array(spacing, ndmin=1)

        points = np.array([i for i, v in np.ndenumerate(arr)], dtype=float)
        points += 0.5
        points *= spacing

        I = np.arange(arr.size).reshape(arr.shape)

        face_joints = [(I[:, :, :-1], I[:, :, 1:]),
                       (I[:, :-1], I[:, 1:]),
                       (I[:-1], I[1:])]

        corner_joints = [(I[:-1, :-1, :-1], I[1:, 1:, 1:]),
                         (I[:-1, :-1, 1:], I[1:, 1:, :-1]),
                         (I[:-1, 1:, :-1], I[1:, :-1, 1:]),
                         (I[1:, :-1, :-1], I[:-1, 1:, 1:])]

        edge_joints = [(I[:, :-1, :-1], I[:, 1:, 1:]),
                       (I[:, :-1, 1:], I[:, 1:, :-1]),
                       (I[:-1, :, :-1], I[1:, :, 1:]),
                       (I[1:, :, :-1], I[:-1, :, 1:]),
                       (I[1:, 1:, :], I[:-1, :-1, :]),
                       (I[1:, :-1, :], I[:-1, 1:, :])]

        if connectivity == 6:
            joints = face_joints
        elif connectivity == 8:
            joints = corner_joints
        elif connectivity == 12:
            joints = edge_joints
        elif connectivity == 14:
            joints = face_joints + corner_joints
        elif connectivity == 18:
            joints = face_joints + edge_joints
        elif connectivity == 20:
            joints = edge_joints + corner_joints
        elif connectivity == 26:
            joints = face_joints + corner_joints + edge_joints
        else:
            raise Exception('Invalid connectivity receieved. Must be 6, 8, '
                            '12, 14, 18, 20 or 26')

        I = np.arange(arr.size).reshape(arr.shape)
        tails, heads = [], []
        for T, H in joints:
            tails.extend(T.flat)
            heads.extend(H.flat)
        pairs = np.vstack([tails, heads]).T

        self['pore.coords'] = points
        self['throat.conns'] = pairs

        super().__init__(Np=points.shape[0], Nt=pairs.shape[0], name=name,
                         project=project)

        self._label_surfaces()

    def _label_surfaces(self):
        r'''
        It applies the default surface labels for a cubic network
        '''
        x, y, z = self['pore.coords'].T
        labels = ['internal', 'front', 'back', 'left', 'right', 'bottom',
                  'top']
        for label in labels:
            if 'pore.'+label not in self.keys():
                self['pore.'+label] = False
        if 'pore.boundary' in self.keys():
            internal = -self['pore.boundary']
        else:
            internal = self['pore.all']
        self['pore.internal'] = internal
        self['pore.front'][x <= x.min()] = True
        self['pore.back'][x >= x.max()] = True
        self['pore.left'][y <= y.min()] = True
        self['pore.right'][y >= y.max()] = True
        self['pore.bottom'][z <= z.min()] = True
        self['pore.top'][z >= z.max()] = True

    def add_boundary_pores(self, labels=['top', 'bottom', 'front', 'back',
                                         'left', 'right']):
        r"""
        Add pores to the faces of the network for use as boundary pores.  Pores
        are offset from the faces by 1/2 a lattice spacing such that they lie
        directly on the boundaries.

        Parameters
        ----------
        labels : list of strings
            Controls which faces the boundary pores are added to.  Options
            include 'top', 'bottom', 'left', 'right', 'front', and 'back'.

        Notes
        -----
        This method uses ``clone_pores`` to clone the surface pores (labeled
        'left','right', etc), then shifts them to the periphery of the domain,
        and gives them the label 'right_face', 'left_face', etc.
        """
        if type(labels) == str:
            labels = [labels]
        x, y, z = self['pore.coords'].T
        Lcx, Lcy, Lcz = self._spacing

        offset = {}
        offset['front'] = offset['left'] = offset['bottom'] = [0, 0, 0]
        offset['back'] = [Lcx*self._shape[0], 0, 0]
        offset['right'] = [0, Lcy*self._shape[1], 0]
        offset['top'] = [0, 0, Lcz*self._shape[2]]

        scale = {}
        scale['front'] = scale['back'] = [0, 1, 1]
        scale['left'] = scale['right'] = [1, 0, 1]
        scale['bottom'] = scale['top'] = [1, 1, 0]

        for label in labels:
            Ps = self.pores(label)
            topotools.clone_pores(network=self, pores=Ps,
                                  labels=label+'_boundary')
            # Translate cloned pores
            ind = self.pores(label+'_boundary')
            coords = self['pore.coords'][ind]
            coords = coords*scale[label] + offset[label]
            self['pore.coords'][ind] = coords

    def _get_spacing(self):
        # Find Network spacing
        P1 = self['throat.conns'][:, 0]
        P2 = self['throat.conns'][:, 1]
        C1 = self['pore.coords'][P1]
        C2 = self['pore.coords'][P2]
        E = np.sqrt(np.sum((C1-C2)**2, axis=1))  # Euclidean distance
        if np.allclose(E, E[0]):
            spacing = E[0]
        else:
            raise Exception('A unique value of spacing could not be inferred')
        return spacing

    spacing = property(fget=_get_spacing)
