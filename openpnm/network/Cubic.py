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
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class Cubic(GenericNetwork):
    r"""
    Simple cubic lattice with connectivity from 6 to 26

    Though simple, the Cubic network offers many advantages such as easy
    visualization and accurate determination of domain area and length in
    transport calculations.

    Parameters
    ----------
    shape : array_like
        The [Nx, Ny, Nz] size of the network in terms of the number of pores in
        each direction

    spacing : array_like, optional
        The spacing between pore centers in each direction. If not given, then
        [1, 1, 1] is assumed.

    connectivity : int, optional
        The number of connections to neighboring pores.  Connections are made
        symmetrically to any combination of face, edge, or corners neighbors.
        The default is 6 to create a simple cubic structure, but options are:

        - 6: Faces only
        - 8: Corners only
        - 12: Edges only
        - 14: Faces and Corners
        - 18: Faces and Edges
        - 20: Edges and Corners
        - 26: Faces, Edges and Corners

        For a more random distribution of connectivity, use a high
        ``connectivity`` (i.e. 26) and then delete a fraction of the throats
        using ``openpnm.topotools.reduce_coordination``.

    name : string
        An optional name for the object to help identify it.  If not given,
        one will be generated.

    project : OpenPNM Project object, optional
        Each OpenPNM object must be part of a *Project*.  If none is supplied
        then one will be created and this Network will be automatically
        assigned to it.  To create a *Project* use ``openpnm.Project()``.

    Attributes
    ----------
    spacing : int or array
        The distance between pore centers.  This value becomes meaningless
        if the topology is manipulated at all (i.e. by adding boundary pores)
        since there is not unique or consistent value.  In such cases an
        exception is thrown.

    shape : array
        The shape of the network.  Like ``spacing`` this values is meaningless
        if the topology is manipulated, so an Exception is thrown.

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[5, 5, 5], spacing=[1, 1, 1])
    >>> pn.Np
    125

    And it can be plotted for quick visualization using:

    >>> fig = op.topotools.plot_connections(network=pn)
    >>> fig = op.topotools.plot_coordinates(network=pn, c='r', s=75, fig=fig)

    .. image:: /../docs/static/images/cubic_network.png
        :align: center

    For larger networks and more control over presentation use `Paraview
    <http://www.paraview.org>`_.
    """
    def __init__(self, shape, spacing=[1, 1, 1], connectivity=6, name=None,
                 project=None):

        arr = np.atleast_3d(np.empty(shape))

        # Store original network shape
        self._shape = sp.shape(arr)
        # Store network spacing
        spacing = sp.around(spacing, decimals=15)
        if spacing.size == 2:
            spacing = sp.concatenate((spacing, [1]))
        self._spacing = sp.ones(3)*sp.array(spacing, ndmin=1)

        points = np.array([i for i, v in np.ndenumerate(arr)], dtype=float)
        points += 0.5

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

        super().__init__(Np=points.shape[0], Nt=pairs.shape[0], name=name,
                         project=project)

        self['pore.coords'] = points
        self['throat.conns'] = pairs
        self['pore.internal'] = True
        self['throat.internal'] = True
        self._label_surface_pores()
        topotools.label_faces(network=self)
        Ps = self['pore.surface']
        self['throat.surface'] = sp.all(Ps[self['throat.conns']], axis=1)
        # Scale network to requested spacing
        self['pore.coords'] *= spacing

    def _label_surface_pores(self):
        r"""
        """
        hits = sp.zeros_like(self.Ps, dtype=bool)
        dims = topotools.dimensionality(self)
        mn = sp.amin(self['pore.coords'], axis=0)
        mx = sp.amax(self['pore.coords'], axis=0)
        for ax in [0, 1, 2]:
            if dims[ax]:
                hits += self['pore.coords'][:, ax] <= mn[ax]
                hits += self['pore.coords'][:, ax] >= mx[ax]
        self['pore.surface'] = hits

    def add_boundary_pores(self, labels=['top', 'bottom', 'front', 'back',
                                         'left', 'right'], spacing=None):
        r"""
        Add pores to the faces of the network for use as boundary pores.

        Pores are offset from the faces by 1/2 a lattice spacing such that
        they lie directly on the boundaries.

        Parameters
        ----------
        labels : string or list of strings
            The labels indicating the pores defining each face where boundary
            pores are to be added (e.g. 'left' or ['left', 'right'])

        spacing : scalar or array_like
            The spacing of the network (e.g. [1, 1, 1]).  This should be given
            since it can be quite difficult to infer from the network, for
            instance if boundary pores have already added to other faces.

        """
        if type(labels) == str:
            labels = [labels]
        x, y, z = self['pore.coords'].T
        if spacing is None:
            spacing = self._get_spacing()
        else:
            spacing = sp.array(spacing)
            if spacing.size == 1:
                spacing = sp.ones(3)*spacing
        Lcx, Lcy, Lcz = spacing

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
            try:
                Ps = self.pores(label)
                topotools.clone_pores(network=self, pores=Ps,
                                      labels=label+'_boundary')
                # Translate cloned pores
                ind = self.pores(label+'_boundary')
                coords = self['pore.coords'][ind]
                coords = coords*scale[label] + offset[label]
                self['pore.coords'][ind] = coords
            except KeyError:
                logger.warning('No pores labelled ' + label + ' were found, ' +
                               'skipping boundary addition')

    def _get_spacing(self):
        # Find Network spacing
        P12 = self['throat.conns']
        C12 = self['pore.coords'][P12]
        mag = np.sqrt(np.sum(np.diff(C12, axis=1)**2, axis=2))
        vec = sp.around(sp.squeeze(np.diff(C12, axis=1))/mag,
                        decimals=10)
        spacing = [0, 0, 0]
        dims = topotools.dimensionality(self)
        # Ensure vectors point in n-dims unique directions
        c = {tuple(row): 1 for row in vec}
        if len(c.keys()) > sum(dims):
            raise Exception('Spacing is undefined when throats point in ' +
                            'more directions than network has dimensions')
        mag = sp.around(mag.squeeze(), decimals=10)
        for ax in [0, 1, 2]:
            if dims[ax]:
                inds = sp.where(vec[:, ax] == vec[:, ax].max())[0]
                temp = sp.unique(mag[inds])
                if np.size(temp) > 1:
                    raise Exception('A unique value of spacing could not be found')
                else:
                    spacing[ax] = temp[0]
        return sp.array(spacing)

    spacing = property(fget=_get_spacing)

    def _get_shape(self):
        L = np.ptp(self['pore.coords'], axis=0)
        mask = L.astype(bool)
        S = self.spacing
        shape = sp.array([1, 1, 1], int)
        shape[mask] = L[mask]/S[mask] + 1
        return shape

    shape = property(fget=_get_shape)

    def to_array(self, values):
        r"""
        Converts the values to a rectangular array with the same shape as the
        network

        Parameters
        ----------
        values : array_like
            An Np-long array of values to convert to

        Notes
        -----
        This method can break on networks that have had boundaries added.  It
        will usually work IF the given values came only from 'internal'
        pores.

        """
        if sp.shape(values)[0] > self.num_pores('internal'):
            raise Exception('The array shape does not match the network')
        Ps = sp.array(self['pore.index'][self.pores('internal')], dtype=int)
        arr = sp.ones(self._shape)*sp.nan
        ind = sp.unravel_index(Ps, self._shape)
        arr[ind[0], ind[1], ind[2]] = values
        return arr

    def from_array(self, array, propname):
        r"""
        Apply data to the network based on a rectangular array filled with
        values.  Each array location corresponds to a pore in the network.

        Parameters
        ----------
        array : array_like
            The rectangular array containing the values to be added to the
            network. This array must be the same shape as the original network.

        propname : string
            The name of the pore property being added.

        """
        array = sp.atleast_3d(array)
        if sp.shape(array) != self._shape:
            raise Exception('The array shape does not match the network')
        temp = array.flatten()
        Ps = sp.array(self['pore.index'][self.pores('internal')], dtype=int)
        propname = 'pore.' + propname.split('.')[-1]
        self[propname] = sp.nan
        self[propname][self.pores('internal')] = temp[Ps]
