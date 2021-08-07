"""
===============================================================================
Cubic: Generate lattice-like networks
===============================================================================

"""
import numpy as np
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
    >>> import matplotlib.pyplot as plt
    >>> pn = op.network.Cubic(shape=[5, 5, 5], spacing=[1, 1, 1])
    >>> pn.Np
    125

    And it can be plotted for quick visualization using:

    >>> fig, ax = plt.subplots()
    >>> _ = op.topotools.plot_connections(network=pn, ax=ax)
    >>> _ = op.topotools.plot_coordinates(network=pn, c='r', s=75, ax=ax)

    .. image:: /../docs/_static/images/cubic_network.png
        :align: center

    For larger networks and more control over presentation use `Paraview
    <http://www.paraview.org>`_.
    """

    def __init__(self, shape, spacing=[1, 1, 1], connectivity=6,
                 name=None, project=None, **kwargs):

        super().__init__(name=name, project=project, **kwargs)
        logger.critical('front and back labels have been switched to obey ' +
                        'the right-hand rule')

        # Take care of 1D/2D networks
        shape = np.array(shape, ndmin=1)
        shape = np.concatenate((shape, [1] * (3 - shape.size))).astype(int)

        arr = np.atleast_3d(np.empty(shape))

        # Store original network shape
        self.settings['shape'] = np.shape(arr)
        # Store network spacing
        spacing = np.float64(spacing)
        if spacing.size == 2:
            spacing = np.concatenate((spacing, [1]))
        spacing = np.ones(3, dtype=float) * np.array(spacing, ndmin=1)
        self.settings['spacing'] = spacing.tolist()

        z = np.tile(np.arange(shape[2]), shape[0] * shape[1])
        y = np.tile(np.repeat(np.arange(shape[1]), shape[2]), shape[0])
        x = np.repeat(np.arange(shape[0]), shape[1] * shape[2])
        points = (np.vstack([x, y, z]).T).astype(float) + 0.5

        idx = np.arange(arr.size).reshape(arr.shape)

        face_joints = [(idx[:, :, :-1], idx[:, :, 1:]),
                       (idx[:, :-1], idx[:, 1:]),
                       (idx[:-1], idx[1:])]

        corner_joints = [(idx[:-1, :-1, :-1], idx[1:, 1:, 1:]),
                         (idx[:-1, :-1, 1:], idx[1:, 1:, :-1]),
                         (idx[:-1, 1:, :-1], idx[1:, :-1, 1:]),
                         (idx[1:, :-1, :-1], idx[:-1, 1:, 1:])]

        edge_joints = [(idx[:, :-1, :-1], idx[:, 1:, 1:]),
                       (idx[:, :-1, 1:], idx[:, 1:, :-1]),
                       (idx[:-1, :, :-1], idx[1:, :, 1:]),
                       (idx[1:, :, :-1], idx[:-1, :, 1:]),
                       (idx[1:, 1:, :], idx[:-1, :-1, :]),
                       (idx[1:, :-1, :], idx[:-1, 1:, :])]

        if connectivity == 6:
            joints = face_joints
        elif connectivity == 6 + 8:
            joints = face_joints + corner_joints
        elif connectivity == 6 + 12:
            joints = face_joints + edge_joints
        elif connectivity == 12 + 8:
            joints = edge_joints + corner_joints
        elif connectivity == 6 + 8 + 12:
            joints = face_joints + corner_joints + edge_joints
        else:
            raise Exception("Invalid connectivity. Must be 6, 14, 18, 20 or 26.")

        tails, heads = np.array([], dtype=int), np.array([], dtype=int)
        for T, H in joints:
            tails = np.concatenate((tails, T.flatten()))
            heads = np.concatenate((heads, H.flatten()))
        pairs = np.vstack([tails, heads]).T

        self["pore.all"] = np.ones([points.shape[0], ], dtype=bool)
        self["throat.all"] = np.ones([pairs.shape[0], ], dtype=bool)
        self["pore.coords"] = points
        self["throat.conns"] = pairs
        self["pore.internal"] = True
        self["throat.internal"] = True
        self._label_surface_pores()
        topotools.label_faces(network=self)
        Ps = self["pore.surface"]
        self["throat.surface"] = np.all(Ps[self["throat.conns"]], axis=1)
        # Scale network to requested spacing
        self["pore.coords"] *= spacing

    def _label_surface_pores(self):
        r"""
        """
        hits = np.zeros_like(self.Ps, dtype=bool)
        dims = topotools.dimensionality(self)
        mn = np.amin(self["pore.coords"], axis=0)
        mx = np.amax(self["pore.coords"], axis=0)
        for ax in [0, 1, 2]:
            if dims[ax]:
                hits += self["pore.coords"][:, ax] <= mn[ax]
                hits += self["pore.coords"][:, ax] >= mx[ax]
        self["pore.surface"] = hits

    def add_boundary_pores(self, labels=["top", "bottom", "front",
                                         "back", "left", "right"],
                           spacing=None):
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
        if isinstance(labels, str):
            labels = [labels]
        x, y, z = self["pore.coords"].T
        if spacing is None:
            spacing = self._get_spacing()
        else:
            spacing = np.array(spacing)
            if spacing.size == 1:
                spacing = np.ones(3) * spacing
        Lcx, Lcy, Lcz = spacing

        offset = {}
        shape = self.settings['shape']
        offset["front"] = offset["left"] = offset["bottom"] = [0, 0, 0]
        offset["right"] = [Lcx * shape[0], 0, 0]
        offset["back"] = [0, Lcy * shape[1], 0]
        offset["top"] = [0, 0, Lcz * shape[2]]

        scale = {}
        scale["left"] = scale["right"] = [0, 1, 1]
        scale["front"] = scale["back"] = [1, 0, 1]
        scale["bottom"] = scale["top"] = [1, 1, 0]

        for label in labels:
            try:
                Ps = self.pores(label)
                topotools.clone_pores(network=self, pores=Ps,
                                      labels=label + "_boundary",
                                      mode='parents')
                # Translate cloned pores
                ind = self.pores(label + "_boundary")
                coords = self["pore.coords"][ind]
                coords = coords * scale[label] + offset[label]
                self["pore.coords"][ind] = coords
            except KeyError:
                logger.warning("No pores labelled " + label
                               + " were found, skipping boundary addition")

    def _get_spacing(self):
        # Find Network spacing
        P12 = self["throat.conns"]
        C12 = self["pore.coords"][P12]
        mag = np.linalg.norm(np.diff(C12, axis=1), axis=2)
        unit_vec = np.around(np.squeeze(np.diff(C12, axis=1)) / mag, decimals=14)
        spacing = [0, 0, 0]
        dims = topotools.dimensionality(self)
        # Ensure vectors point in n-dims unique directions
        c = {tuple(row): 1 for row in unit_vec}
        mag = np.atleast_1d(mag.squeeze()).astype(float)
        if len(c.keys()) > sum(dims):
            raise Exception(
                "Spacing is undefined when throats point in more directions"
                " than network has dimensions."
            )
        for ax in [0, 1, 2]:
            if dims[ax]:
                inds = np.where(unit_vec[:, ax] == unit_vec[:, ax].max())[0]
                temp = np.unique(mag[inds])
                if not np.allclose(temp, temp[0]):
                    raise Exception("A unique value of spacing could not be found.")
                spacing[ax] = temp[0]
        self.settings['spacing'] = spacing
        return np.array(spacing)

    spacing = property(fget=_get_spacing)

    def _get_shape(self):
        L = np.ptp(self["pore.coords"], axis=0)
        mask = L.astype(bool)
        S = self.spacing
        shape = np.array([1, 1, 1], int)
        shape[mask] = L[mask] / S[mask] + 1
        self.settings['shape'] = shape.tolist()
        return shape

    shape = property(fget=_get_shape)

    @property
    def _shape(self):
        return self.settings['shape']

    @property
    def _spacing(self):
        return np.array(self.settings['spacing'])

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
        if np.shape(values)[0] > self.num_pores("internal"):
            raise Exception("The array shape does not match the network")
        Ps = np.array(self["pore.index"][self.pores("internal")], dtype=int)
        arr = np.ones(self.settings['shape']) * np.nan
        ind = np.unravel_index(Ps, self.settings['shape'])
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
        array = np.atleast_3d(array)
        if np.shape(array) != self._shape:
            raise Exception("The array shape does not match the network")
        temp = array.flatten()
        Ps = np.array(self["pore.index"][self.pores("internal")], dtype=int)
        propname = "pore." + propname.split(".")[-1]
        self[propname] = np.nan
        self[propname][self.pores("internal")] = temp[Ps]
