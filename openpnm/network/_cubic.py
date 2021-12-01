import numpy as np
from openpnm.network import GenericNetwork
from openpnm import topotools
from openpnm.utils import logging, SettingsAttr
from auto_all import start_all, end_all
logger = logging.getLogger(__name__)


start_all()

class Cubic(GenericNetwork):
    r"""
    Simple cubic lattice with connectivity from 6 to 26

    Though simple, the Cubic network offers many advantages such as easy
    visualization and accurate determination of domain area and length in
    transport calculations.

    Parameters
    ----------
    shape : array_like
        The [Nx, Ny, Nz] size of the network in terms of the number of
        pores in each direction
    spacing : array_like, optional
        The spacing between pore centers in each direction. If not given,
        then [1, 1, 1] is assumed.
    connectivity : int, optional
        The number of connections to neighboring pores. Connections are
        made symmetrically to any combination of face, edge, or corners
        neighbors. The default is 6 to create a simple cubic structure,
        but options are:

            6
                Faces only
            14
                Faces and Corners
            18
                Faces and Edges
            20
                Edges and Corners
            26
                Faces, Edges and Corners

        For a more random distribution of connectivity, use a high
        connectivity (i.e. 26) and then delete a fraction of the throats
        using ``openpnm.topotools.reduce_coordination``.
    name : str
        An optional name for the object to help identify it. If not given,
        one will be generated.

    Examples
    --------
    .. plot::

       import openpnm as op
       import matplotlib.pyplot as plt

       pn = op.network.Cubic(shape=[5, 5, 5], spacing=[1, 1, 1])

       fig, ax = plt.subplots(figsize=(5, 5))
       op.topotools.plot_connections(network=pn, ax=ax)
       op.topotools.plot_coordinates(network=pn, c='r', s=75, ax=ax)

       plt.show()

    For larger networks and more control over presentation use `Paraview
    <http://www.paraview.org>`_.

    """

    def __init__(self, shape, spacing=[1, 1, 1], connectivity=6, **kwargs):
        super().__init__(**kwargs)

        # Take care of 1D/2D networks
        shape = np.array(shape, ndmin=1)
        shape = np.concatenate((shape, [1] * (3 - shape.size))).astype(int)

        arr = np.atleast_3d(np.empty(shape))

        spacing = np.float64(spacing)
        if spacing.size == 2:
            spacing = np.concatenate((spacing, [1]))
        spacing = np.ones(3, dtype=float) * np.array(spacing, ndmin=1)

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
        labels : str or list[str]
            The labels indicating the pores defining each face where
            boundary pores are to be added (e.g. 'left' or
            ['left', 'right'])
        spacing : scalar or array_like
            The spacing of the network (e.g. [1, 1, 1]). This should be
            given since it can be quite difficult to infer from the
            network, for instance if boundary pores have already added to
            other faces.

        """
        if isinstance(labels, str):
            labels = [labels]
        x, y, z = self["pore.coords"].T
        if spacing is None:
            spacing = topotools.get_spacing(self)
        else:
            spacing = np.array(spacing)
            if spacing.size == 1:
                spacing = np.ones(3) * spacing
        Lcx, Lcy, Lcz = spacing

        offset = {}
        shape = topotools.get_shape(self)
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

end_all()
