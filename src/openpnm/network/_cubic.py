import logging
import numpy as np
from openpnm.network import Network
from openpnm.utils import Docorator
from openpnm import topotools
from openpnm import _skgraph as skgr
skgr.settings.node_prefix = 'pore'
skgr.settings.edge_prefix = 'throat'

docstr = Docorator()
logger = logging.getLogger(__name__)
__all__ = ['Cubic']


@docstr.dedent
class Cubic(Network):
    r"""
    Simple cubic lattice with connectivity from 6 to 26

    Though simple, the Cubic network offers many advantages such as easy
    visualization and accurate determination of domain area and length in
    transport calculations.

    Parameters
    ----------
    shape : array_like
        The [Nx, Ny, Nz] size of the network in terms of the number of pores
        in each direction.  For a 2D network
    spacing : array_like, optional
        The spacing between pore centers in each direction. If not given,
        then [1, 1, 1] is assumed.
    connectivity : int, optional
        The number of connections to neighboring pores.  Connections are made
        symmetrically to any combination of face, edge, or corners neighbors.
        The default is 6 to create a simple cubic structure, but options are:

        ===========  =====================================================
        Value        Result
        ===========  =====================================================
        6            Faces only
        14           Faces and Corners
        18           Faces and Edges
        20           Edges and Corners
        26           Faces, Edges and Corners
        ===========  =====================================================

        For a more random distribution of connectivity, use a high
        ``connectivity`` (i.e. 26) and then delete a fraction of the throats
        using ``openpnm.topotools.reduce_coordination``.  Also note that
        corners-only and edges-only are not permitted since they create
        disconnected networks. If you require one of these topologies you
        can specify 14 or 18, then use ``openpnm.topotools.trim`` to remove
        the face-to-face connections, which can be identified by looking
        for throats with a length equal to the network spacing.

    %(Network.parameters)s

    """

    def __init__(self, shape, spacing=[1, 1, 1], connectivity=6, **kwargs):
        super().__init__(**kwargs)
        net = skgr.generators.cubic(shape=shape, spacing=spacing,
                                    connectivity=connectivity,
                                    node_prefix='pore', edge_prefix='throat')
        self.update(net)
        self._post_init()
        self["pore.surface"] = skgr.tools.find_surface_nodes_cubic(self)
        Ps = self["pore.surface"]
        self["throat.surface"] = np.all(Ps[self["throat.conns"]], axis=1)
        self.update(skgr.generators.tools.label_faces_cubic(self))

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
            spacing = skgr.tools.get_cubic_spacing(self)
        else:
            spacing = np.array(spacing)
            if spacing.size == 1:
                spacing = np.ones(3) * spacing
        Lcx, Lcy, Lcz = spacing

        offset = {}
        shape = skgr.tools.get_cubic_shape(self)
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
                                      labels=f"{label}_boundary",
                                      mode='parents')
                # Translate cloned pores
                ind = self.pores(label + "_boundary")
                coords = self["pore.coords"][ind]
                coords = coords * scale[label] + offset[label]
                self["pore.coords"][ind] = coords
            except KeyError:
                msg = f"No pores labelled {label} found, skipping boundary addition"
                logger.warning(msg)
