import numpy as np
from openpnm.network import GenericNetwork
from openpnm import topotools
from openpnm.network.generators import cubic, tools
from openpnm.utils import logging

logger = logging.getLogger(__name__)


class Cubic(GenericNetwork):
    r"""
    Simple cubic lattice with connectivity from 6 to 26

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

    """

    def __init__(self, shape, spacing=[1, 1, 1], connectivity=6,
                 name=None, project=None, **kwargs):
        super().__init__(name=name, project=project, **kwargs)
        temp = cubic(shape=shape,
                     spacing=spacing,
                     connectivity=connectivity)
        temp = tools.add_all_label(temp)
        self.update(temp)
        # Add some labels
        topotools.label_faces(network=self)
        Ps = self.pores('surface')
        Ts = self.find_neighbor_throats(pores=Ps, mode='xnor')
        self.set_label(throats=Ts, label='surface')
        self['throat.internal'] = True
        self['pore.internal'] = True

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
            spacing = tools.get_spacing(network=self)
        else:
            spacing = np.array(spacing)
            if spacing.size == 1:
                spacing = np.ones(3) * spacing
        Lcx, Lcy, Lcz = spacing

        offset = {}
        shape = tools.get_shape(self)
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
