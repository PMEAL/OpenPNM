import numpy as np
from openpnm import topotools
from openpnm.network import Network
from openpnm.utils import Workspace
from openpnm.utils import Docorator
from openpnm._skgraph.generators import fcc


docstr = Docorator()
ws = Workspace()
__all__ = ['FaceCenteredCubic']


@docstr.dedent
class FaceCenteredCubic(Network):
    r"""
    Face-Centered Cubic lattice which is a simple cubic lattice with additional
    pores located on each face of the cubic unit cells. These pores are
    connected to each other and their diagonal neighbors.

    Parameters
    ----------
    shape : array_like
        The number of pores forming the corners of the unit cell in each
        direction
    spacing : array_like (optional)
        The spacing between pores that form the corners of the unit cells

    %(Network.parameters)s

    See Also
    --------
    Cubic
    BodyCenteredCubic

    Notes
    -----
    The pores are labelled as beloning to 'corner_sites' and 'face_sites'
    Throats are labelled by the which type of pores they connect, e.g.
    'throat.corner_to_body'.

    """

    def __init__(self, shape, mode='sc', spacing=1, **kwargs):
        super().__init__(**kwargs)
        shape = np.array(shape)
        if np.any(shape < 2):
            raise Exception('FCC lattice networks must have at least 2 '
                            'pores in all directions')
        net = fcc(shape=shape, spacing=spacing,
                  node_prefix='pore', edge_prefix='throat')
        self.update(net)
        # Add labels
        Ts1 = np.all(self['pore.corner'][self.conns], axis=1)
        self['throat.corner_to_corner'] = Ts1
        Ts2 = np.all(self['pore.face'][self.conns], axis=1)
        self['throat.face_to_face'] = Ts2
        self['throat.corner_to_face'] = ~(Ts1 + Ts2)

        topotools.label_faces(self)
        Ps = self.pores(['xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax'])
        Ps = self.to_mask(pores=Ps)
        self['pore.surface'] = Ps

        # Finally scale network to specified spacing
        self['pore.coords'] *= np.array(spacing)

    def add_boundary_pores(self, labels, spacing):
        r"""
        Add boundary pores to the specified faces of the network

        Pores are offset from the faces by 1/2 of the given ``spacing``,
        such that they lie directly on the boundaries.

        Parameters
        ----------
        labels : str or list[str]
            The labels indicating the pores defining each face where
            boundary pores are to be added (e.g. 'left' or
            ['left', 'right'])
        spacing : scalar or array_like
            The spacing of the network (e.g. [1, 1, 1]).  This must be
            given since it can be quite difficult to infer from the
            network, for instance if boundary pores have already added
            to other faces.

        """
        spacing = np.array(spacing)
        if spacing.size == 1:
            spacing = np.ones(3)*spacing
        for item in labels:
            Ps = self.pores(item)
            coords = np.absolute(self['pore.coords'][Ps])
            axis = np.count_nonzero(np.diff(coords, axis=0), axis=0) == 0
            offset = np.array(axis, dtype=int)/2
            if np.amin(coords) == np.amin(coords[:, np.where(axis)[0]]):
                offset = -1*offset
            topotools.add_boundary_pores(network=self, pores=Ps, offset=offset,
                                         apply_label=item + '_boundary')
