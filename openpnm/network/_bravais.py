import logging
import numpy as np
from openpnm import topotools
from openpnm.network import GenericNetwork
from openpnm.utils import Workspace
from openpnm.utils import Docorator
from openpnm._skgraph.generators import cubic, fcc, bcc


docstr = Docorator()
logger = logging.getLogger(__name__)
ws = Workspace()
__all__ = ['Bravais']


@docstr.dedent
class Bravais(GenericNetwork):
    r"""
    Crystal lattice types including fcc, bcc, sc, and hcp

    These arrangements not only allow more dense packing than the standard
    Cubic for higher porosity materials, but also have more interesting
    non-straight connections between the various pore sites.

    More information on Bravais lattice notation can be `found on
    wikipedia <https://en.wikipedia.org/wiki/Bravais_lattice>`_.

    Parameters
    ----------
    shape : array_like
        The number of pores in each direction. This value is a bit
        ambiguous for the more complex unit cells used here, but generally
        refers to the the number for 'corner' sites
    spacing : array_like (optional)
        The spacing between pores in all three directions. Like the
        ``shape`` this is a bit ambiguous but refers to the spacing
        between corner sites. Essentially it controls the dimensions of
        the unit cell. If a scalar is given it is applied to all
        directions. The default is 1.
    mode : str
        The type of lattice to create. Options are:
            ===========  =====================================================
            mode         meaning
            ===========  =====================================================
            'sc'         Simple cubic (Same as ``Cubic``)
            'bcc'        body-centered cubic lattice
            'fcc'        Face-centered cubic lattice
            'hcp'        Hexagonal close packed (Not Implemented Yet)
            ===========  =====================================================

    %(GenericNetwork.parameters)s

    See Also
    --------
    Cubic
    CubicDual

    Notes
    -----
    The pores are labelled as beloning to 'corner_sites' and 'body_sites' in
    bcc or 'face_sites' in fcc.  Throats are labelled by the which type of
    pores they connect, e.g. 'throat.corner_to_body'.

    Limitations:

    * In principle Bravais lattice can also have a skew to them, but this is
      not implemented yet.
    * Support for 2D networks has not been added yet.
    * Hexagonal Close Packed (hcp) has not been implemented yet, but is on
      the todo list.

    """

    def __init__(self, shape, mode='sc', spacing=1, **kwargs):
        super().__init__(**kwargs)
        shape = np.array(shape)
        if np.any(shape < 2):
            raise Exception('Bravais lattice networks must have at least 2 '
                            'pores in all directions')
        if mode == 'bcc':
            net = bcc(shape=shape, spacing=spacing)
            self['pore.coords'] = net['node.coords']
            self['throat.conns'] = net['edge.conns']
            self['pore.all'] = np.ones(net['node.coords'].shape[0], dtype=bool)
            self['throat.all'] = np.ones(net['edge.conns'].shape[0], dtype=bool)
            # Deal with labels
            self['pore.corner_sites'] = net['node.corner']
            self['pore.body_sites'] = net['node.body']
            Ts = self.find_neighbor_throats(pores=self.pores('body_sites'),
                                            mode='exclusive_or')
            self['throat.corner_to_body'] = False
            self['throat.corner_to_body'][Ts] = True
            Ts = self.find_neighbor_throats(pores=self.pores('corner_sites'),
                                            mode='xnor')
            self['throat.corner_to_corner'] = False
            self['throat.corner_to_corner'][Ts] = True
            Ts = self.find_neighbor_throats(pores=self.pores('body_sites'),
                                            mode='xnor')
            self['throat.body_to_body'] = False
            self['throat.body_to_body'][Ts] = True

        elif mode == 'fcc':
            net = fcc(shape=shape, spacing=spacing)
            self['pore.coords'] = net['node.coords']
            self['throat.conns'] = net['edge.conns']
            self['pore.all'] = np.ones(net['node.coords'].shape[0], dtype=bool)
            self['throat.all'] = np.ones(net['edge.conns'].shape[0], dtype=bool)
            # Deal with labels
            self['pore.face_sites'] = net['node.face']
            self['pore.corner_sites'] = net['node.corner']
            Ts = self.find_neighbor_throats(pores=self.pores('corner_sites'),
                                            mode='xnor')
            self['throat.corner_to_corner'] = False
            self['throat.corner_to_corner'][Ts] = True
            Ts = self.find_neighbor_throats(pores=self.pores('face_sites'))
            self['throat.corner_to_face'] = False
            self['throat.corner_to_face'][Ts] = True

        elif mode == 'hcp':
            raise NotImplementedError('hcp is not implemented yet')

        elif mode == 'sc':
            net = cubic(shape=shape, spacing=1)
            self['pore.coords'] = net['node.coords']
            self['throat.conns'] = net['edge.conns']
            self['pore.all'] = np.ones(net['node.coords'].shape[0], dtype=bool)
            self['throat.all'] = np.ones(net['edge.conns'].shape[0], dtype=bool)
            self['pore.corner_sites'] = True
            self['throat.corner_to_corner'] = True

        else:
            raise Exception('Unrecognized lattice type: ' + mode)

        # Finally scale network to specified spacing
        topotools.label_faces(self)
        Ps = self.pores(['left', 'right', 'top', 'bottom', 'front', 'back'])
        Ps = self.to_mask(pores=Ps)
        self['pore.surface'] = Ps
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
