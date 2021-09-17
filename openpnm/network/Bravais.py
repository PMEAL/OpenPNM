from openpnm.network import GenericNetwork
from openpnm.network.generators import cubic, fcc, bcc
from openpnm import topotools
from openpnm.utils import logging, Workspace
import numpy as np
logger = logging.getLogger(__name__)
ws = Workspace()


class Bravais(GenericNetwork):
    r"""
    Crystal lattice types including fcc, bcc, sc, and hcp

    These arrangements not only allow more dense packing than the standard
    Cubic for higher porosity materials, but also have more interesting
    non-straight connections between the various pore sites.

    More information on Bravais lattice notation can be `found on wikipedia
    <https://en.wikipedia.org/wiki/Bravais_lattice>`_.

    Parameters
    ----------
    shape : array_like
        The number of pores in each direction.  This value is a bit ambiguous
        for the more complex unit cells used here, but generally refers to the
        the number for 'corner' sites
    spacing : array_like (optional)
        The spacing between pores in all three directions.  Like the ``shape``
        this is a bit ambiguous but refers to the spacing between corner sites.
        Essentially it controls the dimensions of the unit cell.  It a scalar
        is given it is applied to all directions.  The default is 1.
    mode : string
        The type of lattice to create.  Options are:

        - 'sc' : Simple cubic (Same as ``Cubic``)
        - 'bcc' : Body-centered cubic lattice
        - 'fcc' : Face-centered cubic lattice
        - 'hcp' : Hexagonal close packed (Not Implemented Yet)

    name : string
        An optional name for the object to help identify it.  If not given,
        one will be generated.
    project : OpenPNM Project object, optional
        Each OpenPNM object must be part of a Project.  If none is supplied
        then one will be created and this Network will be automatically
        assigned to it.  To create a Project use ``openpnm.Project()``.

    See Also
    --------
    Cubic
    CubicDual

    Notes
    -----
    The pores are labelled as beloning to 'corner_sites' and 'body_sites' in
    bcc and additionally 'face_sites' in fcc.  Throats are labelled by which
    type of pores they connect, e.g. 'throat.corner_to_body'.

    """

    def __init__(self, shape, mode='sc', spacing=1, **kwargs):
        super().__init__(**kwargs)
        shape = np.array(shape)
        if np.any(shape < 2):
            raise Exception('Bravais lattice networks must have at least 2 '
                            'pores in all directions')
        if mode == 'bcc':
            # Make a basic cubic for the coner pores
            net1 = cubic(shape=shape)
            # Create a smaller cubic for the body pores, and shift it
            net2 = cubic(shape=shape-1)
            net2['pore.coords'] += 0.5
            # Stitch them together
            net3 = join(net1, net2, L_max=0.99)
            net3['pore.all'] = np.ones(net3['pore.coords'].shape[0], dtype=bool)
            net3['throat.all'] = np.ones(net3['throat.conns'].shape[0], dtype=bool)
            self.update(net3)

            # Deal with labels
            Ps = np.any(np.mod(self['pore.coords'], 1) == 0, axis=1)
            self['pore.body_sites'] = Ps
            self['pore.corner_sites'] = ~Ps

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
            shape = np.array(shape)
            # Create base cubic network of corner sites
            net1 = cubic(shape=shape)
            # Create 3 networks to become face sites
            net2 = cubic(shape=shape - [1, 1, 0])
            net3 = cubic(shape=shape - [1, 0, 1])
            net4 = cubic(shape=shape - [0, 1, 1])
            # Offset
            net2['pore.coords'] += np.array([0.5, 0.5, 0])
            net3['pore.coords'] += np.array([0.5, 0, 0.5])
            net4['pore.coords'] += np.array([0, 0.5, 0.5])
            net2 = join(net2, net3, L_max=0.75)
            net2 = join(net2, net4, L_max=0.75)
            net1 = join(net1, net2, L_max=0.75)
            net1['pore.all'] = np.ones(net1['pore.coords'].shape[0])
            net1['throat.all'] = np.ones(net1['throat.conns'].shape[0])
            self.update(net1)
            # Deal with labels
            self.clear(mode='labels')
            Ps = np.any(np.mod(self['pore.coords'], 1) == 0, axis=1)
            self['pore.face_sites'] = Ps
            self['pore.corner_sites'] = ~Ps
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
            self.update(net)
            self.clear(mode='labels')
            self['pore.corner_sites'] = True
            self['throat.corner_to_corner'] = True

        else:
            raise Exception('Unrecognized lattice type: ' + mode)

        # Finally scale network to specified spacing
        topotools.label_faces(self)
        Ps = self.pores(['left', 'right', 'top', 'bottom', 'front', 'back'])
        Ps = self.tomask(pores=Ps)
        self['pore.surface'] = Ps
        self['pore.internal'] = ~Ps
        self['pore.coords'] *= np.array(spacing)


def join(net1, net2, L_max=0.99):
    # Perform neighbor query
    from scipy.spatial import KDTree
    t1 = KDTree(net1['pore.coords'])
    t2 = KDTree(net2['pore.coords'])
    pairs = t1.query_ball_tree(t2, r=0.99)
    # Combine existing network data
    net3 = {}
    Np1 = net1['pore.coords'].shape[0]
    Np2 = net2['pore.coords'].shape[0]
    net3['pore.coords'] = np.vstack((net1.pop('pore.coords'),
                                     net2.pop('pore.coords')))
    net3['throat.conns'] = np.vstack((net1.pop('throat.conns'),
                                      net2.pop('throat.conns') + Np1))
    # Convert kdtree result into new connections
    nnz = sum([len(row) for row in pairs])
    conns = np.zeros((nnz, 2), dtype=int)
    i = 0
    for j, row in enumerate(pairs):
        for col in row:
            conns[i, :] = j, col + Np1
            i += 1
    # Add new connections to network
    net3['throat.conns'] = np.vstack((net3.pop('throat.conns'), conns))
    # Finally, expand any other data arrays on given networks
    keys = set(net1.keys()).union(net2.keys())
    for item in keys:
        temp1 = net1.pop(item, None)
        temp2 = net2.pop(item, None)
        if temp1 is None:
            temp1 = np.zeros(Np1, dtype=temp2.dtype)
        else:
            temp2 = np.zeros(Np2, dtype=temp1.dtype)
        net3[item] = np.hstack((temp1, temp2)).T
    return net3
