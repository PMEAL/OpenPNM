from openpnm.network import GenericNetwork, Cubic
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
        - 'hcp' : Hexagonal close packed (Note Implemented Yet)

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
    bcc or 'face_sites' in fcc.  Throats are labelled by the which type of
    pores they connect, e.g. 'throat.corner_to_body'.

    Limitations:

    * Bravais lattice can also have a skew to them, but this is not
      implemented yet.
    * Support for 2D networks has not been added yet.
    * Hexagonal Close Packed (hcp) has not been implemented yet, but is on
      the todo list.

    Examples
    --------
    >>> import openpnm as op
    >>> sc = op.network.Bravais(shape=[3, 3, 3], mode='sc')
    >>> bcc = op.network.Bravais(shape=[3, 3, 3], mode='bcc')
    >>> fcc = op.network.Bravais(shape=[3, 3, 3], mode='fcc')
    >>> sc.Np, bcc.Np, fcc.Np
    (27, 35, 63)

    Since these three networks all have the same domain size, it is clear that
    both 'bcc' and 'fcc' have more pores per unit volume.  This is particularly
    helpful for modeling higher porosity materials.

    They all have the same number corner sites, which corresponds to the
    [3, 3, 3] shape that was specified:

    >>> sc.num_pores('corner*'), bcc.num_pores('cor*'), fcc.num_pores('cor*')
    (27, 27, 27)

    Visualization of these three networks can be done quickly using the
    functions in topotools.  Firstly, merge them all into a single network
    for convenience:

    >>> bcc['pore.coords'][:, 0] += 3
    >>> fcc['pore.coords'][:, 0] += 6
    >>> op.topotools.merge_networks(sc, [bcc, fcc])
    >>> fig = op.topotools.plot_connections(sc)

    .. image:: /../docs/_static/images/bravais_networks.png
        :align: center

    For larger networks and more control over presentation use `Paraview
    <http://www.paraview.org>`_.

    """

    def __init__(self, shape, mode='sc', spacing=1, **kwargs):
        super().__init__(**kwargs)
        shape = np.array(shape)
        if np.any(shape < 2):
            raise Exception('Bravais lattice networks must have at least 2 '
                            'pores in all directions')
        if mode == 'bcc':
            # Make a basic cubic for the coner pores
            net1 = Cubic(shape=shape)
            net1['pore.net1'] = True
            # Create a smaller cubic for the body pores, and shift it
            net2 = Cubic(shape=shape-1)
            net2['pore.net2'] = True
            net2['pore.coords'] += 0.5
            # Stitch them together
            topotools.stitch(net1, net2, net1.Ps, net2.Ps, len_max=0.99)
            self.update(net1)
            ws.close_project(net1.project)

            # Deal with labels
            Ps1 = self['pore.net2']
            self.clear(mode='labels')
            self['pore.corner_sites'] = ~Ps1
            self['pore.body_sites'] = Ps1
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
            net1 = Cubic(shape=shape)
            # Create 3 networks to become face sites
            net2 = Cubic(shape=shape - [1, 1, 0])
            net3 = Cubic(shape=shape - [1, 0, 1])
            net4 = Cubic(shape=shape - [0, 1, 1])
            net2['pore.coords'] += np.array([0.5, 0.5, 0])
            net3['pore.coords'] += np.array([0.5, 0, 0.5])
            net4['pore.coords'] += np.array([0, 0.5, 0.5])
            # Remove throats from net2 (trim doesn't work when removing ALL)
            for n in [net2, net3, net4]:
                n.clear(element='throat', mode='all')
                n.update({'throat.all': np.array([], dtype=bool)})
                n.update({'throat.conns': np.ndarray([0, 2], dtype=bool)})
            # Join networks 2, 3 and 4 into one with all face sites
            topotools.stitch(net2, net3, net2.Ps, net3.Ps, method='radius',
                             len_max=0.75)
            topotools.stitch(net2, net4, net2.Ps, net4.Ps, method='radius',
                             len_max=0.75)
            # Join face sites network with the corner sites network
            topotools.stitch(net1, net2, net1.Ps, net2.Ps, method='radius',
                             len_max=0.75)
            self.update(net1)
            ws.close_project(net1.project)
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
            net = Cubic(shape=shape, spacing=1)
            self.update(net)
            ws.close_project(net.project)
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

    def add_boundary_pores(self, labels, spacing):
        r"""
        Add boundary pores to the specified faces of the network

        Pores are offset from the faces by 1/2 of the given ``spacing``, such
        that they lie directly on the boundaries.

        Parameters
        ----------
        labels : string or list of strings
            The labels indicating the pores defining each face where boundary
            pores are to be added (e.g. 'left' or ['left', 'right'])

        spacing : scalar or array_like
            The spacing of the network (e.g. [1, 1, 1]).  This must be given
            since it can be quite difficult to infer from the network,
            for instance if boundary pores have already added to other faces.

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
