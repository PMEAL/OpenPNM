import scipy as sp
from openpnm.network import GenericNetwork, Cubic
from openpnm import topotools
from openpnm.utils import logging, Workspace
logger = logging.getLogger(__name__)
ws = Workspace()


class CubicDual(GenericNetwork):
    r"""
    Body centered cubic lattice plus face centered nodes on the surfaces

    This network is essentially a *'bcc'* lattice, *except* that the seconary
    network (body-centered pores) has pores on each face of the domain, which
    breaks the body-centric arranagement.  This allows boundary conditions to
    be applied to the seconary network for transport simuations.

    Parameters
    ----------
    shape : list of ints
        The size and shape of the primary cubic network in terms of the
        number of pores in each direction.  Secondary nodes will be added at
        centers of each unit cell.

    spacing : list of floats
        The distance between pores of the primary network in each of the
        principal directions

    label_1 : string
        The label to apply to the primary cubic lattices, which defaults to
        'primary'

    label_2 : string
        The label to apply to the secondary cubic lattices, which defaults to
        'seconary'

    project : OpenPNM Project object (optional)
        If not provided one will be generated and the network will be assigned
        to it.  It can be retrieved from ``net.project``.

    name : string
        An optional name for the object to help identify it.  If not given,
        one will be generated.

    See Also
    --------
    Bravais

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.CubicDual(shape=[3, 3, 3])
    >>> pn.num_pores('pore.primary')  # Normal cubic network is present
    27
    >>> pn.Np  # But more pores are present from seconary network
    59

    And it can be plotted for quick visualization using:

    >>> fig = op.topotools.plot_connections(network=pn,
    ...                                     throats=pn.throats('primary'),
    ...                                     color='b')
    >>> fig = op.topotools.plot_connections(network=pn,
    ...                                     throats=pn.throats('secondary'),
    ...                                     color='r')
    >>> fig = op.topotools.plot_coordinates(network=pn, c='r', s=75, fig=fig)

    .. image:: /../docs/static/images/cubic_dual_network.png
        :align: center

    For larger networks and more control over presentation use `Paraview
    <http://www.paraview.org>`_.

    """
    def __init__(self, shape, spacing=1, label_1='primary',
                 label_2='secondary', **kwargs):
        super().__init__(**kwargs)
        spacing = sp.array(spacing)
        shape = sp.array(shape)
        # Deal with non-3D shape arguments
        shape = sp.pad(shape, [0, 3-shape.size], mode='constant',
                       constant_values=1)
        net = Cubic(shape=shape, spacing=1)
        net['throat.'+label_1] = True
        net['pore.'+label_1] = True
        single_dim = shape == 1
        shape[single_dim] = 2
        dual = Cubic(shape=shape-1, spacing=1)
        faces = [['front', 'back'], ['left', 'right'], ['top', 'bottom']]
        faces = [faces[i] for i in sp.where(~single_dim)[0]]
        faces = sp.array(faces).flatten().tolist()
        dual.add_boundary_pores(faces)
        # Add secondary network name as a label
        dual['pore.'+label_2] = True
        dual['throat.'+label_2] = True
        # Shift coordinates prior to stitching
        dual['pore.coords'] += 0.5*(~single_dim)
        topotools.stitch(net, dual, P_network=net.Ps, P_donor=dual.Ps,
                         len_max=1)
        net['throat.interconnect'] = net['throat.stitched']
        del net['throat.stitched']
        # Clean-up labels
        net['pore.surface'] = False
        net['throat.surface'] = False
        for face in faces:
            # Remove face label from secondary network since it's internal now
            Ps = net.pores(labels=[face, label_2], mode='intersection')
            net['pore.'+face][Ps] = False
            Ps = net.pores(labels=[face+'_boundary'])
            net['pore.'+face][Ps] = True
            Ps = net.pores(face)
            net['pore.surface'][Ps] = True
            Ts = net.find_neighbor_throats(pores=Ps, mode='intersection')
            net['throat.surface'][Ts] = True
            net['throat.'+face] = net.tomask(throats=Ts)
        [net.pop(item) for item in net.labels() if 'boundary' in item]
        # Label non-surface pores and throats as internal
        net['pore.internal'] = ~net['pore.surface']
        Ts = net.find_neighbor_throats(pores=net['pore.internal'])
        net['throat.internal'] = False
        net['throat.internal'][Ts] = True
        # Transfer all dictionary items from 'net' to 'self'
        [self.update({item: net[item]}) for item in net]
        ws.close_project(net.project)
        # Finally, scale network to requested spacing
        net['pore.coords'] *= spacing

    def add_boundary_pores(self, labels=['top', 'bottom', 'front', 'back',
                                         'left', 'right'], spacing=None):
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
        spacing = sp.array(spacing)
        if spacing.size == 1:
            spacing = sp.ones(3)*spacing
        for item in labels:
            Ps = self.pores(item)
            coords = sp.absolute(self['pore.coords'][Ps])
            axis = sp.count_nonzero(sp.diff(coords, axis=0), axis=0) == 0
            offset = sp.array(axis, dtype=int)*spacing/2
            if sp.amin(coords) == sp.amin(coords[:, sp.where(axis)[0]]):
                offset = -1*offset
            topotools.add_boundary_pores(network=self, pores=Ps, offset=offset,
                                         apply_label=item + '_boundary')
