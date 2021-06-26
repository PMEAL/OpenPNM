import numpy as np
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
    >>> import matplotlib.pyplot as plt
    >>> pn = op.network.CubicDual(shape=[3, 3, 3])
    >>> pn.num_pores('pore.primary')  # Normal cubic network is present
    27
    >>> pn.Np  # But more pores are present from seconary network
    59

    And it can be plotted for quick visualization using:

    >>> fig, ax = plt.subplots()
    >>> _ = op.topotools.plot_connections(network=pn,
    ...                                   throats=pn.throats('primary'),
    ...                                   color='b', ax=ax)
    >>> _ = op.topotools.plot_connections(network=pn,
    ...                                   throats=pn.throats('secondary'),
    ...                                   color='r', ax=ax)
    >>> _ = op.topotools.plot_coordinates(network=pn, c='r', s=75, ax=ax)

    .. image:: /../docs/_static/images/cubic_dual_network.png
        :align: center

    For larger networks and more control over presentation use `Paraview
    <http://www.paraview.org>`_.

    """
    def __init__(self, shape, spacing=1, label_1='primary',
                 label_2='secondary', **kwargs):
        super().__init__(**kwargs)
        spacing = np.array(spacing)
        shape = np.array(shape)
        # Deal with non-3D shape arguments
        shape = np.pad(shape, [0, 3-shape.size], mode='constant',
                       constant_values=1)
        net = Cubic(shape=shape, spacing=1)
        net['throat.'+label_1] = True
        net['pore.'+label_1] = True
        single_dim = shape == 1
        shape[single_dim] = 2
        dual = Cubic(shape=shape-1, spacing=1)
        faces = [['left', 'right'], ['front', 'back'], ['top', 'bottom']]
        faces = [faces[i] for i in np.where(~single_dim)[0]]
        faces = np.array(faces).flatten().tolist()
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
            Ps = net.pores(labels=[face, label_2], mode='xnor')
            net['pore.'+face][Ps] = False
            Ps = net.pores(labels=[face+'_boundary'])
            net['pore.'+face][Ps] = True
            Ps = net.pores(face)
            net['pore.surface'][Ps] = True
            Ts = net.find_neighbor_throats(pores=Ps, mode='xnor')
            net['throat.surface'][Ts] = True
            net['throat.'+face] = net.tomask(throats=Ts)
        for item in net.labels():
            if 'boundary' in item:
                net.pop(item)
        # Label non-surface pores and throats as internal
        net['pore.internal'] = True
        net['throat.internal'] = True
        # Transfer all dictionary items from 'net' to 'self'
        for item in net:
            self.update({item: net[item]})
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
        spacing = np.array(spacing)
        if spacing.size == 1:
            spacing = np.ones(3)*spacing
        for item in labels:
            Ps = self.pores(item)
            coords = np.absolute(self['pore.coords'][Ps])
            axis = np.count_nonzero(np.diff(coords, axis=0), axis=0) == 0
            offset = np.array(axis, dtype=int)*spacing/2
            if np.amin(coords) == np.amin(coords[:, np.where(axis)[0]]):
                offset = -1*offset
            topotools.add_boundary_pores(network=self, pores=Ps, offset=offset,
                                         apply_label=item + '_boundary')
