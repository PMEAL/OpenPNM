# -*- coding: utf-8 -*-
"""
===============================================================================
CubicDUal: Generate lattice-like networks with a
===============================================================================

"""
import scipy as sp
import OpenPNM.Utilities.misc as misc
from OpenPNM.Network import GenericNetwork
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)


class CubicDual(GenericNetwork):
    r"""
    Generates a cubic network of the specified size and shape.

    Parameters
    ----------
    name : string
        A unique name for the network

    shape : list of ints
        The size and shape of the principal cubic network in terms of the
        number of pores in each direction.  Secondary nodes will be added at
        the center of each unit cell defined by a cube 8 pores in the primary
        network.

    spacing : list of floats
        The distance between pores in each of the principal directions

    label_1 and label_2 : strings
        The labels to apply to the main cubic lattice and the interpenetrating
        cubic lattice (i.e. the dual network).  The defaults are 'corners' and
        'centers', which refers to the position on the unit cell.

    Examples
    --------
    >>>
    """
    def __init__(self, shape=None, spacing=[1, 1, 1], label_1='corners',
                 label_2='centers', **kwargs):
        super().__init__(**kwargs)
        import OpenPNM as op
        spacing = sp.array(spacing)
        shape = sp.array(shape)
        net = op.Network.Cubic(shape=shape, spacing=[1, 1, 1])
        net['throat.'+label_1] = True
        net['pore.'+label_1] = True
        dual = op.Network.Cubic(shape=shape-1)
        dual['pore.'+label_2] = True
        dual['throat.'+label_2] = True
        dual['pore.coords'] += 0.5
        op.Network.tools.stitch(net, dual, P_network=net.Ps,
                                P_donor=dual.Ps, len_max=1)
        net['throat.interconnect'] = net['throat.stitched']
        del net['throat.stitched']
        net['pore.surface'] = False
        net['throat.surface'] = False
        # Create center pores on each surface
        offset = {'back': [0.5, 0, 0], 'front': [-0.5, 0, 0],
                  'right': [0, 0.5, 0], 'left': [0, -0.5, 0],
                  'top': [0, 0, 0.5], 'bottom': [0, 0, -0.5]}
        surface_labels = ['top', 'bottom', 'left', 'right', 'front', 'back']
        for item in surface_labels:
            # Clone 'center' pores and shift to surface
            Ps = net.pores(labels=[item, label_2], mode='intersection')
            net.clone_pores(pores=Ps, apply_label=[label_2, 'surface', item])
            Ps = net.pores(labels=['surface', item], mode='intersection')
            net['pore.coords'][Ps] += offset[item]
            Ts = net.find_neighbor_throats(pores=Ps)
            # Label pores and throats
            net['pore.surface'][Ps] = True
            net['throat.'+label_2][Ts] = True
        Ps = net.pores(labels=surface_labels)
        net['pore.surface'][Ps] = True

        net['pore.coords'] *= spacing
        [self.update({item: net[item]}) for item in net]
        del self.workspace[net.name]

    def add_boundary_pores(self, pores, offset, apply_label='boundary'):
        r"""
        This method uses ``clone_pores`` to clone the input pores, then shifts
        them the specified amount and direction, then applies the given label.

        Parameters
        ----------
        pores : array_like
            List of pores to offset.  If no pores are specified, then it
            assumes that all surface pores are to be cloned.

        offset : 3 x 1 array
            The distance in vector form which the cloned boundary pores should
            be offset.

        apply_label : string
            This label is applied to the boundary pores.  Default is
            'boundary'.

        Examples
        --------
        >>> import OpenPNM as op
        >>> pn = op.Network.Cubic(shape=[5, 5, 5])
        >>> print(pn.Np)  # Confirm initial Network size
        125
        >>> Ps = pn.pores('top')  # Select pores on top face
        >>> pn.add_boundary_pores(pores=Ps, offset=[0, 0, 1])
        >>> print(pn.Np)  # Confirm addition of 25 new pores
        150
        >>> 'pore.boundary' in pn.labels()  # Default label is created
        True
        """
        # Parse the input pores
        Ps = sp.array(pores, ndmin=1)
        if Ps.dtype is bool:
            Ps = self.toindices(Ps)
        if sp.size(pores) == 0:  # Handle an empty array if given
            return sp.array([], dtype=sp.int64)
        # Clone the specifed pores
        self.clone_pores(pores=Ps)
        newPs = self.pores('pore.clone')
        del self['pore.clone']
        # Offset the cloned pores
        self['pore.coords'][newPs] += offset
        # Apply labels to boundary pores (trim leading 'pores' if present)
        label = apply_label.split('.')[-1]
        label = 'pore.' + label
        logger.debug('The label \''+label+'\' has been applied')
        self[label] = False
        self[label][newPs] = True
