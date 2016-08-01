# -*- coding: utf-8 -*-
"""
===============================================================================
CubicDual: Generate a cubic lattice with an interpentrating dual network
===============================================================================

"""
import scipy as sp
import OpenPNM.Utilities.misc as misc
from OpenPNM.Network import GenericNetwork
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)


class CubicDual(GenericNetwork):
    r"""
    Generates a cubic network of the specified size and shape, then inserts
    as second cubic network to the corners of its lattice cells.  These two
    networks are further connected by throats enabling material to be
    exchanged between them.

    Parameters
    ----------
    name : string
        A unique name for the network

    shape : list of ints
        The size and shape of the primary cubic network in terms of the
        number of pores in each direction.  Secondary nodes will be added at
        corners of each unit cell so the dual network will generaly have a
        size of ''shape'' + 1.

    spacing : list of floats
        The distance between pores of the primary network in each of the
        principal directions

    label_1 and label_2 : strings
        The labels to apply to the primary and secondary cubic lattices, which
        defaults to 'primary' and 'secondary' respectively.

    Examples
    --------
    >>>
    """
    def __init__(self, shape=None, spacing=[1, 1, 1], label_1='corners',
                 label_2='centers', **kwargs):
        super().__init__(**kwargs)
        import OpenPNM as op
        spacing = sp.array(1)
        shape = sp.array([5, 5, 5])
        label_1 = 'primary'
        label_2 = 'secondary'
        net = op.Network.Cubic(shape=shape, spacing=[1, 1, 1])
        net.add_boundaries()
        net['throat.'+label_1] = True
        net['pore.'+label_1] = True
        dual = op.Network.Cubic(shape=shape+1)
        dual['pore.'+label_2] = True
        dual['throat.'+label_2] = True
        dual['pore.coords'] -= 0.5
        op.Network.tools.stitch(net, dual, P_network=net.Ps,
                                P_donor=dual.Ps, len_max=1)
        net['throat.interconnect'] = net['throat.stitched']
        del net['throat.stitched']
        net['pore.surface'] = False
        net['throat.surface'] = False


        # Clean-ups
        net['pore.coords'] *= spacing
        Ps = net.pores(labels=surface_labels)
        net['pore.surface'][Ps] = True
        net['pore.internal'] = ~net['pore.surface']
        Ts = net.find_neighbor_throats(pores=net['pore.internal'])
        net['throat.internal'] = False
        net['throat.internal'][Ts] = True
        del net['pore._clone']
        del net['throat._clone']
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
