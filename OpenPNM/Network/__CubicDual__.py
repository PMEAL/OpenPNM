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
    def __init__(self, shape=None, spacing=[1, 1, 1], label_1='primary',
                 label_2='secondary', **kwargs):
        super().__init__(**kwargs)
        import OpenPNM as op
        spacing = sp.array(1)
        shape = sp.array([5, 5, 5])
        label_1 = 'primary'
        label_2 = 'secondary'
        spacing = sp.array(spacing)
        shape = sp.array(shape)
        net = op.Network.Cubic(shape=shape, spacing=[1, 1, 1])
        net.add_boundaries()
        net['throat.'+label_1] = True
        net['pore.'+label_1] = True
        dual = op.Network.Cubic(shape=shape+1)
        dual['pore.'+label_2] = True
        dual['throat.'+label_2] = True
        dual['pore.coords'] -= 0.5
        op.Network.tools.stitch(net, dual, P_network=net.Ps, P_donor=dual.Ps,
                                len_max=1)
        net['throat.interconnect'] = net['throat.stitched']
        net['pore.coords'] *= spacing

        # Clean-up
        net['pore.surface'] = False
        net['throat.surface'] = False
        surface_labels = ['top', 'bottom', 'front', 'back', 'left', 'right']
        for face in surface_labels:
            Ps = net.pores(labels=[face, label_1], mode='intersection')
            net['pore.'+face][Ps] = False
            Ps = net.pores(labels=[face+'_boundary'])
            net['pore.surface'][Ps] = True
            net['pore.'+face][Ps] = True

        [net.pop(item) for item in list(net.keys()) if 'boundary' in item]
        # Label all remaining 'face' pores as 'surface'
        Ps = net.pores(labels=surface_labels)
        net['pore.surface'][Ps] = True
        Ts = net.find_neighbor_throats(pores=net.pores('surface'),
                                       mode='intersection')
        net['throat.surface'][Ts] = True
        # Label non-surface pores and throats as internal
        net['pore.internal'] = ~net['pore.surface']
        Ts = net.find_neighbor_throats(pores=net['pore.internal'])
        net['throat.internal'] = False
        net['throat.internal'][Ts] = True
        # Remove unused labels
        del net['throat.stitched']
        # Transfer all dictionary items from 'net' to 'self'
        [self.update({item: net[item]}) for item in net]
        del self.workspace[net.name]
