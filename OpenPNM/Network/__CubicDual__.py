# -*- coding: utf-8 -*-
"""
===============================================================================
CubicDual: Generate a cubic lattice with an interpentrating dual network
===============================================================================

"""
import scipy as sp
from OpenPNM.Network.tools import stitch
from OpenPNM.Network import GenericNetwork, Cubic
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
        spacing = sp.array(spacing)
        shape = sp.array(shape)
        # Deal with non-3D shape arguments
        shape = sp.pad(shape, [0, 3-shape.size], mode='constant', constant_values=1)
        net = Cubic(shape=shape, spacing=[1, 1, 1])
        net['throat.'+label_1] = True
        net['pore.'+label_1] = True
        single_dim = shape == 1
        shape[single_dim] = 2
        dual = Cubic(shape=shape-1, spacing=[1, 1, 1])
        faces = [['front', 'back'], ['left', 'right'], ['top', 'bottom']]
        faces = [faces[i] for i in sp.where(~single_dim)[0]]
        faces = sp.array(faces).flatten().tolist()
        dual.add_boundaries(faces)
        # Add secondary network name as a label
        dual['pore.'+label_2] = True
        dual['throat.'+label_2] = True
        # Shift coordinates prior to stitching
        dual['pore.coords'] += 0.5*(~single_dim)
        stitch(net, dual, P_network=net.Ps, P_donor=dual.Ps, len_max=1)
        net['throat.interconnect'] = net['throat.stitched']
        del net['throat.stitched']
        net['pore.coords'] *= spacing
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
        del self.workspace[net.name]
