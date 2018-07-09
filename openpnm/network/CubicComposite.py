# -*- coding: utf-8 -*-
"""
===============================================================================
CubicComposite: Generate lattice-like networks with different regions
===============================================================================

"""
import numpy as np
from openpnm.network import GenericNetwork
from openpnm.network import Cubic
from openpnm.topotools import extend, stitch, label_faces

shapes = [[10, 10, 10], [20, 20, 5], [5, 10, 8]]
spacings = [1, 0.5, [2, 1, 0.8]]
origins = [[0, 0, 0], [0, 0, 10], [0, 0, 12.5]]
labels = None


class CubicComposite(GenericNetwork):
    r"""
    """
    def __init__(self, shapes, spacings, origins, labels=None, name=None,
                 project=None):

        # If no labels given, then use default
        if labels is None:
            labels = [None for i in shapes]

        # Generate cubic networks of specified properties and store in a list
        nets = []
        for n in range(len(shapes)):
            net = Cubic(shape=shapes[n], spacing=spacings[n], name=labels[n])
#            label_faces(net)
            net['pore.coords'] += origins[n]
            nets.append(net)

        # Begin the process of merging the networks, but no stitching yet
        main_net = nets.pop(0)
        for net in nets:
            extend(network=main_net, pore_coords=net['pore.coords'],
                   throat_conns=net['throat.conns'] + main_net.Np)

        for net in nets:
            P1 = main_net.pores('surface')
            P1 = net.pores('surface')
            stitch(network=main_net)
