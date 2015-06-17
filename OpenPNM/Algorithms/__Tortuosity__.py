# -*- coding: utf-8 -*-
"""
========================================================================
Tortuosity: Network Tortuosity Algorithm
========================================================================
This algorithm uses Djkstra's algorithm to get the shortest path between
two points folliwng the network, and the the direct distance between the same
points.  The ratio of these is returned as the 'tortuosity'

TODO: It currently uses the 'throat.length' to weight the network connections
but this should probably use diffusive conductance.
"""

import scipy as sp
import scipy.sparse.csgraph as spgr
from OpenPNM.Algorithms import GenericAlgorithm
import OpenPNM.Network
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)


class Tortuosity(GenericAlgorithm):
    r"""
    Determines the tortuosity of the network using a shortest path search algorithm.

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        logger.debug('Create Tortuosity Object')

    def estimate_time(self):
        pn_temp = OpenPNM.Network.TestNet()
        graph = pn_temp.create_adjacency_matrix(sprsfmt='csr')
        self._net.tic()
        path = spgr.shortest_path(csgraph=graph, method='D', directed=False)
        t = self._net.toc(quiet=True)
        N = 125
        k = 6
        O = t / (N * (N * k + N * sp.log(N)))
        N = self._net.num_pores()
        k = sp.median(self._net.num_neighbors(pores=self._net.pores()))
        t_est = O * (N * (N * k + N * sp.log(N)))
        print('Based on the network size and PC performance, this algorithm \
              will require: ', t_est, ' seconds')
        return

    def run(self, phase=None):
        logger.warning('This algorithm can take some time...')
        graph = self._net.create_adjacency_matrix(data=self._net['throat.length'],
                                                  sprsfmt='csr')

        if phase is not None:
            self._phase = phase
            if 'throat.occupancy' in self._phase.props():
                temp = self._net['throat.length'] * \
                    (self._phase['throat.occupancy'] == 1)
                graph = self._net.create_adjacency_matrix(data=temp,
                                                          sprsfmt='csr',
                                                          prop='temp')
        path = spgr.shortest_path(csgraph=graph, method='D', directed=False)

        Px = sp.array(self._net['pore.coords'][:, 0], ndmin=2)
        Py = sp.array(self._net['pore.coords'][:, 1], ndmin=2)
        Pz = sp.array(self._net['pore.coords'][:, 2], ndmin=2)

        Cx = sp.square(Px.T - Px)
        Cy = sp.square(Py.T - Py)
        Cz = sp.square(Pz.T - Pz)
        Ds = sp.sqrt(Cx + Cy + Cz)

        temp = path / Ds

        temp[sp.isnan(temp)] = 0
        temp[sp.isinf(temp)] = 0

        return temp
