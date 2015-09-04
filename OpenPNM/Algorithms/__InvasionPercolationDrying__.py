# -*- coding: utf-8 -*-
"""
===============================================================================
InvasionPercolationDrying: IP for drying of a wetting phase
===============================================================================

"""
import heapq as hq
import scipy as sp
from collections import namedtuple
from OpenPNM.Algorithms import GenericAlgorithm
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)
Tinfo = namedtuple('queued_throat', ('order', 'number'))


class InvasionPercolationDrying(GenericAlgorithm):
    r"""
    A classic/basic invasion percolation algorithm optimized for speed.

    Parameters
    ----------
    network : OpenPNM Network object
        The Network upon which the invasion should occur.

    Notes
    ----
    n/a

    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.queues = []

    def setup(self, invading_phase, defending_phase,
              throat_prop='throat.capillary_pressure'):
        r"""
        Set up the required parameters for the algorithm

        Parameters
        ----------
        invading_phase : OpenPNM Phase object
            The phase to be injected into the Network.  The Phase must have the
            capillary entry pressure values for the system.

        throat_prop : string
            The name of the throat property containing the capillary entry
            pressure.  The default is 'throat.capillary_pressure'.

        """
        self._phase = invading_phase
        # Setup arrays and info
        self['throat.entry_pressure'] = defending_phase[throat_prop]
        # Indices into t_entry giving a sorted list
        self['throat.sorted'] = sp.argsort(self['throat.entry_pressure'], axis=0)
        self['throat.order'] = sp.zeros_like(self['throat.sorted'])
        self['throat.order'][self['throat.sorted']] = sp.arange(0, self._net.Nt)
        self['throat.invaded'] = -sp.ones((self._net.Nt,))
        self['pore.invaded'] = -sp.ones((self._net.Np,))
        self['pore.inlets'] = False
        self['pore.outlets'] = False
        self['pore.volume'] = sp.copy(self._net['pore.volume'])
        self['throat.volume'] = sp.copy(self._net['throat.volume'])
        self._tcount = 0

    def set_inlets(self, pores):
        r"""

        Parameters
        ----------
        pores : array_like
            The list of inlet pores from which the invading Phase can enter the
            Network
        """
        self['pore.inlets'][pores] = True

        self['pore.invaded'][pores] = 0

        self.qregen()

    def set_outlets(self, pores):
        r"""

        Parameters
        ----------
        pores : array_like
            The list of outlet pores from which the defending Phase can leave
            the Network
        """
        self['pore.outlets'][pores] = True

    def run(self, n_steps=None):
        if n_steps is None:
            n_steps = sp.inf
        count = 0
        while count < n_steps:
            count += 1
            for q in self.queues:
                T = self.qpop(q)
                self['throat.invaded'][T] = 1
                Ps = self._net['throat.conns'][T]
                # Remove invaded pore(s)
                P = Ps[self['pore.invaded'][Ps] < 0]
                if len(P) > 0:
                    # Note newly invaded pore
                    self['pore.invaded'][P] = 1
                    # Add newly connected throats to queue
                    Ts = self._net.find_neighbor_throats(pores=P)
                    # Remove already invaded throats from Ts
                    Ts = Ts[self['throat.invaded'][Ts] < 0]
                    self.qpush(queue=q, throats=Ts)
            # Remove empty queues from list of queues
            self.queues = [x for x in self.queues if x]

    def run2(self, queue, volume):
        while len(queue) > 0:
            # Check if throat is fully drained, if not do partial step
            T = queue[0][1]  # Get top throat, but don't pop it yet
            if volume < self['throat.volume'][T]:
                self['throat.volume'][T] -= volume
                v = str(self['throat.volume'][T])
                print('Throat '+str(T)+' partially drained to '+v)
                break
            volume -= self['throat.volume'][T]  # Remove throat vol from total
            self['throat.volume'][T] = 0  # Set throat vol to 0
            print('Throat '+str(T)+' completely drained, volume remaining is '+str(volume))
            # Find uninvaded neighbor pores
            Ps = self._net['throat.conns'][T]
            P = Ps[self['pore.invaded'][Ps] < 0]
            if len(P) > 0:
                # Check if pore volume is fully drained, do a partial step
                if volume < self['pore.volume'][P]:
                    self['pore.volume'][P] -= volume
                    v = str(self['pore.volume'][P[0]])
                    print('Pore '+str(P[0])+' partially drained to '+v)
                    break
                volume -= self['pore.volume'][P]
                self['pore.volume'][P] = 0
                print('Pore '+str(P[0])+' completely drained, volume remaining is '+str(volume))
            # Since both throat and pore can be drained, pop throat from queue
            self.qpop(queue)
            # Set both throat and pore to invaded
            self['throat.invaded'][T] = 1
            self['pore.invaded'][P] = 1
            # Find newly connected throats
            Ts = self._net.find_neighbor_throats(pores=P)
            Ts = Ts[self['throat.invaded'][Ts] < 0]
            # Add new throats to the queue
            self.qpush(queue=queue, throats=Ts)

    def qregen(self):
        self.queues = []
        Pdef = self['pore.invaded'] == -1
        clusters = self._net.find_clusters2(mask=Pdef)
        self['pore.queue_number'] = -1
        label = 0
        for c in sp.unique(clusters)[1:]:
            Ps = sp.where(clusters == c)[0]
            Ts = self._net.find_neighbor_throats(pores=Ps,
                                                 mode='not_intersection')
            q = self.qmake(Ts)
            self.queues.append(q)
            self['pore.queue_number'][Ps] = label
            label += 1

    def qmake(self, throats):
        Ts = self._parse_locations(throats)
        order = self['throat.order']
        queue = []
        [hq.heappush(queue, Tinfo(order[T], T)) for T in Ts]
        return queue

    def qpop(self, queue):
        Tup = hq.heappop(queue)
        T = Tup[1]
        return T

    def qpush(self, queue, throats):
        Ts = self._parse_locations(throats)
        order = self['throat.order']
        [hq.heappush(queue, Tinfo(order[T], T)) for T in Ts]


























