# -*- coding: utf-8 -*-
"""
===============================================================================
InvasionPercolationBasic: Simple IP
===============================================================================

"""

import scipy as sp
import bisect
from collections import deque
from OpenPNM.Algorithms import GenericAlgorithm
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)


class InvasionPercolation(GenericAlgorithm):
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

    def run(self, phase, inlets, throat_prop='throat.capillary_pressure'):
        r"""
        Perform the algorithm

        Parameters
        ----------
        phase : OpenPNM Phase object
            The phase to be injected into the Network.  The Phase must have the
            capillary entry pressure values for the system.

        inlets : array_like
            The list of inlet pores from which the Phase can enter the Network

        throat_prop : string
            The name of the throat property containing the capillary entry
            pressure.  The default is 'throat.capillary_pressure'.

        """
        import heapq as hq
        queue = []
        hq.heapify(queue)
        self._phase = phase
        net = self._net
        # Setup arrays and info
        t_entry = phase[throat_prop]
        # Indices into t_entry giving a sorted list
        t_sorted = sp.argsort(t_entry, axis=0)
        t_order = sp.zeros_like(t_sorted)
        t_order[t_sorted] = sp.arange(0, net.Nt)
        t_inv = -sp.ones_like(net.Ts)
        p_inv = -sp.ones_like(net.Ps)
        p_inv[inlets] = 0
        # Perform initial analysis on input pores
        Ts = net.find_neighbor_throats(pores=inlets)
        [hq.heappush(queue, T) for T in t_order[Ts]]
        tcount = 1
        while len(queue) > 0:
            # Find throat at the top of the queue
            t = hq.heappop(queue)
            # Extract actual throat number
            t_next = t_sorted[t]
            t_inv[t_next] = tcount
            # If throat is duplicated
            while len(queue) > 0 and queue[0] == t:
                # Note: Preventing duplicate entries below might save some time here
                t = hq.heappop(queue)
            # Find pores connected to newly invaded throat
            Ps = net['throat.conns'][t_next]
            # Remove already invaded pores from Ps
            Ps = Ps[p_inv[Ps] < 0]
            if len(Ps) > 0:
                p_inv[Ps] = tcount
                Ts = net.find_neighbor_throats(pores=Ps)  # Find connected throats
                Ts = Ts[t_inv[Ts] < 0]  # Remove already invaded throats from Ts
                [hq.heappush(queue, T) for T in t_order[Ts]]
            tcount += 1
        self['throat.invasion_sequence'] = t_inv
        self['pore.invasion_sequence'] = p_inv

    def return_results(self, pores=[], throats=[]):
        r"""
        Places the results of the IP simulation into the Phase object.

        Parameters
        ----------
        pores and throats : array_like
            The list of pores and throats whose values should be returned to
            the Phase object.  Default is all of them.

        Returns
        -------
        invasion_sequence : array_like
            The sequence in which each pore and throat is invaded  This depends
            on the inlet locations.  All inlets are invaded at step 0.  It is
            possible to recontruct an animation of the invasion process, in
            Paraview for instance, using this sequence information.

        """
        pores = sp.array(pores, ndmin=1)
        throats = sp.array(throats, ndmin=1)
        if len(pores) == 0:
            pores = self.Ps
        if len(throats) == 0:
            throats = self.Ts
        self._phase['throat.invasion_sequence'] = sp.nan
        self._phase['pore.invasion_sequence'] = sp.nan
        self._phase['throat.invasion_sequence'][throats] = \
            self['throat.invasion_sequence'][throats]
        self._phase['pore.invasion_sequence'][pores] = \
            self['pore.invasion_sequence'][pores]

    def apply_flow(self, flowrate):
        r"""
        Convert the invaded sequence into an invaded time for a given flow rate
        considering the volume of invaded pores and throats.

        Parameters
        ----------
        flowrate : float
            The flow rate of the injected fluid

        Returns
        -------
        Creates a throat array called 'invasion_time' in the Algorithm
        dictionary

        """
        P12 = self._net['throat.conns']
        a = self['throat.invasion_sequence']
        b = sp.argsort(self['throat.invasion_sequence'])
        P12_inv = self['pore.invasion_sequence'][P12]
        # Find if the connected pores were invaded with or before each throat
        P1_inv = P12_inv[:, 0] == a
        P2_inv = P12_inv[:, 1] == a
        c = sp.column_stack((P1_inv, P2_inv))
        d = sp.sum(c, axis=1, dtype=bool)  # List of Pores invaded with each throat
        # Find volume of these pores
        P12_vol = sp.zeros((self.Nt,))
        P12_vol[d] = self._net['pore.volume'][P12[c]]
        # Add invaded throat volume to pore volume (if invaded)
        T_vol = P12_vol + self._net['throat.volume']
        # Cumulative sum on the sorted throats gives cumulated inject volume
        e = sp.cumsum(T_vol[b] / flowrate)
        t = sp.zeros((self.Nt,))
        t[b] = e  # Convert back to original order
        self._phase['throat.invasion_time'] = t
