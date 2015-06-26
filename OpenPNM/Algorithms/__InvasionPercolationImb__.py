# -*- coding: utf-8 -*-
"""
===============================================================================
InvasionPercolationImb: Simple IP for imbibition
===============================================================================

"""

import scipy as sp
import bisect
from collections import deque
from OpenPNM.Algorithms import GenericAlgorithm
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)


class InvasionPercolationImb(GenericAlgorithm):
    r"""
    A classic/basic invasion percolation algorithm optimized for speed.
    Imbibition rules that pores and throats are treated the same
    IP returns smallest Pc first IP for imbibition returns largest

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

    def run(self, phase, inlets, throat_prop='throat.capillary_pressure',
            pore_prop='pore.capillary_pressure', imbibition=True):
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
        p_entry = phase[pore_prop]
        # switch algorithm to sort highest first
        t_seq = sp.zeros(len(t_entry), dtype=int)
        p_seq = sp.zeros(len(p_entry), dtype=int)
        t_inv_Pc = sp.zeros(len(t_entry))
        p_inv_Pc = sp.zeros(len(p_entry))
        # Perform initial analysis on input pores
        p_seq[inlets] = 1
        Ts = net.find_neighbor_throats(pores=inlets)
        [hq.heappush(queue, (t_entry[T], T, 't')) for T in Ts]
        sequence = 1
        inv_Pc = -sp.inf
        while len(queue) > 0:
            # Find element at the top of the queue
            (Pc, index, data_type) = hq.heappop(queue)
            if Pc > inv_Pc:
                inv_Pc = Pc
            # Extract actual throat number
            if data_type == 't':
                t_seq[index] = sequence
                t_inv_Pc[index] = inv_Pc
                sequence += 1
                newPs = net['throat.conns'][index]
                for P in newPs:
                    if p_seq[P] == 0:
                        hq.heappush(queue, (p_entry[P], P, 'p'))
            elif data_type == 'p':
                p_seq[index] = sequence
                p_inv_Pc[index] = inv_Pc
                sequence += 1
                newTs = net.find_neighbor_throats(pores=index)
                for T in newTs:
                    if t_seq[T] == 0:
                        hq.heappush(queue, (t_entry[T], T, 't'))

        self['throat.invasion_sequence'] = t_seq
        self['pore.invasion_sequence'] = p_seq
        self['throat.invasion_Pc'] = t_inv_Pc
        self['pore.invasion_Pc'] = p_inv_Pc

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
        self._phase['throat.invasion_Pc'] = sp.nan
        self._phase['pore.invasion_Pc'] = sp.nan
        self._phase['throat.invasion_Pc'][throats] = \
            self['throat.invasion_Pc'][throats]
        self._phase['pore.invasion_Pc'][pores] = \
            self['pore.invasion_Pc'][pores]

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

    def evaluate_trapping(self, outlets):
        r"""
        Finds trapped pores and throats after a full ordinary
        percolation drainage has been run
        Parameters
        ----------
        outlets : array_like
            A list of pores that define the wetting phase outlets.
            Disconnection from these outlets results in trapping.
        """
        try:
            inv_points = sp.unique(self['pore.invasion_Pc'])  # Get points used in OP
        except:
            logger.error('Invasion percolation has not been run!')
            raise Exception('Aborting algorithm')
        self._p_trap = sp.zeros_like(self['pore.invasion_Pc'], dtype=float)
        self._t_trap = sp.zeros_like(self['throat.invasion_Pc'], dtype=float)
        tind = self._net.throats()
        conns = self._net.find_connected_pores(tind)
        for inv_val in inv_points[0:-1]:
            # Find clusters of defender pores
            Pinvaded = self['pore.invasion_Pc'] <= inv_val
            Cstate = sp.sum(Pinvaded[conns], axis=1)
            Tinvaded = self['throat.invasion_Pc'] <= inv_val
            Cstate = Cstate + Tinvaded
            # 0 = all open, 1=1 pore filled
            # 2=2 pores filled
            # 3=2 pores + 1 throat filled
            clusters = self._net.find_clusters(Cstate == 0)
            # Clean up clusters (invaded = -1, defended >=0)
            clusters = clusters*(~Pinvaded) - (Pinvaded)
            # Identify clusters connected to outlet sites
            out_clusters = sp.unique(clusters[outlets])
            trapped_pores = ~sp.in1d(clusters, out_clusters)
            self._p_trap[(self._p_trap == 0)[trapped_pores]] = inv_val
            trapped_throats = self._net.find_neighbor_throats(trapped_pores,
                                                              mode='intersection')
            if len(trapped_throats) > 0:
                self._t_trap[(self._t_trap == 0)[trapped_throats]] = inv_val
        self['pore.invasion_Pc'][self._p_trap > 0] = sp.inf
        self['throat.invasion_Pc'][self._t_trap > 0] = sp.inf
