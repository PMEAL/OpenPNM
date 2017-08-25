# -*- coding: utf-8 -*-
"""
===============================================================================
InvasionPercolationBasic: Simple IP
===============================================================================

"""
import heapq as hq
import scipy as sp
import numpy as np
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

    def setup(self, phase, throat_prop='throat.capillary_pressure', **kwargs):
        r"""
        Set up the required parameters for the algorithm

        Parameters
        ----------
        phase : OpenPNM Phase object
            The phase to be injected into the Network.  The Phase must have the
            capillary entry pressure values for the system.

        throat_prop : string
            The name of the throat property containing the capillary entry
            pressure.  The default is 'throat.capillary_pressure'.

        """
        self._phase = phase
        # Setup arrays and info
        self['throat.entry_pressure'] = phase[throat_prop]
        # Indices into t_entry giving a sorted list
        self['throat.sorted'] = sp.argsort(self['throat.entry_pressure'], axis=0)
        self['throat.order'] = sp.zeros_like(self['throat.sorted'])
        self['throat.order'][self['throat.sorted']] = sp.arange(0, self._net.Nt)
        self['throat.invaded'] = -sp.ones((self._net.Nt,))
        self['pore.invaded'] = -sp.ones((self._net.Np,))
        self._tcount = 0

    def set_inlets(self, pores=None, **kwargs):
        r"""

        Parameters
        ----------
        pores : array_like
            The list of inlet pores from which the Phase can enter the Network
        """
        if 'inlets' in kwargs.keys():
            pores = kwargs['inlets']
        self['pore.invaded'][pores] = 0

        # Perform initial analysis on input pores
        Ts = self._net.find_neighbor_throats(pores=pores)
        self.queue = []
        [hq.heappush(self.queue, T) for T in self['throat.order'][Ts]]

    def run(self, n_steps=None, **kwargs):
        r"""
        Perform the algorithm

        Parameters
        ----------
        n_steps : int
            The number of throats to invaded during this step

        """
        if 'throat.entry_pressure' not in self.keys():
            self.setup(**kwargs)
        if sp.all(self['pore.invaded'] == -1):
            self.set_inlets(**kwargs)

        if n_steps is None:
            n_steps = sp.inf

        queue = self.queue
        if len(queue) == 0:
            logger.warn('queue is empty, this network is fully invaded')
            return
        t_sorted = self['throat.sorted']
        t_order = self['throat.order']
        t_inv = self['throat.invaded']
        p_inv = self['pore.invaded']

        count = 0
        while (len(queue) > 0) and (count < n_steps):
            # Find throat at the top of the queue
            t = hq.heappop(queue)
            # Extract actual throat number
            t_next = t_sorted[t]
            t_inv[t_next] = self._tcount
            # If throat is duplicated
            while len(queue) > 0 and queue[0] == t:
                # Note: Preventing duplicate entries below might save some time here
                t = hq.heappop(queue)
            # Find pores connected to newly invaded throat
            Ps = self._net['throat.conns'][t_next]
            # Remove already invaded pores from Ps
            Ps = Ps[p_inv[Ps] < 0]
            if len(Ps) > 0:
                p_inv[Ps] = self._tcount
                Ts = self._net.find_neighbor_throats(pores=Ps)
                Ts = Ts[t_inv[Ts] < 0]  # Remove already invaded throats from Ts
                [hq.heappush(queue, T) for T in t_order[Ts]]
            count += 1
            self._tcount += 1
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

    def apply_trapping(self, outlets):
        """
        Apply trapping based on algorithm described by Y. Masson [1].
        It is applied as a post-process and runs the percolation algorithm in
        reverse assessing the occupancy of pore neighbors. Consider the
        following scenario when running standard IP without trapping,
        3 situations can happen after each invasion step:
            The number of defending clusters stays the same and clusters can
            shrink
            A cluster of size one is suppressed
            A cluster is split into multiple clusters
        In reverse the following opposite situations can happen:
            The number of defending clusters stays the same and clusters can
            grow
            A cluster of size one is created
            Mutliple clusters merge into one cluster
        With trapping the reversed rules are adjusted so that:
            Only clusters that do not connect to a sink can grow and merge.
            At the point that a neighbor connected to a sink is touched the
            trapped cluster stops growing as this is the point of trapping in
            forward invasion time.

        Logger info displays the invasion sequence and pore index and a message
        with condition number based on the modified trapping rules and the
        assignment of the pore to a given cluster.

        Initially all invaded pores are given cluster label -1
        Outlets / Sinks are given -2
        New clusters that grow into fully trapped clusters are either
        identified at the point of breakthrough or grow from nothing if the
        full invasion sequence is run, they are assigned numbers from 0 up.

        Ref:
        [1] Masson, Y., 2016. A fast two-step algorithm for invasion
        percolation with trapping. Computers & Geosciences, 90, pp.41-48

        Parameters
        ----------
        outlets : list or array of pore indices for defending fluid to escape
        through

        Returns
        -------
        Creates a throat array called 'pore.clusters' in the Algorithm
        dictionary. Any positive number is a trapped cluster

        Also creates 2 boolean arrays Np and Nt long called '<element>.trapped'
        """
        # First see if network is fully invaded
        invaded_ps = self['pore.invasion_sequence'] > -1
        if ~np.all(invaded_ps):
            # Put defending phase into clusters
            clusters = self._net.find_clusters2(~invaded_ps)
            # Identify clusters that are connected to an outlet and set to -2
            # -1 is the invaded fluid
            # -2 is the defender fluid able to escape
            # All others now trapped clusters which grow as invasion is reversed
            out_clusters = sp.unique(clusters[outlets])
            for c in out_clusters:
                if c >= 0:
                    clusters[clusters == c] = -2
        else:
            # Go from end
            clusters = np.ones(self._net.Np, dtype=int)*-1
            clusters[outlets] = -2

        # Turn into a list for indexing
        inv_seq = np.vstack((self['pore.invasion_sequence'].astype(int),
                             np.arange(0, self._net.Np, dtype=int))).T
        # Reverse sort list
        inv_seq = inv_seq[inv_seq[:, 0].argsort()][::-1]
        next_cluster_num = np.max(clusters)+1
        # For all the steps after the inlets are set up to break-through
        # Reverse the sequence and assess the neighbors cluster state
        stopped_clusters = np.zeros(self._net.Np, dtype=bool)
        all_neighbors = self._net.find_neighbor_pores(self._net.pores(),
                                                      flatten=False)
        for un_seq, pore in inv_seq:
            if pore not in outlets:  # Don't bother with outlets
                nc = clusters[all_neighbors[pore]]  # Neighboring clusters
                unique_ns = np.unique(nc[nc != -1])  # Unique Neighbors
                seq_pore = "S:"+str(un_seq)+" P:"+str(pore)
                if np.all(nc == -1):
                    # This is the start of a new trapped cluster
                    clusters[pore] = next_cluster_num
                    next_cluster_num += 1
                    msg = (seq_pore+" C:1 new cluster number: " +
                           str(clusters[pore]))
                    logger.info(msg)
                elif len(unique_ns) == 1:
                    # Grow the only connected neighboring cluster
                    if not stopped_clusters[unique_ns[0]]:
                        clusters[pore] = unique_ns[0]
                        msg = (seq_pore+" C:2 joins cluster number: " +
                               str(clusters[pore]))
                        logger.info(msg)
                    else:
                        clusters[pore] = -2
                elif -2 in unique_ns:
                    # We have reached a sink neighbor, stop growing cluster
                    msg = (seq_pore+" C:3 joins sink cluster")
                    logger.info(msg)
                    clusters[pore] = -2
                    # Stop growth and merging
                    stopped_clusters[unique_ns[unique_ns > -1]] = True
                else:
                    # We might be able to do some merging
                    # Check if any stopped clusters are neighbors
                    if np.any(stopped_clusters[unique_ns]):
                        msg = (seq_pore+" C:4 joins sink cluster")
                        logger.info(msg)
                        clusters[pore] = -2
                        # Stop growing all neighboring clusters
                        stopped_clusters[unique_ns] = True
                    else:
                        # Merge multiple un-stopped trapped clusters
                        new_num = unique_ns[0]
                        clusters[pore] = new_num
                        for c in unique_ns:
                            clusters[clusters == c] = new_num
                            msg = (seq_pore + " C:5 merge clusters: " +
                                   str(c) + " into "+str(new_num))
                            logger.info(msg)

        # And now return clusters
        self['pore.clusters'] = clusters
        logger.info("Number of trapped clusters" +
                    str(np.sum(np.unique(clusters) >= 0)))
        self['pore.trapped'] = self['pore.clusters'] > -1
        trapped_ts = self._net.find_neighbor_throats(self['pore.trapped'])
        self['throat.trapped'] = np.zeros([self._net.Nt], dtype=bool)
        self['throat.trapped'][trapped_ts] = True
        self['pore.invasion_sequence'][self['pore.trapped']] = np.inf
        self['throat.invasion_sequence'][self['throat.trapped']] = np.inf
