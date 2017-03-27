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

    def set_occupancy(self, sequence):
        self._phase['throat.occupancy'] = self['throat.invasion_sequence'] <= sequence
        self._phase['pore.occupancy'] = self['pore.invasion_sequence'] <= sequence
    
    
    def apply_trapping(self, outlets, bt=False):
        """
        Apply trapping based on algorithm described by Y. Masson [1].
        It is applied as a post-process and runs the percolation algorithm in 
        reverse assessing the occupancy of un-invaded pore neighbors.
        3 situations can happen:
            The number of clusters stays the same
            A new cluster of size one is created
            Multiple clusters are merged together
        Ref:
        [1] Masson, Y., 2016. A fast two-step algorithm for invasion
        percolation with trapping. Computers & Geosciences, 90, pp.41-48
        """
        if bt:
            #Go from breakthrough
            #First assess sequence at which break-through was acheived
            bt_seq = np.min(self['pore.invasion_sequence'][outlets])
            print("Break-through Sequence: ",bt_seq)
            #Set occupancy
            self.set_occupancy(bt_seq)
            #Put defending phase into clusters
            clusters = self._net.find_clusters2(~self._phase['pore.occupancy'])
            #Identify clusters that are connected to an outlet and set to -2
            #-1 is the invaded fluid
            #-2 is the defender fluid able to escape
            #All others now trapped clusters which grow as invasion is reversed
            out_clusters = sp.unique(clusters[outlets])
            for c in out_clusters:
                if c >=0:
                    clusters[clusters==c] = -2
        else:
            #Go from end
            clusters = np.ones(self._net.Np)*-1
            clusters[outlets] = -2
            bt_seq = np.max(self['pore.invasion_sequence'])

        inv_list = list(self['pore.invasion_sequence'])
        
        next_cluster_num = np.max(clusters)+1
        #For all the steps after the inlets are set up to break-through
        #Reverse the sequence and assess the neighbors cluster state
        for uninvasion_sequence in np.arange(1, bt_seq+1)[::-1]:
            if uninvasion_sequence in inv_list:
                pore = inv_list.index(uninvasion_sequence)
                neighbors = self._net.find_neighbor_pores(pore)
                if np.all(clusters[neighbors] == -1) :
                    #This is the start of a new trapped cluster
                    clusters[pore] = next_cluster_num
                    next_cluster_num +=1
                    logger.info("C: 1, P: "+str(pore)+
                                " new cluster number: "+str(clusters[pore]))
                elif (np.all(clusters[neighbors] == clusters[neighbors][0]) and
                             clusters[neighbors][0] > -1):
                    #This means pore belongs to this cluster
                    clusters[pore] = clusters[neighbors][0]
                    logger.info("C: 2, P: "+str(pore)+
                                " joins cluster number: "+str(clusters[pore]))
                else:
                    #There are a mixture of neighboring clusters so merge them
                    new_num = None
                    nc = np.unique(clusters[neighbors])
                    for c in nc:
                        if c >= 0:
                            if new_num == None:
                                new_num = c
                            else:
                                clusters[clusters == c] = new_num
                                logger.info("C: 3, P: "+str(pore)+
                                            " merge clusters: "+str(c)+" into "+
                                            str(new_num))
                            clusters[pore] = new_num
                #Now check whether a neighbor is connected to a sink
                if -2 in clusters[neighbors]:
                    #Whoopie we found an outlet so can escape
                    logger.info("C: 4, P: "+str(pore)+ " can escape")
                    clusters[pore] = -2

                
        #And now return clusters
        self['pore.clusters']=clusters
        logger.info("Number of trapped clusters",
                    np.sum(np.unique(clusters)>=0))
