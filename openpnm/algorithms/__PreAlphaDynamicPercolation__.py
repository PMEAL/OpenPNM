# -*- coding: utf-8 -*-
"""
===============================================================================
InvasionPercolationBasic: Simple IP
===============================================================================

"""
import heapq as hq
import scipy as sp
import numpy as np
from openpnm.algorithms import GenericAlgorithm
import logging
import matplotlib.pyplot as plt
import time
logger = logging.getLogger(__name__)


class MixedPercolation(GenericAlgorithm):
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

    def setup(self, phase, def_phase, **kwargs):
        r"""
        Set up the required parameters for the algorithm

        Parameters
        ----------
        phase : OpenPNM Phase object
            The phase to be injected into the Network.  The Phase must have the
            capillary entry pressure values for the system.

        prop : string
            The name of the property containing the capillary entry
            pressure.  This should be applied to both throats and pores.

        """
        self._phase = phase
        self._def = def_phase
        # Setup arrays and info
        if np.shape(np.shape(phase['throat.capillary_pressure']))[0] > 1:
            self._bi_directional = True
            # If bi-directional throats, get the first one for now, switched
            # later when algorithm runs and works out which one to apply
            tcp = phase['throat.capillary_pressure']
            self['throat.entry_pressure'] = tcp[:, 0]
        else:
            self['throat.entry_pressure'] = phase['throat.capillary_pressure']
            self._bi_directional = False
        self['pore.entry_pressure'] = phase['pore.capillary_pressure']
        self.reset_invasion_info()
        self._key_words = kwargs
        # Need to setup cooperative pore filling seperately
        self._coop_fill = False

    def reset_invasion_info(self):
        r'''
        Resets all the invasion data
        '''
        for elem in ['pore', 'throat']:
            for prop in ['inv_Pc',
                         'inv_sat',
                         'occupancy',
                         'action']:
                try:
                    del self[elem+'.'+prop]
                except:
                    pass
            for prop in ['inv_seq',
                         'inv_Pc',
                         'cluster',
                         'action']:
                try:
                    self[elem+'.'+prop] = -1
                except:
                    pass
        # Masks for tracking pores and throats at the interface
        # Basically a quick way of getting to all the elements in the queues
        self._interface_Ts = np.zeros(self.Nt, dtype=bool)
        self._interface_Ps = np.zeros(self.Np, dtype=bool)

    def set_inlets(self, pores=None, clusters=None):
        r"""

        Parameters
        ----------
        pores : array_like
            The list of inlet pores from which the Phase can enter the Network
        """
        if 'inlets' in self._key_words.keys():
            pores = self._key_words['inlets']

        if clusters is None:
            clusters = []
            clusters.append(pores)
        self.queue = {}
        try:
            inlet_inv_seq = self._key_words['inlet_inv_seq']
        except:
            inlet_inv_seq = -1
        for i, cluster in enumerate(clusters):
            self.queue[i] = []
            # Perform initial analysis on input pores
            self['pore.inv_seq'][cluster] = inlet_inv_seq
            self['pore.cluster'][cluster] = i
            if np.size(cluster) > 1:
                for elem_id in cluster:
                    self._add_ts2q(elem_id, self.queue[i], action=0)
            elif np.size(cluster) == 1:
                self._add_ts2q(cluster, self.queue[i], action=0)
            else:
                logger.warning("Some inlet clusters have no pores")
        try:
            if self._key_words['snap_off']:
                self.apply_snap_off()
        except:
            pass
        try:
            if self._key_words['partial']:
                self.apply_partial_sat()
        except:
            self._phase['pore.occuancy'] = False
            self._phase['throat.occupancy'] = False
            self._def['pore.occupancy'] = True
            self._def['throat.occupancy'] = True

    def _add_ts2q(self, pore, queue, action=-1):
        """
        Helper method to add throats to the queue
        """
        net = self.project.network
        elem_type = 'throat'
        # Find throats connected to newly invaded pore
        Ts = net.find_neighbor_throats(pores=pore)
        # Remove already invaded throats from Ts
        Ts = Ts[self['throat.inv_seq'][Ts] <= 0]
        tcp = self._phase['throat.capillary_pressure']
        if len(Ts) > 0:
            self._interface_Ts[Ts] = True
            for T in Ts:
                if self._bi_directional:
                    # Get index of pore being invaded next and apply correct
                    # entry pressure
                    pmap = net['throat.conns'][T] != pore
                    pind = list(pmap).index(True)
                    self['throat.entry_pressure'][T] = tcp[T][pind]
                data = []
                # Pc
                data.append(self['throat.entry_pressure'][T])
                # Element Index
                data.append(T)
                # Element Type (Pore of Throat)
                data.append(elem_type)
                # Invasion Action - 0=Burst, 1=Coop, 2=Snap, 3=Touch??
                data.append(action)
                hq.heappush(queue, data)

    def _add_ps2q(self, throat, queue, action=-1):
        """
        Helper method to add pores to the queue
        """
        net = self.project.network
        elem_type = 'pore'
        # Find pores connected to newly invaded throat
        Ps = net['throat.conns'][throat]
        # Remove already invaded pores from Ps
        Ps = Ps[self['pore.inv_seq'][Ps] <= 0]
        if len(Ps) > 0:
            self._interface_Ps[Ps] = True
            for P in Ps:
                data = []
                # Pc
                data.append(self["pore.entry_pressure"][P])
                # Element Index
                data.append(P)
                # Element Type (Pore of Throat)
                data.append(elem_type)
                # Invasion Action - 0=Burst, 1=Coop, 2=Snap, 3=Touch??
                data.append(action)
                hq.heappush(queue, data)

    def run(self, max_pressure=None, outlets=None, **kwargs):
        r"""
        Perform the algorithm

        Parameters
        ----------
        n_steps : int
            The number of throats to invaded during this step

        """
        if 'throat.entry_pressure' not in self.keys():
            logger.error("Setup method must be run first")
        if 'inlets' in self._key_words.keys():
            logger.info("Setting inlet pores at shared pressure")
            self.set_inlets(pores=self._key_words['inlets'])
        elif 'clusters' in self._key_words.keys():
            logger.info("Setting inlet clusters at individual pressures")
            self.set_inlets(clusters=self._key_words['clusters'])
        else:
            logger.error("Either 'inlets' or 'clusters' must be passed to" +
                         " setup method")

        if max_pressure is None:
            max_pressure = sp.inf
        if len(self.queue.items()) == 0:
            logger.warn('queue is empty, this network is fully invaded')
            return

        max_p_reached = [False]*len(self.queue.items())
        count = -1
        invasion_running = [True]*len(self.queue.items())
        high_Pc = np.ones(len(self.queue.items()))*-np.inf
        while np.any(invasion_running) and not np.all(max_p_reached):
            # Loop over clusters
            for c_num in self.queue.keys():
                if invasion_running[c_num]:
                    queue = self.queue[c_num]
                    pressure, elem_id, elem_type, action = hq.heappop(queue)
                    if elem_type == 'pore':
                        self._interface_Ps[elem_id] = False
                    else:
                        self._interface_Ts[elem_id] = False
                    if pressure > max_pressure:
                        max_p_reached[c_num] = True
                    else:
                        elem_cluster = self[elem_type+'.cluster'][elem_id]
                        elem_cluster = elem_cluster.astype(int)
                        # Cluster is the uninvaded cluster
                        if elem_cluster == -1:
                            count += 1
                            # Record highest Pc cluster has reached
                            if high_Pc[c_num] < pressure:
                                high_Pc[c_num] = pressure
                            # The newly invaded element is available for
                            # invasion
                            self[elem_type+'.inv_seq'][elem_id] = count
                            self[elem_type+'.cluster'][elem_id] = c_num
                            self[elem_type+'.inv_Pc'][elem_id] = high_Pc[c_num]
                            self[elem_type+'.action'][elem_id] = action
                            if elem_type == 'throat':
                                self._add_ps2q(elem_id, queue, action=0)
                            elif elem_type == 'pore':
                                self._add_ts2q(elem_id, queue, action=0)
                                if self._coop_fill:
                                    self._check_coop(elem_id, queue)
                        # Element is part of existing cluster
                        elif (elem_cluster != c_num and
                              invasion_running[elem_cluster]):
                            # The newly invaded element is part of an invading
                            # cluster. Merge the clusters using the existing
                            # cluster number
                            logger.info("Merging cluster "+str(c_num) +
                                        " into cluster "+str(elem_cluster) +
                                        " at sequence "+str(count))

                        elif (elem_cluster != c_num and not
                              invasion_running[elem_cluster]):
                            # The newly invaded element is part of cluster that
                            # has stopped invading
                            logger.info("Cluster " + str(c_num) +
                                        " terminated")

                        elif self[elem_type+'.cluster'][elem_id] == c_num:
                            # Self intersecting or repeating elements"
                            pass
                        else:
                            logger.warning("Clusters " + str(c_num) + " and " +
                                           str(elem_cluster) + " performed " +
                                           " strange operation!")

                    if len(queue) == 0 or max_p_reached[c_num]:
                        # If the cluster contains no more entries invasion has
                        # finished
                        invasion_running[c_num] = False
            self._invade_isolated_Ts()
            if outlets is not None:
                # terminated clusters
                tcs = np.unique(self['pore.cluster'][outlets]).astype(int)
                tcs = tcs[tcs >= 0]
                if len(tcs) > 0:
                    for tc in tcs:
                        if invasion_running[tc] is True:
                            invasion_running[tc] = False
                            logger.info("Cluster " + str(tc) + " reached " +
                                        " outlet at sequence " + str(count))

    def return_results(self, pores=[], throats=[], Pc=None):
        r"""
        Places the results of the IP simulation into the Phase object.

        Parameters
        ----------
        pores and throats : array_like
            The list of pores and throats whose values should be returned to
            the Phase object.  Default is all of them.

        """
        pores = sp.array(pores, ndmin=1)
        throats = sp.array(throats, ndmin=1)
        if len(pores) == 0:
            pores = self.Ps
        if len(throats) == 0:
            throats = self.Ts
        self._phase['throat.inv_seq'] = sp.nan
        self._phase['pore.inv_seq'] = sp.nan
        self._phase['throat.inv_seq'][throats] = \
            self['throat.inv_seq'][throats]
        self._phase['pore.inv_seq'][pores] = \
            self['pore.inv_seq'][pores]
        self._phase['throat.cluster'] = self['throat.cluster']
        self._phase['pore.cluster'] = self['pore.cluster']
        if Pc is not None:
            if 'pore.inv_Pc' in self.keys():
                self['throat.occupancy'] = self['throat.inv_Pc'] <= Pc
                self['pore.occupancy'] = self['pore.inv_Pc'] <= Pc
                self._phase['throat.occupancy'] = self['throat.occupancy']
                self._phase['pore.occupancy'] = self['pore.occupancy']
                self._def['throat.occupancy'] = ~self['throat.occupancy']
                self._def['pore.occupancy'] = ~self['pore.occupancy']
            else:
                logger.warning("Occupancy not updated, please run " +
                               "extract_drainage() method to populate" +
                               " the inv_Pc data")
        if "pore.inv_Pc" in self.props():
            self._phase['pore.inv_Pc'] = \
                self['pore.inv_Pc']
            self._phase['throat.inv_Pc'] = \
                self['throat.inv_Pc']
        if "pore.inv_sat" in self.props():
            self._phase['pore.inv_sat'] = \
                self['pore.inv_sat']
            self._phase['throat.inv_sat'] = \
                self['throat.inv_sat']
        if "pore.trapped" in self.labels():
            self._phase['pore.trapped'] = self['pore.trapped']
        if "pore.trapping_sequence" in self.props():
            self._phase['pore.trapping_sequence'] = \
                self['pore.trapping_sequence']
        if "throat.trapped" in self.labels():
            self._phase['throat.trapped'] = self['throat.trapped']
        if "throat.trapping_sequence" in self.props():
            self._phase['throat.trapping_sequence'] = \
                self['throat.trapping_sequence']

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
        net = self.project.network
        P12 = net['throat.conns']
        aa = self['throat.inv_seq']
        bb = sp.argsort(self['throat.inv_seq'])
        P12_inv = self['pore.inv_seq'][P12]
        # Find if the connected pores were invaded with or before each throat
        P1_inv = P12_inv[:, 0] <= aa
        P2_inv = P12_inv[:, 1] <= aa
        cc = sp.column_stack((P1_inv, P2_inv))
        # List of Pores invaded with each throat
        dd = sp.sum(cc, axis=1, dtype=bool)
        # Find volume of these pores
        P12_vol = sp.sum(net['pore.volume'][P12]*cc, axis=1)*dd
        # Add invaded throat volume to pore volume (if invaded)
        T_vol = P12_vol + net['throat.volume']
        # Cumulative sum on the sorted throats gives cumulated inject volume
        ee = sp.cumsum(T_vol[bb] / flowrate)
        t = sp.zeros((self.Nt,))
        t[bb] = ee  # Convert back to original order
        self._phase['throat.invasion_time'] = t

    def plot_drainage_curve(self, fig=None, inv_points=None, npts=100,
                            lpf=False, trapping_outlets=None):
        r"""
        Plot a simple drainage curve
        """
        net = self.project.network
        if "pore.inv_Pc" not in self.props():
            logger.error("Cannot plot drainage curve. Please run " +
                         " algorithm first")
        if inv_points is None:
            ok_Pc = self['throat.inv_Pc'][~sp.isnan(self['throat.inv_Pc'])]
            inv_points = np.unique(ok_Pc)
        sat_p = np.zeros(len(inv_points))
        sat_t = np.zeros(len(inv_points))
        inv_p = self['pore.inv_Pc']
        inv_t = self['throat.inv_Pc']
        # Handle trapped pores and throats
        if np.sum(self['pore.inv_seq'] == -1) > 0:
            inv_p[self['pore.inv_seq'] == -1] = 2*inv_points.max()
        if np.sum(self['throat.inv_seq'] == -1) > 0:
            inv_t[self['throat.inv_seq'] == -1] = 2*inv_points.max()
        num_p = np.zeros(len(inv_points), dtype=int)
        num_t = np.zeros(len(inv_points), dtype=int)
        for i, Pc in enumerate(inv_points):
            if lpf:
                frac = self.evaluate_late_pore_filling(Pc,
                                                       Swp_init=0.25,
                                                       eta=2.5)
                p_vol = net['pore.volume']*frac
            else:
                p_vol = net['pore.volume']
            sat_p[i] = np.sum(p_vol[inv_p <= Pc])
            sat_t[i] = np.sum(net['throat.volume'][inv_t <= Pc])
            num_p[i] = np.sum(inv_p <= Pc)
            num_t[i] = np.sum(inv_t <= Pc)

        pvol = np.sum(net['pore.volume'])
        tvol = np.sum(net['throat.volume'])
        tot_vol = pvol + tvol
        tot_sat = sat_p + sat_t
        # Normalize
        sat_p /= tot_vol
        sat_t /= tot_vol
        tot_sat /= tot_vol
        if fig is None:
            fig = plt.figure()
        a = fig.add_subplot(111)
        a.plot(inv_points, sat_p, 'r*-', label='pore')
        a.plot(inv_points, sat_t, 'b*-', label='throat')
        a.plot(inv_points, tot_sat, 'g*-', label='total')
        a.legend(bbox_to_anchor=(0, 1.02, 1, 0.102),
                 loc=3, ncol=3, borderaxespad=0)
        a.set_xlabel("Capillary Pressure [Pa]")
        a.set_ylabel("Saturation")
        return fig, tot_sat, num_p, num_t

    def evaluate_late_pore_filling(self, Pc, Swp_init=0.75, eta=3.0,
                                   wetting_phase=False):
        r"""
        Compute the volume fraction of the phase in each pore given an initial
        wetting phase fraction (Swp_init) and a growth exponent (eta)
        returns the fraction of the pore volume occupied by wetting or
        non-wetting phase.
        Assumes Non-wetting phase displaces wetting phase
        """
        net = self.project.network
        try:
            Swp = np.ones(len(self['pore.inv_Pc']))
            mask = self['pore.inv_Pc'] < Pc
            Swp[mask] = Swp_init*(self['pore.inv_Pc'][mask]/Pc)**eta
            Swp[np.isnan(Swp)] = 1.0
        except:
            Swp = np.ones(net.Np)
        Snwp = 1-Swp
        if wetting_phase:
            return Swp
        else:
            return Snwp

    def apply_trapping(self, outlets, partial=False):
        """
        Apply trapping based on algorithm described by Y. Masson [1].
        It is applied as a post-process and runs the percolation algorithm in
        reverse assessing the occupancy of pore neighbors.
        3 situations can happen on invasion without trapping:
            The number of defending clusters stays the same and clusters can
            shrink
            A cluster of size one is suppressed
            A cluster is split into multiple clusters
        In reverse the following situations can happen:
            The number of defending clusters stays the same and clusters can
            grow
            A cluster of size one is created
            Mutliple clusters merge into one cluster
        With trapping the reversed rules are adjusted so that:
            Only clusters that do not connect to a sink can grow and merge.
            At the point that a neighbor connected to a sink is touched the
            trapped cluster stops growing as this is the point of trapping in
            forward invasion time.

        Logger info displays the invasion Sequence and pore index and a message
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

        partial : Boolean indicating whether partially filled network

        Returns
        -------
        Creates a throat array called 'pore.clusters' in the Algorithm
        dictionary. Any positive number is a trapped cluster

        Also creates 2 boolean arrays Np and Nt long called '<element>.trapped'
        """
        net = self.project.network
        if partial:
            # Set occupancy
            invaded_ps = self['pore.inv_seq'] > -1
            # Put defending phase into clusters
            clusters = net.find_clusters2(~invaded_ps)
            # Identify clusters that are connected to an outlet and set to -2
            # -1 is the invaded fluid
            # -2 is the defender fluid able to escape
            # All others now trapped clusters which grow as invasion is
            # reversed
            out_clusters = sp.unique(clusters[outlets])
            for c in out_clusters:
                if c >= 0:
                    clusters[clusters == c] = -2
        else:
            # Go from end
            clusters = np.ones(net.Np, dtype=int)*-1
            clusters[outlets] = -2

        # Turn into a list for indexing
        inv_seq = np.vstack((self['pore.inv_seq'].astype(int),
                             np.arange(0, net.Np, dtype=int))).T
        # Reverse sort list
        inv_seq = inv_seq[inv_seq[:, 0].argsort()][::-1]
        next_cluster_num = np.max(clusters)+1
        # For all the steps after the inlets are set up to break-through
        # Reverse the sequence and assess the neighbors cluster state
        stopped_clusters = np.zeros(net.Np, dtype=bool)
        all_neighbors = net.find_neighbor_pores(net.pores(), flatten=False)
        for un_seq, pore in inv_seq:
            if pore not in outlets and un_seq > -1:  # Don't include outlets
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
        clusters[outlets] = -2
        num_trap = np.sum(np.unique(clusters) >= 0)
        if num_trap > 0:
            logger.info("Number of trapped clusters" + str(num_trap))
            self['pore.trapped'] = clusters > -1
            num_tPs = np.sum(self['pore.trapped'])
            logger.info("Number of trapped pores: " + str(num_tPs))
            self['pore.inv_seq'][self['pore.trapped']] = -1
            self['throat.trapped'] = np.zeros([net.Nt], dtype=bool)
            for c in np.unique(clusters[clusters >= 0]):
                c_ts = net.find_neighbor_throats(clusters == c,
                                                 mode='intersection')
                self['throat.trapped'][c_ts] = True
            num_tTs = np.sum(self['throat.trapped'])
            logger.info("Number of trapped throats: " + str(num_tTs))
            self['throat.inv_seq'][self['throat.trapped']] = -1
            # Assumes invasion has run to the end
            self._phase['pore.occupancy'] = ~self['pore.trapped']
            self._def['pore.occupancy'] = self['pore.trapped']
            self._phase['throat.occupancy'] = ~self['throat.trapped']
            self._def['throat.occupancy'] = self['throat.trapped']
        else:
            logger.info("No trapped clusters found")

    def apply_snap_off(self, snap_off='throat.snap_off', queue=None):
        r"""
        Add all the throats to the queue with snap off pressure
        """
        net = self.project.network
        if queue is None:
            queue = self.queue[0]
        try:
            Pc_snap_off = self._phase[snap_off]
            logger.info("Adding snap off pressures to queue")
            for T in net.throats():
                if not np.isnan(Pc_snap_off[T]):
                    hq.heappush(queue, [Pc_snap_off[T], T, 'throat', 2])
        except:
            logger.warning("Phase " + self._phase.name + " doesn't have " +
                           "property " + snap_off)

    def apply_partial_sat(self, queue=None):
        r"""
        Method to start invasion from a partially saturated state
        """
        net = self.project.network
        if queue is None:
            invading_cluster = 0
            queue = self.queue[invading_cluster]
        occ_type = self._phase['pore.occupancy'].dtype
        occupied = np.array([1], dtype=occ_type)
        occ_Ps = self._phase['pore.occupancy'] == occupied
        occ_Ts = self._phase['throat.occupancy'] == occupied
        low_val = -np.inf
        if np.sum(occ_Ps) > 0:
            logger.warn("Applying partial saturation to " +
                        str(np.sum(occ_Ps)) + " pores")
            self['pore.inv_seq'][occ_Ps] = 0
            for P in net.pores()[occ_Ps]:
                self._add_ts2q(P, queue, action=0)
                self['pore.cluster'][P] = invading_cluster
                self['pore.inv_Pc'][P] = low_val
        if np.sum(occ_Ts) > 0:
            logger.warn("Applying partial saturation to " +
                        str(np.sum(occ_Ts)) + " throats")
        self['throat.inv_seq'][occ_Ts] = 0
        for T in net.throats()[occ_Ts]:
            self['throat.cluster'][T] = invading_cluster
            self['throat.inv_Pc'][T] = low_val

    def _perpendicular_vector(self, vec):
        return np.cross(vec, [1, 1, 1])

    def _circ_points(self, a, b, c):
        t = np.arange(0, 2*np.pi, 0.1)
        pts = []
        for i in range(3):
            pts.append(c[i] + a[i]*np.sin(t) + b[i]*np.cos(t))
        return np.asarray(pts).T

    def _check_intersection(self, c1, c2, r1, r2, pore_center, pore_rad):
        r"""
        Helper function to determine whether the intersetion between two
        Spheres formed by growing menisci in the throats lies within the pore
        body that they connect to
        """
        intersection = np.zeros(len(r1), dtype=bool)
        vec_n = c2 - c1
        dist = np.linalg.norm(vec_n, axis=1)
        vec_a = self._perpendicular_vector(vec_n)
        dist_a = np.linalg.norm(vec_a, axis=1)
        vec_a *= 1/(np.vstack((dist_a, dist_a, dist_a)).T)
        vec_b = np.cross(vec_a, vec_n)
        dist_b = np.linalg.norm(vec_b, axis=1)
        vec_b *= 1/(np.vstack((dist_b, dist_b, dist_b)).T)
        x = (dist**2 - r2**2 + r1**2)/(2*dist)
        # intersection centre
        p = c1 + (vec_n.T*(x/dist)).T
        sq = 4 * dist**2 * r1**2 - (dist**2 - r2**2 + r1**2)**2
        # Could try to vectorize this but it would be pretty complicated!!!
        for i in range(len(r1)):
            if sq[i] > 0:
                # An intersection is found
                rad = (1/(2*dist[i]))*np.sqrt(sq[i])
                # Points around circle
                circle = self._circ_points(vec_a[i]*rad, vec_b[i]*rad, p[i])
                # If any points on the circle defining the intersection
                # Are within the pore radius distance from the pore center
                # The intersection is inside the pore - Not exact as pores
                # Are irregular shapes
                c2p = np.linalg.norm(circle-pore_center, axis=1)
                if np.any(c2p < pore_rad):
                    intersection[i] = True

        return intersection

    def setup_coop_filling(self, inv_points=None, capillary_model='purcell',
                           pores=None, radius=None):
        r"""
        Evaluate the cooperative pore filling condition that the combined
        filling angle in next neighbor throats cannot exceed the geometric
        angle between their throat planes.
        This is used when the invading fluid has access to multiple throats
        connected to a pore
        """
        net = self.project.network
        phases_physics = self.project.find_physics(phase=self._phase)
        from openpnm.models import physics as pm
        self._coop_fill = True
        if inv_points is None:
            inv_points = np.arange(0, 30100, 500)
        if pores is None:
            pores = range(self.Np)
        # Throat-Throat cooperative filling pressure
        self.tt_Pc = np.ndarray([self.Nt, self.Nt], dtype=float)
        self.tt_Pc.fill(sp.nan)
        # Dictionary with keys of capillary pressure
        self.coop_data = {}
        start = time.time()
        tfill_angle = 'throat.filling_angle'
        tmen_rad = 'throat.meniscus_radius'
        tmen_cen = 'throat.meniscus_center'
        try:
            # The following properties will all be there for Voronoi
            p_centroids = net['pore.centroid']
            t_centroids = net['throat.centroid']
            p_diam = net['pore.indiameter']
            t_norms = net['throat.normal']
        except:
            # Chances are this isn't Voronoi so calculate or replace all
            p_centroids = net['pore.coords']
            temp = net['pore.coords'][net['throat.conns']]
            t_centroids = np.mean(temp, axis=1)
            p_diam = net['pore.diameter']
            t_norms = net['throat.normal']

        if capillary_model == 'purcell':
            angle_model = pm.capillary_pressure.purcell_filling_angle
            radius_model = pm.capillary_pressure.purcell_meniscus_radius
            center_model = pm.capillary_pressure.purcell_meniscus_center
        elif capillary_model == 'sinusoidal':
            angle_model = pm.capillary_pressure.sinusoidal
            radius_model = pm.capillary_pressure.sinusoidal
            center_model = pm.capillary_pressure.sinusoidal
        else:
            logger.exception('capillary model '+capillary_model+' not valid')
        for phys in phases_physics:
            for Pc in inv_points:
                # Dictionary with keys of pore id
                pore_data = {}
                regen_mode='normal'
                if capillary_model == 'purcell':
                    phys.add_model(propname=tfill_angle,
                                   model=angle_model,
                                   r_toroid=radius,
                                   Pc=Pc,
                                   regen_mode=regen_mode)
                    phys.add_model(propname=tmen_rad,
                                   model=radius_model,
                                   r_toroid=radius,
                                   filling_angle=tfill_angle,
                                   regen_mode=regen_mode)
                    phys.add_model(propname=tmen_cen,
                                   model=center_model,
                                   r_toroid=radius,
                                   filling_angle=tfill_angle,
                                   regen_mode=regen_mode)
                elif capillary_model == 'sinusoidal':
                    phys.add_model(propname='throat.men_data',
                                   model=angle_model,
                                   mode='men',
                                   target_Pc=Pc,
                                   regen_mode=regen_mode)
                    phys[tfill_angle] = phys['throat.men_data.alpha']
                    phys[tmen_rad] = phys['throat.men_data.rad']
                    phys[tmen_cen] = phys['throat.men_data.cen']

                for pore in pores:
                    # Dictionary with keys of throat id
                    throat_data = {}
                    p_cen = p_centroids[pore]
                    p_rad = p_diam[pore]/2
                    throats = net.find_neighbor_throats(pores=pore, flatten=True)
                    throat_centres = t_centroids[throats]
                    throat_normals = t_norms[throats]
                    unit = np.linalg.norm(throat_normals, axis=1)
                    throat_normals /= np.vstack((unit, unit, unit)).T
                    v = p_cen - throat_centres
                    sign = np.sign(np.sum(v*throat_normals, axis=1))
                    cen = phys[tmen_cen][throats]
                    c3 = np.vstack((cen*sign, cen*sign, cen*sign)).T
                    men_cen = throat_centres + c3*throat_normals
                    pairs = []
                    for i, T in enumerate(throats):
                        men_data = {}
                        men_data['cen'] = men_cen[i]
                        men_data['rad'] = phys[tmen_rad][T]
                        men_data['offset'] = phys[tmen_cen][T]
                        men_data['alpha'] = phys[tfill_angle][T]
                        throat_data[T] = men_data
                    for ni in range(len(throats)):
                        for nj in range(len(throats))[ni+1:]:
                            pairs.append([ni, nj])
                    pairs = np.asarray(pairs)
                    if len(pairs) > 0:
                        c1 = men_cen[pairs[:, 0]]
                        c2 = men_cen[pairs[:, 1]]
                        c2c = (c1-c2)
                        dist = np.linalg.norm(c2c, axis=1)
                        t1 = throats[pairs[:, 0]]
                        t2 = throats[pairs[:, 1]]
                        r1 = phys[tmen_rad][t1]
                        r2 = phys[tmen_rad][t2]
                        # nans may exist if pressure is outside the range
                        # set these to zero to be ignored by next step without
                        # causing RuntimeWarning
                        r1[sp.isnan(r1)] = 0
                        r2[sp.isnan(r2)] = 0
                        check_pos = np.logical_and(r1 > 0, r2 > 0)
                        # simple initial distance check on sphere rads
                        check_rads = (r1+r2) >= dist
                        # check whether the filling angle is ok at this Pc
                        check_alpha_t1 = ~sp.isnan(phys[tfill_angle][t1])
                        check_alpha_t2 = ~sp.isnan(phys[tfill_angle][t2])
                        check_alpha = check_alpha_t1*check_alpha_t2
                        # check whether this throat pair already has a coop value
                        check_nans = sp.isnan(self.tt_Pc[t1, t2])
                        mask = check_pos*check_alpha*check_nans*check_rads
                        # if all checks pass
                        if np.any(mask):
                            # Check if intersecting circle lies within pore
                            inter = self._check_intersection(c1=c1[mask],
                                                             c2=c2[mask],
                                                             r1=r1[mask],
                                                             r2=r2[mask],
                                                             pore_center=p_cen,
                                                             pore_rad=p_rad)
                            if np.any(inter):
                                self.tt_Pc[t1[mask][inter], t2[mask][inter]] = Pc
                                self.tt_Pc[t2[mask][inter], t1[mask][inter]] = Pc
                            # intersection
                        # pairs of mensici
                    # pore
                    pore_data[pore] = throat_data
                # Pc
                self.coop_data[Pc] = pore_data
            # Physics
        logger.info("Coop filling finished in " +
                    str(np.around(time.time()-start, 2)) + " s")

    def _check_coop(self, pore, queue):
        r"""
        Method run in loop after every pore invasion. All connecting throats
        are now given access to the invading phase. Two throats with access to
        the invading phase can cooperatively fill any pores that they are both
        connected to, common pores.
        The invasion of the throats connected to the common pore is handled
        elsewhere.
        """
        net = self.project.network
        for throat in net.find_neighbor_throats(pores=pore):
            # A pore has just been invaded, all it's throats now have
            # An interface residing inside them
            if self['throat.inv_seq'][throat] == -1:
                # If the throat is not the invading throat that gave access
                # to this pore, get the pores that this throat connects with
                a = set(net['throat.conns'][throat])
                # Get a list of pre-calculated coop filling pressures for all
                # Throats this throat can coop fill with
                ts_Pc = self.tt_Pc[throat]
                ts = np.argwhere(~np.isnan(ts_Pc))
                # If there are any potential coop filling throats
                if len(ts) > 0:
                    # For each throat find the common pore and the uncommon
                    # pores
                    for t in ts.flatten():
                        # Find common pore (cP) and uncommon pores (uPs)
                        b = set(net['throat.conns'][t])
                        cP = list(a.intersection(b))
                        uPs = list(a.symmetric_difference(b))
                        # If the common pore is not invaded but the others are
                        # The potential coop filling event can now happen
                        # Add the coop pressure to the queue
                        if ((np.all(self['pore.inv_seq'][uPs] > -1)) and
                           (self['pore.inv_seq'][cP] == -1)):
                            # Coop pore filling fills the common pore
                            # The throats that gave access are not invaded now
                            # However, isolated throats between invaded pores
                            # Are taken care of elsewhere...
                            # This could result in cyclycally calling the
                            # functio.
                            # If all elements are added to the queue and this
                            # function does not have access to the invasion seq
                            hq.heappush(queue, [ts_Pc[t], list(cP), 'pore', 1])

    def _invade_isolated_Ts(self):
        r"""
        Throats that are uninvaded connected to pores that are both invaded
        should be invaded too.
        """
        net = self.project.network
        Ts = net['throat.conns'].copy()
        invaded_Ps = self['pore.inv_seq'] > -1
        uninvaded_Ts = self['throat.inv_seq'] == -1
        isolated_Ts = np.logical_and(invaded_Ps[Ts[:, 0]],
                                     invaded_Ps[Ts[:, 1]])
        isolated_Ts = np.logical_and(isolated_Ts, uninvaded_Ts)
        inv_Pc = self['pore.inv_Pc']
        if np.any(isolated_Ts):
            Pc = np.max(inv_Pc[Ts], axis=1)
            seq = np.max(inv_Pc[Ts], axis=1)
            max_array = Ts[:, 0]
            second_higher = inv_Pc[Ts][:, 1] > inv_Pc[Ts][:, 0]
            max_array[second_higher] = Ts[:, 1][second_higher]
            action = self['pore.action'][max_array]
            self['throat.inv_Pc'][isolated_Ts] = Pc[isolated_Ts]
            self['throat.inv_seq'][isolated_Ts] = seq[isolated_Ts]
            self['throat.action'][isolated_Ts] = action[isolated_Ts]
