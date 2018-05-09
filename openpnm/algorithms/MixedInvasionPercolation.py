# -*- coding: utf-8 -*-
"""
===============================================================================
InvasionPercolationBasic: Simple IP
===============================================================================

"""
import heapq as hq
import scipy as sp
import numpy as np
from openpnm.algorithms import GenericPercolation
import logging
import matplotlib.pyplot as plt
logger = logging.getLogger(__name__)


class MixedInvasionPercolation(GenericPercolation):
    r"""
    An implemetation of invasion percolation which can invade bonds, sites or a
    mixture of both. Inlets can be treated as individual injection points that
    share a common pressure or have their own and progess independently.
    Inlets can also be single pores or clusters.

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
        self.settings.update({'pore_volume': 'pore.volume',
                              'throat_volume': 'throat.volume',
                              'pore_entry_pressure': 'pore.capillary_pressure',
                              'throat_entry_pressure': 'throat.capillary_pressure',
                              'mode': 'mixed',
                              'partial_saturation': False,
                              'snap_off': False})

    def setup(self, phase, def_phase):
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
        self['throat.entry_pressure'] = phase[self.settings['throat_entry_pressure']]
        self['pore.entry_pressure'] = phase[self.settings['pore_entry_pressure']]
        self.reset()

    def reset(self):
        r"""
        Resets the various data arrays on the object back to their original
        state. This is useful for repeating a simulation at different inlet
        conditions, or invasion points for instance.
        """
        self['pore.invasion_pressure'] = np.inf
        self['throat.invasion_pressure'] = np.inf
        self['pore.invasion_sequence'] = -1
        self['throat.invasion_sequence'] = -1
        self['pore.trapped'] = np.inf
        self['throat.trapped'] = np.inf
        self['pore.inlets'] = False
        self['pore.outlets'] = False
        self['pore.residual'] = False
        self['throat.residual'] = False
        for elem in ['pore', 'throat']:
            for prop in ['invasion_saturation',
                         'occupancy',
                         'action']:
                try:
                    del self[elem+'.'+prop]
                except:
                    pass
            for prop in ['cluster',
                         'action']:
                try:
                    self[elem+'.'+prop] = -1
                except:
                    pass
        # Masks for tracking pores and throats at the interface
        # Basically a quick way of getting to all the elements in the queues
        self._interface_Ts = np.zeros(self.Nt, dtype=bool)
        self._interface_Ps = np.zeros(self.Np, dtype=bool)

    def set_inlets(self, pores=None, clusters=None, starting_sequence=-1):
        r"""

        Parameters
        ----------
        pores : array_like
            The list of inlet pores from which the Phase can enter the Network
        
        clusters : list of lists - can be just one list but each list defines
            a cluster of pores that share a common invasion pressure.
        
        starting_sequence: int default: -1
            The invasion sequence to set the inlets at, -1 is helpful for
            thresholding from a freshly started simulation. Other numbers may
            be used for starting from a partially saturated state which may or
            may not have been determined from running this algorithm multiple
            times.
            
        Like Basic Invasion Percolation a queue of 
        """
        if pores is not None:
            logger.info("Setting inlet pores at shared pressure")
            clusters = []
            clusters.append(pores)
        elif clusters is not None:
            logger.info("Setting inlet clusters at individual pressures")
        else:
            logger.error("Either 'inlets' or 'clusters' must be passed to" +
                         " setup method")

        self.queue = {}
        inlet_inv_seq = starting_sequence

        for i, cluster in enumerate(clusters):
            self.queue[i] = []
            # Perform initial analysis on input pores
            self['pore.invasion_sequence'][cluster] = inlet_inv_seq
            self['pore.cluster'][cluster] = i
            self['pore.invasion_pressure'][cluster] = -np.inf
            if np.size(cluster) > 1:
                for elem_id in cluster:
                    self._add_ts2q(elem_id, self.queue[i], action=0)
            elif np.size(cluster) == 1:
                self._add_ts2q(cluster, self.queue[i], action=0)
            else:
                logger.warning("Some inlet clusters have no pores")
        if self.settings['snap_off']:
            self._apply_snap_off()

        if self.settings['partial_saturation']:
            self.apply_partial_sat()
        else:
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
        Ts = Ts[self['throat.invasion_sequence'][Ts] <= 0]
        if len(Ts) > 0:
            self._interface_Ts[Ts] = True
            for T in Ts:
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
        Ps = Ps[self['pore.invasion_sequence'][Ps] <= 0]
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

    def run(self, max_pressure=None):
        r"""
        Perform the algorithm

        Parameters
        ----------
        n_steps : int
            The number of throats to invaded during this step

        """
        if 'throat.entry_pressure' not in self.keys():
            logger.error("Setup method must be run first")

        if max_pressure is None:
            self.max_pressure = sp.inf
        else:
            self.max_pressure = max_pressure
        if len(self.queue.items()) == 0:
            logger.warn('queue is empty, this network is fully invaded')
            return

        self.max_p_reached = [False]*len(self.queue.items())
        self.count = -1
        self.invasion_running = [True]*len(self.queue.items())
        self.high_Pc = np.ones(len(self.queue.items()))*-np.inf
        outlets = self['pore.outlets']
        terminate_clusters = np.sum(outlets) > 0
        while np.any(self.invasion_running) and not np.all(self.max_p_reached):
            # Loop over clusters
            for c_num in self.queue.keys():
                if self.invasion_running[c_num]:
                    self._invade_cluster(c_num)
                    queue = self.queue[c_num]
                    if len(queue) == 0 or self.max_p_reached[c_num]:
                        # If the cluster contains no more entries invasion has
                        # finished
                        self.invasion_running[c_num] = False
            self._invade_isolated_Ts()
            if terminate_clusters:
                # terminated clusters
                tcs = np.unique(self['pore.cluster'][outlets]).astype(int)
                tcs = tcs[tcs >= 0]
                if len(tcs) > 0:
                    for tc in tcs:
                        if self.invasion_running[tc] is True:
                            self.invasion_running[tc] = False
                            logger.info("Cluster " + str(tc) + " reached " +
                                        " outlet at sequence " + str(self.count))

    def _invade_cluster(self, c_num):
        queue = self.queue[c_num]
        pressure, elem_id, elem_type, action = hq.heappop(queue)
        if elem_type == 'pore':
            self._interface_Ps[elem_id] = False
        else:
            self._interface_Ts[elem_id] = False
        if pressure > self.max_pressure:
            self.max_p_reached[c_num] = True
        else:
            elem_cluster = self[elem_type+'.cluster'][elem_id]
            elem_cluster = elem_cluster.astype(int)
            # Cluster is the uninvaded cluster
            if elem_cluster == -1:
                self.count += 1
                # Record highest Pc cluster has reached
                if self.high_Pc[c_num] < pressure:
                    self.high_Pc[c_num] = pressure
                # The newly invaded element is available for
                # invasion
                self[elem_type+'.invasion_sequence'][elem_id] = self.count
                self[elem_type+'.cluster'][elem_id] = c_num
                self[elem_type+'.invasion_pressure'][elem_id] = self.high_Pc[c_num]
                self[elem_type+'.action'][elem_id] = action
                if elem_type == 'throat':
                    self._add_ps2q(elem_id, queue, action=0)
                elif elem_type == 'pore':
                    self._add_ts2q(elem_id, queue, action=0)
#                    if self._coop_fill:
#                        self._check_coop(elem_id, queue)
            # Element is part of existing cluster
            elif (elem_cluster != c_num and
                  self.invasion_running[elem_cluster]):
                # The newly invaded element is part of an invading
                # cluster. Merge the clusters using the existing
                # cluster number
                logger.info("Merging cluster "+str(c_num) +
                            " into cluster "+str(elem_cluster) +
                            " at sequence "+str(self.count))

            elif (elem_cluster != c_num and not
                  self.invasion_running[elem_cluster]):
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
        self._phase['throat.invasion_sequence'] = sp.nan
        self._phase['pore.invasion_sequence'] = sp.nan
        self._phase['throat.invasion_sequence'][throats] = \
            self['throat.invasion_sequence'][throats]
        self._phase['pore.invasion_sequence'][pores] = \
            self['pore.invasion_sequence'][pores]
        self._phase['throat.cluster'] = self['throat.cluster']
        self._phase['pore.cluster'] = self['pore.cluster']
        if Pc is not None:
            if 'pore.invasion_pressure' in self.keys():
                self['throat.occupancy'] = self['throat.invasion_pressure'] <= Pc
                self['pore.occupancy'] = self['pore.invasion_pressure'] <= Pc
                self._phase['throat.occupancy'] = self['throat.occupancy']
                self._phase['pore.occupancy'] = self['pore.occupancy']
                self._def['throat.occupancy'] = ~self['throat.occupancy']
                self._def['pore.occupancy'] = ~self['pore.occupancy']
            else:
                logger.warning("Occupancy not updated, please run " +
                               "extract_drainage() method to populate" +
                               " the invasion_pressure data")
        if "pore.invasion_pressure" in self.props():
            self._phase['pore.invasion_pressure'] = \
                self['pore.invasion_pressure']
            self._phase['throat.invasion_pressure'] = \
                self['throat.invasion_pressure']
        if "pore.invasion_saturation" in self.props():
            self._phase['pore.invasion_saturation'] = \
                self['pore.invasion_saturation']
            self._phase['throat.invasion_saturation'] = \
                self['throat.invasion_saturation']
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
        aa = self['throat.invasion_sequence']
        bb = sp.argsort(self['throat.invasion_sequence'])
        P12_inv = self['pore.invasion_sequence'][P12]
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
        if "pore.invasion_pressure" not in self.props():
            logger.error("Cannot plot drainage curve. Please run " +
                         " algorithm first")
        if inv_points is None:
            mask = ~sp.isnan(self['throat.invasion_pressure'])
            ok_Pc = self['throat.invasion_pressure'][mask]
            inv_points = np.unique(ok_Pc)
        sat_p = np.zeros(len(inv_points))
        sat_t = np.zeros(len(inv_points))
        inv_p = self['pore.invasion_pressure']
        inv_t = self['throat.invasion_pressure']
        # Handle trapped pores and throats
        if np.sum(self['pore.invasion_sequence'] == -1) > 0:
            inv_p[self['pore.invasion_sequence'] == -1] = 2*inv_points.max()
        if np.sum(self['throat.invasion_sequence'] == -1) > 0:
            inv_t[self['throat.invasion_sequence'] == -1] = 2*inv_points.max()
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
            Swp = np.ones(len(self['pore.invasion_pressure']))
            mask = self['pore.invasion_pressure'] < Pc
            Swp[mask] = Swp_init*(self['pore.invasion_pressure'][mask]/Pc)**eta
            Swp[np.isnan(Swp)] = 1.0
        except:
            Swp = np.ones(net.Np)
        Snwp = 1-Swp
        if wetting_phase:
            return Swp
        else:
            return Snwp

    def apply_trapping(self, partial=False):
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
        partial : Boolean indicating whether partially filled network

        Returns
        -------
        Creates a throat array called 'pore.clusters' in the Algorithm
        dictionary. Any positive number is a trapped cluster

        Also creates 2 boolean arrays Np and Nt long called '<element>.trapped'
        """
        net = self.project.network
        inlets = self['pore.inlets']
        outlets = self['pore.outlets']
        if np.sum(outlets) == 0:
            raise Exception('Outlets must be set using the set_outlets method' +
                            ' before applying trapping')
        if partial:
            # Set occupancy
            invaded_ps = self['pore.invasion_sequence'] > -1
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
        inv_seq = np.vstack((self['pore.invasion_sequence'].astype(int),
                             np.arange(0, net.Np, dtype=int))).T
        # Reverse sort list
        inv_seq = inv_seq[inv_seq[:, 0].argsort()][::-1]
        next_cluster_num = np.max(clusters)+1
        # For all the steps after the inlets are set up to break-through
        # Reverse the sequence and assess the neighbors cluster state
        stopped_clusters = np.zeros(net.Np, dtype=bool)
        all_neighbors = net.find_neighbor_pores(net.pores(), flatten=False)
        for un_seq, pore in inv_seq:
            if ~outlets[pore] and un_seq > -1:  # Don't include outlets
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
            logger.info("Number of trapped clusters " + str(num_trap))
            self['pore.trapped'] = clusters > -1
            num_tPs = np.sum(self['pore.trapped'])
            logger.info("Number of trapped pores: " + str(num_tPs))
            self['pore.invasion_sequence'][self['pore.trapped']] = -1
            self['throat.trapped'] = np.zeros([net.Nt], dtype=bool)
            for c in np.unique(clusters[clusters >= 0]):
                c_ts = net.find_neighbor_throats(clusters == c,
                                                 mode='intersection')
                self['throat.trapped'][c_ts] = True
            num_tTs = np.sum(self['throat.trapped'])
            logger.info("Number of trapped throats: " + str(num_tTs))
            self['throat.invasion_sequence'][self['throat.trapped']] = -1
            # Assumes invasion has run to the end
            self._phase['pore.occupancy'] = ~self['pore.trapped']
            self._def['pore.occupancy'] = self['pore.trapped']
            self._phase['throat.occupancy'] = ~self['throat.trapped']
            self._def['throat.occupancy'] = self['throat.trapped']
        else:
            logger.info("No trapped clusters found")

    def _apply_snap_off(self, snap_off='throat.snap_off', queue=None):
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
            logger.info("Applying partial saturation to " +
                        str(np.sum(occ_Ps)) + " pores")
            self['pore.invasion_sequence'][occ_Ps] = 0
            for P in net.pores()[occ_Ps]:
                self._add_ts2q(P, queue, action=0)
                self['pore.cluster'][P] = invading_cluster
                self['pore.invasion_pressure'][P] = low_val
        if np.sum(occ_Ts) > 0:
            logger.info("Applying partial saturation to " +
                        str(np.sum(occ_Ts)) + " throats")
        self['throat.invasion_sequence'][occ_Ts] = 0
        for T in net.throats()[occ_Ts]:
            self['throat.cluster'][T] = invading_cluster
            self['throat.invasion_pressure'][T] = low_val

    def _invade_isolated_Ts(self):
        r"""
        Throats that are uninvaded connected to pores that are both invaded
        should be invaded too.
        """
        net = self.project.network
        Ts = net['throat.conns'].copy()
        invaded_Ps = self['pore.invasion_sequence'] > -1
        uninvaded_Ts = self['throat.invasion_sequence'] == -1
        isolated_Ts = np.logical_and(invaded_Ps[Ts[:, 0]],
                                     invaded_Ps[Ts[:, 1]])
        isolated_Ts = np.logical_and(isolated_Ts, uninvaded_Ts)
        inv_Pc = self['pore.invasion_pressure']
        inv_seq = self['pore.invasion_sequence']
        if np.any(isolated_Ts):
            max_array = Ts[:, 0]
            second_higher = inv_seq[Ts][:, 1] > inv_seq[Ts][:, 0]
            max_array[second_higher] = Ts[:, 1][second_higher]
            mPc = inv_Pc[max_array]
            mSeq = inv_seq[max_array]
            mAct = self['pore.action'][max_array]
            mClu = self['pore.cluster'][max_array]
            self['throat.invasion_pressure'][isolated_Ts] = mPc[isolated_Ts]
            self['throat.invasion_sequence'][isolated_Ts] = mSeq[isolated_Ts]
            self['throat.action'][isolated_Ts] = mAct[isolated_Ts]
            self['throat.cluster'][isolated_Ts] = mClu[isolated_Ts]
