"""
===============================================================================
MixedInvasionPercolation: IP allowing pores and throats to invade separately
===============================================================================

"""
import logging
import heapq as hq
import numpy as np
from collections import namedtuple
from openpnm.algorithms import GenericAlgorithm
from openpnm.topotools import find_clusters, site_percolation

logger = logging.getLogger(__name__)


class MixedInvasionPercolation(GenericAlgorithm):
    r"""
    An implemetation of invasion percolation which can invade bonds,
    sites or a mixture of both. Inlets can be treated as individual
    injection points that share a common pressure or have their own and
    progess independently. Inlets can also be single pores or clusters.

    Parameters
    ----------
    network : GenericNetwork
        The network upon which the invasion should occur.

    Notes
    -----
    n/a

    """

    def __init__(self, settings={}, **kwargs):
        def_set = {
            "pore_entry_pressure": "pore.entry_pressure",
            "throat_entry_pressure": "throat.entry_pressure",
            "snap_off": "",
            "invade_isolated_Ts": False,
            "late_pore_filling": "",
            "late_throat_filling": "",
            "gui": {
                "setup": {
                    "pore_entry_pressure": "",
                    "throat_entry_pressure": "",
                    "snap_off": "",
                    "invade_isolated_Ts": "",
                },
                "set_inlets": {"pores": None, "clusters": None},
                "set_outlets": {"pores": None, "overwrite": False},
                "apply_flow": {"flowrate": None},
                "apply_trapping": {"partial": False},
                "set_residual": {"pores": None, "overwrite": False},
            },
        }
        super().__init__(**kwargs)
        self.settings.update(def_set)
        self.settings.update(settings)

    def setup(
        self,
        phase=None,
        pore_entry_pressure="pore.entry_pressure",
        throat_entry_pressure="throat.entry_pressure",
        snap_off="",
        invade_isolated_Ts=False,
        late_pore_filling="",
        late_throat_filling="",
    ):
        r"""
        Used to specify necessary arguments to the simulation.  This method is
        useful for resetting the algorithm or applying more explicit control.

        Parameters
        ----------
        phase : OpenPNM Phase object
            The Phase object containing the physical properties of the invading
            fluid.

        pore_entry_pressure : string
            The dictionary key on the Phase object where the pore entry
            pressure values are stored.  The default is
            'pore.entry_pressure'.

        throat_entry_pressure : string
            The dictionary key on the Phase object where the throat entry
            pressure values are stored.  The default is
            'throat.entry_pressure'.

        snap_off : string
            The dictionary key on the Phase object where the throat snap-off
            pressure values are stored.

        invade_isolated_Ts : boolean
            If True, isolated throats are invaded at the higher invasion
            pressure of their connected pores.

        late_pore_filling : string
            The name of the model used to determine late pore filling as
            a function of applied pressure.

        late_throat_filling : string
            The name of the model used to determine late throat filling as
            a function of applied pressure.

        """
        if phase:
            self.settings["phase"] = phase.name
        if throat_entry_pressure:
            self.settings["throat_entry_pressure"] = throat_entry_pressure
            phase = self.project.find_phase(self)
        self["throat.entry_pressure"] = phase[self.settings["throat_entry_pressure"]]
        if len(np.shape(self["throat.entry_pressure"])) > 1:
            self._bidirectional = True
        else:
            self._bidirectional = False
        if pore_entry_pressure:
            self.settings["pore_entry_pressure"] = pore_entry_pressure
            phase = self.project.find_phase(self)
        self["pore.entry_pressure"] = phase[self.settings["pore_entry_pressure"]]
        if snap_off:
            self.settings["snap_off"] = snap_off
        if invade_isolated_Ts:
            self.settings["invade_isolated_Ts"] = invade_isolated_Ts
        if late_pore_filling:
            self.settings["late_pore_filling"] = late_pore_filling
        if late_throat_filling:
            self.settings["late_throat_filling"] = late_throat_filling
        self.reset()

    def reset(self):
        r"""
        Resets the various data arrays on the object back to their original
        state. This is useful for repeating a simulation at different inlet
        conditions, or invasion points for instance.
        """
        self["pore.invasion_pressure"] = np.inf
        self["throat.invasion_pressure"] = np.inf
        self["pore.invasion_sequence"] = -1
        self["throat.invasion_sequence"] = -1
        self["pore.invasion_saturation"] = -1
        self["throat.invasion_saturation"] = -1
        self["pore.cluster"] = -1
        self["throat.cluster"] = -1
        self["pore.trapped"] = np.inf
        self["throat.trapped"] = np.inf
        self["pore.inlets"] = False
        self["pore.outlets"] = False
        self["pore.residual"] = False
        self["throat.residual"] = False
        for elem in ["pore", "throat"]:
            for prop in ["occupancy"]:
                try:
                    del self[elem + "." + prop]
                except KeyError:
                    pass

        # Masks for tracking pores and throats at the interface
        # Basically a quick way of getting to all the elements in the queues
        self._interface_Ts = np.zeros(self.Nt, dtype=bool)
        self._interface_Ps = np.zeros(self.Np, dtype=bool)
        if hasattr(self, "invasion_running"):
            del self.invasion_running

    def set_inlets(self, pores=None, clusters=None):
        r"""

        Parameters
        ----------
        pores : array_like
            The list of inlet pores from which the Phase can enter the Network

        clusters : list of lists - can be just one list but each list defines
            a cluster of pores that share a common invasion pressure.

        Like Basic Invasion Percolation a queue of
        """
        if pores is not None:
            logger.info("Setting inlet pores at shared pressure")
            clusters = []
            clusters.append(pores)
        elif clusters is not None:
            logger.info("Setting inlet clusters at individual pressures")
        else:
            logger.error("Either 'inlets' or 'clusters' must be passed to" + " setup method")
        self.queue = []
        check_coop = self.settings["cooperative_pore_filling"] and hasattr(self, "tt_Pc")
        for i, cluster in enumerate(clusters):
            self.queue.append([])
            # Perform initial analysis on input pores
            self["pore.invasion_sequence"][cluster] = 0
            self["pore.cluster"][cluster] = i
            self["pore.invasion_pressure"][cluster] = -np.inf
            if np.size(cluster) > 1:
                for elem_id in cluster:
                    self._add_ts2q(elem_id, self.queue[i])
                    if check_coop:
                        self._check_coop(elem_id, self.queue[i])
            elif np.size(cluster) == 1:
                self._add_ts2q(cluster, self.queue[i])
                if check_coop:
                    self._check_coop(cluster, self.queue[i])
            else:
                logger.warning("Some inlet clusters have no pores")
        if self.settings["snap_off"]:
            self._apply_snap_off()

    def set_outlets(self, pores=[], overwrite=False):
        r"""
        Set the locations through which defender exits the network.
        This is only necessary if 'trapping' was set to True when ``setup``
        was called.

        Parameters
        ----------
        pores : array_like
            Locations where the defender can exit the network.  Any defender
            that does not have access to these sites will be trapped.

        overwrite : boolean
            If ``True`` then all existing outlet locations will be removed and
            then the supplied locations will be added.  If ``False`` (default),
            then supplied locations are added to any already existing outlet
            locations.

        """
        if self.settings["trapping"] is False:
            logger.warning(
                "Setting outlets is meaningless unless trapping "
                + "was set to True during setup"
            )
        Ps = self._parse_indices(pores)
        if np.sum(self["pore.inlets"][Ps]) > 0:
            raise Exception("Some outlets are already defined as inlets")
        if overwrite:
            self["pore.outlets"] = False
        self["pore.outlets"][Ps] = True

    def _add_ts2q(self, pore, queue):
        """
        Helper method to add throats to the cluster queue
        """
        net = self.project.network
        elem_type = "throat"
        # Find throats connected to newly invaded pore
        Ts = net.find_neighbor_throats(pores=pore)
        # Remove already invaded throats from Ts
        Ts = Ts[self["throat.invasion_sequence"][Ts] <= 0]
        tcp = self["throat.entry_pressure"]
        if len(Ts) > 0:
            self._interface_Ts[Ts] = True
            for T in Ts:
                data = []
                # Pc
                if self._bidirectional:
                    # Get index of pore being invaded next and apply correct
                    # entry pressure
                    pmap = net["throat.conns"][T] != pore
                    pind = list(pmap).index(True)
                    data.append(tcp[T][pind])
                else:
                    data.append(tcp[T])
                # Element Index
                data.append(T)
                # Element Type (Pore of Throat)
                data.append(elem_type)
                hq.heappush(queue, data)

    def _add_ps2q(self, throat, queue):
        """
        Helper method to add pores to the cluster queue
        """
        net = self.project.network
        elem_type = "pore"
        # Find pores connected to newly invaded throat
        Ps = net["throat.conns"][throat]
        # Remove already invaded pores from Ps
        Ps = Ps[self["pore.invasion_sequence"][Ps] <= 0]
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
                hq.heappush(queue, data)

    def run(self, max_pressure=None):
        r"""
        Perform the algorithm

        Parameters
        ----------
        max_pressure : float
            The maximum pressure applied to the invading cluster. Any pores and
            throats with entry pressure above this value will not be invaded.

        """
        if "throat.entry_pressure" not in self.keys():
            logger.error("Setup method must be run first")

        if max_pressure is None:
            self.max_pressure = np.inf
        else:
            self.max_pressure = max_pressure
        if len(self.queue) == 0:
            logger.warn("queue is empty, this network is fully invaded")
            return
        # track whether each cluster has reached the maximum pressure
        self.max_p_reached = [False] * len(self.queue)
        # starting invasion sequence
        self.count = 0
        # highest pressure reached so far - used for porosimetry curve
        self.high_Pc = np.ones(len(self.queue)) * -np.inf
        outlets = self["pore.outlets"]
        terminate_clusters = np.sum(outlets) > 0
        if not hasattr(self, "invasion_running"):
            self.invasion_running = [True] * len(self.queue)
        else:
            # created by set_residual
            pass
        while np.any(self.invasion_running) and not np.all(self.max_p_reached):
            # Loop over clusters
            for c_num in np.argwhere(self.invasion_running).flatten():
                self._invade_cluster(c_num)
                queue = self.queue[c_num]
                if len(queue) == 0 or self.max_p_reached[c_num]:
                    # If the cluster contains no more entries invasion has
                    # finished
                    self.invasion_running[c_num] = False
            if self.settings["invade_isolated_Ts"]:
                self._invade_isolated_Ts()
            if terminate_clusters:
                # terminated clusters
                tcs = np.unique(self["pore.cluster"][outlets]).astype(int)
                tcs = tcs[tcs >= 0]
                if len(tcs) > 0:
                    for tc in tcs:
                        if self.invasion_running[tc] is True:
                            self.invasion_running[tc] = False
                            logger.info(
                                "Cluster "
                                + str(tc)
                                + " reached "
                                + " outlet at sequence "
                                + str(self.count)
                            )

    def _invade_cluster(self, c_num):
        queue = self.queue[c_num]
        pressure, elem_id, elem_type = hq.heappop(queue)
        if elem_type == "pore":
            self._interface_Ps[elem_id] = False
        else:
            self._interface_Ts[elem_id] = False
        if pressure > self.max_pressure:
            self.max_p_reached[c_num] = True
        else:
            elem_cluster = self[elem_type + ".cluster"][elem_id]
            elem_cluster = elem_cluster.astype(int)
            # Cluster is the uninvaded cluster
            if elem_cluster == -1:
                self.count += 1
                # Record highest Pc cluster has reached
                if self.high_Pc[c_num] < pressure:
                    self.high_Pc[c_num] = pressure
                # The newly invaded element is available for
                # invasion
                self[elem_type + ".invasion_sequence"][elem_id] = self.count
                self[elem_type + ".cluster"][elem_id] = c_num
                self[elem_type + ".invasion_pressure"][elem_id] = self.high_Pc[c_num]
                if elem_type == "throat":
                    self._add_ps2q(elem_id, queue)
                elif elem_type == "pore":
                    self._add_ts2q(elem_id, queue)
                    if self.settings["cooperative_pore_filling"] and hasattr(
                        self, "tt_Pc"
                    ):
                        self._check_coop(elem_id, queue)
            # Element is part of existing cluster that is still invading
            elif elem_cluster != c_num and self.invasion_running[elem_cluster]:
                # The newly invaded element is part of an invading
                # cluster. Merge the clusters using the existing
                # cluster number)
                self._merge_cluster(c2keep=c_num, c2empty=elem_cluster)
                logger.info(
                    "Merging cluster "
                    + str(elem_cluster)
                    + " into cluster "
                    + str(c_num)
                    + " at sequence "
                    + str(self.count)
                )
            # Element is part of residual cluster - now invasion can start
            elif elem_cluster != c_num and len(self.queue[elem_cluster]) > 0:
                # The newly invaded element is part of an invading
                # cluster. Merge the clusters using the existing
                # cluster number)
                self._merge_cluster(c2keep=c_num, c2empty=elem_cluster)
                logger.info(
                    "Merging residual cluster "
                    + str(elem_cluster)
                    + " into cluster "
                    + str(c_num)
                    + " at sequence "
                    + str(self.count)
                )
            else:
                pass

    def _merge_cluster(self, c2keep, c2empty):
        r"""
        Little helper function to merger clusters but only add the uninvaded
        elements
        """
        while len(self.queue[c2empty]) > 0:
            temp = [_pc, _id, _type] = hq.heappop(self.queue[c2empty])
            if self[_type + ".invasion_sequence"][_id] == -1:
                hq.heappush(self.queue[c2keep], temp)
        self.invasion_running[c2empty] = False

    def results(self, Pc):
        r"""
        Places the results of the IP simulation into the Phase object.

        Parameters
        ----------
        Pc : float
            Capillary Pressure at which phase configuration was reached

        """
        if Pc is None:
            results = {
                "pore.invasion_sequence": self["pore.invasion_sequence"],
                "throat.invasion_sequence": self["throat.invasion_sequence"],
            }
        else:
            phase = self.project.find_phase(self)
            net = self.project.network
            inv_p = self["pore.invasion_pressure"].copy()
            inv_t = self["throat.invasion_pressure"].copy()
            # Handle trapped pores and throats by moving their pressure up
            # to be ignored
            if np.sum(self["pore.invasion_sequence"] == -1) > 0:
                inv_p[self["pore.invasion_sequence"] == -1] = Pc + 1
            if np.sum(self["throat.invasion_sequence"] == -1) > 0:
                inv_t[self["throat.invasion_sequence"] == -1] = Pc + 1
            p_inv = inv_p <= Pc
            t_inv = inv_t <= Pc

            if self.settings["late_pore_filling"]:
                # Set pressure on phase to current capillary pressure
                phase["pore.pressure"] = Pc
                # Regenerate corresponding physics model
                for phys in self.project.find_physics(phase=phase):
                    phys.regenerate_models(self.settings["late_pore_filling"])
                # Fetch partial filling fraction from phase object (0->1)
                frac = phase[self.settings["late_pore_filling"]]
                p_vol = net["pore.volume"] * frac
            else:
                p_vol = net["pore.volume"]
            if self.settings["late_throat_filling"]:
                # Set pressure on phase to current capillary pressure
                phase["throat.pressure"] = Pc
                # Regenerate corresponding physics model
                for phys in self.project.find_physics(phase=phase):
                    phys.regenerate_models(self.settings["late_throat_filling"])
                # Fetch partial filling fraction from phase object (0->1)
                frac = phase[self.settings["late_throat_filling"]]
                t_vol = net["throat.volume"] * frac
            else:
                t_vol = net["throat.volume"]
            results = {
                "pore.occupancy": p_inv * p_vol,
                "throat.occupancy": t_inv * t_vol,
            }
        return results

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
        P12 = net["throat.conns"]
        aa = self["throat.invasion_sequence"]
        bb = np.argsort(self["throat.invasion_sequence"])
        P12_inv = self["pore.invasion_sequence"][P12]
        # Find if the connected pores were invaded with or before each throat
        P1_inv = P12_inv[:, 0] <= aa
        P2_inv = P12_inv[:, 1] <= aa
        cc = np.column_stack((P1_inv, P2_inv))
        # List of Pores invaded with each throat
        dd = np.sum(cc, axis=1, dtype=bool)
        # Find volume of these pores
        P12_vol = np.sum(net["pore.volume"][P12] * cc, axis=1) * dd
        # Add invaded throat volume to pore volume (if invaded)
        T_vol = P12_vol + net["throat.volume"]
        # Cumulative sum on the sorted throats gives cumulated inject volume
        ee = np.cumsum(T_vol[bb] / flowrate)
        t = np.zeros((self.Nt,))
        t[bb] = ee  # Convert back to original order
        phase = self.project.find_phase(self)
        phase["throat.invasion_time"] = t

    def get_intrusion_data(self, inv_points=None):
        r"""
        Returns the intrusion data.
        """
        net = self.project.network

        if "pore.invasion_pressure" not in self.props():
            logger.error("You must run the algorithm first.")
            return None

        if inv_points is None:
            mask = ~np.isnan(self["throat.invasion_pressure"])
            ok_Pc = self["throat.invasion_pressure"][mask]
            inv_points = np.unique(ok_Pc)
        sat_p = np.zeros(len(inv_points))
        sat_t = np.zeros(len(inv_points))

        for i, Pc in enumerate(inv_points):
            res = self.results(Pc=Pc)
            sat_p[i] = np.sum(res["pore.occupancy"])
            sat_t[i] = np.sum(res["throat.occupancy"])

        pvol = np.sum(net["pore.volume"])
        tvol = np.sum(net["throat.volume"])
        tot_vol = pvol + tvol
        tot_sat = sat_p + sat_t

        # Normalize
        sat_p /= tot_vol
        sat_t /= tot_vol
        tot_sat /= tot_vol
        pc_curve = namedtuple("pc_curve", ("Pcap", "S_pore", "S_throat", "S_tot"))
        data = pc_curve(inv_points, sat_p, sat_t, tot_sat)

        return data

    def plot_intrusion_curve(self, ax=None, inv_points=None, num_markers=25):
        r"""
        Plot a simple drainage curve.
        """
        import matplotlib.pyplot as plt

        data = self.get_intrusion_data(inv_points)
        if data is None:
            raise Exception("You must run the algorithm first.")
        if ax is None:
            fig, ax = plt.subplots()
        markevery = max(data.Pcap.size // num_markers, 1)
        ax.plot(data.Pcap, data.S_pore, "r*-", label="pore", markevery=markevery)
        ax.plot(data.Pcap, data.S_throat, "b*-", label="throat", markevery=markevery)
        ax.plot(data.Pcap, data.S_tot, "g*-", label="total", markevery=markevery)
        ax.legend(bbox_to_anchor=(0, 1.02, 1, 0.102), loc=3, ncol=3, borderaxespad=0)
        ax.set_xlabel("Capillary pressure [Pa]")
        ax.set_ylabel("Saturation")
        ax.set_ybound(lower=0.0, upper=1.0)

    def apply_trapping(self, partial=False):
        """
        Apply trapping based on algorithm described by Y. Masson [1].

        Parameters
        ----------
        partial : boolean
            Indicating whether partially filled network

        Notes
        -----
        It is applied as a post-process and runs the percolation algorithm in
        reverse assessing the occupancy of pore neighbors. 3 situations can
        happen on invasion without trapping:

        * The number of defending clusters stays the same and clusters can shrink
        * A cluster of size one is suppressed
        * A cluster is split into multiple clusters

        In reverse the following situations can happen:

        * The number of defending clusters stays the same and clusters can grow
        * A cluster of size one is created
        * Mutliple clusters merge into one cluster

        With trapping the reversed rules are adjusted so that:

        * Only clusters that do not connect to a sink can grow and merge.
        * At the point that a neighbor connected to a sink is touched, the
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

        References
        ----------
        [1] Masson, Y., 2016. A fast two-step algorithm for invasion
        percolation with trapping. Computers & Geosciences, 90, pp.41-48

        Returns
        -------
        Creates a throat array called 'pore.clusters' in the Algorithm
        dictionary. Any positive number is a trapped cluster. Also creates 2
        boolean arrays Np and Nt long called '<element>.trapped'

        """
        net = self.project.network
        outlets = self["pore.outlets"]
        if np.sum(outlets) == 0:
            raise Exception(
                "Outlets must be set using the set_outlets method"
                + " before applying trapping"
            )
        if partial:
            # Set occupancy
            invaded_ps = self["pore.invasion_sequence"] > -1
            # Put defending phase into clusters
            clusters = find_clusters(network=net, pmask=~invaded_ps)
            # Identify clusters that are connected to an outlet and set to -2
            # -1 is the invaded fluid
            # -2 is the defender fluid able to escape
            # All others now trapped clusters which grow as invasion is
            # reversed
            out_clusters = np.unique(clusters[outlets])
            for c in out_clusters:
                if c >= 0:
                    clusters[clusters == c] = -2
        else:
            # Go from end
            clusters = np.ones(net.Np, dtype=int) * -1
            clusters[outlets] = -2

        # Turn into a list for indexing
        inv_seq = np.vstack(
            (
                self["pore.invasion_sequence"].astype(int),
                np.arange(0, net.Np, dtype=int),
            )
        ).T
        # Reverse sort list
        inv_seq = inv_seq[inv_seq[:, 0].argsort()][::-1]
        next_cluster_num = np.max(clusters) + 1
        # For all the steps after the inlets are set up to break-through
        # Reverse the sequence and assess the neighbors cluster state
        stopped_clusters = np.zeros(net.Np, dtype=bool)
        all_neighbors = net.find_neighbor_pores(
            net.pores(), flatten=False, include_input=True
        )
        for un_seq, pore in inv_seq:
            if ~outlets[pore] and un_seq > -1:  # Don't include outlets
                nc = clusters[all_neighbors[pore]]  # Neighboring clusters
                unique_ns = np.unique(nc[nc != -1])  # Unique Neighbors
                seq_pore = "S:" + str(un_seq) + " P:" + str(pore)
                if np.all(nc == -1):
                    # This is the start of a new trapped cluster
                    clusters[pore] = next_cluster_num
                    next_cluster_num += 1
                    msg = seq_pore + " C:1 new cluster number: " + str(clusters[pore])
                    logger.info(msg)
                elif len(unique_ns) == 1:
                    # Grow the only connected neighboring cluster
                    if not stopped_clusters[unique_ns[0]]:
                        clusters[pore] = unique_ns[0]
                        msg = (
                            seq_pore
                            + " C:2 joins cluster number: "
                            + str(clusters[pore])
                        )
                        logger.info(msg)
                    else:
                        clusters[pore] = -2
                elif -2 in unique_ns:
                    # We have reached a sink neighbor, stop growing cluster
                    msg = seq_pore + " C:3 joins sink cluster"
                    logger.info(msg)
                    clusters[pore] = -2
                    # Stop growth and merging
                    stopped_clusters[unique_ns[unique_ns > -1]] = True
                else:
                    # We might be able to do some merging
                    # Check if any stopped clusters are neighbors
                    if np.any(stopped_clusters[unique_ns]):
                        msg = seq_pore + " C:4 joins sink cluster"
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
                            msg = (
                                seq_pore
                                + " C:5 merge clusters: "
                                + str(c)
                                + " into "
                                + str(new_num)
                            )
                            logger.info(msg)

        # And now return clusters
        clusters[outlets] = -2
        num_trap = np.sum(np.unique(clusters) >= 0)
        if num_trap > 0:
            logger.info("Number of trapped clusters " + str(num_trap))
            self["pore.trapped"] = clusters > -1
            num_tPs = np.sum(self["pore.trapped"])
            logger.info("Number of trapped pores: " + str(num_tPs))
            self["pore.invasion_sequence"][self["pore.trapped"]] = -1
            self["throat.trapped"] = np.zeros([net.Nt], dtype=bool)
            for c in np.unique(clusters[clusters >= 0]):
                c_ts = net.find_neighbor_throats(clusters == c, mode="xnor")
                self["throat.trapped"][c_ts] = True
            num_tTs = np.sum(self["throat.trapped"])
            logger.info("Number of trapped throats: " + str(num_tTs))
            self["throat.invasion_sequence"][self["throat.trapped"]] = -1
            # Assumes invasion has run to the end
            phase = self.project.find_phase(self)
            phase["pore.occupancy"] = ~self["pore.trapped"]
            phase["throat.occupancy"] = ~self["throat.trapped"]
        else:
            logger.info("No trapped clusters found")

    def _apply_snap_off(self, queue=None):
        r"""
        Add all the throats to the queue with snap off pressure
        This is probably wrong!!!! Each one needs to start a new cluster.
        """
        net = self.project.network
        phase = self.project.find_phase(self)
        snap_off = self.settings["snap_off"]
        if queue is None:
            queue = self.queue[0]
        try:
            Pc_snap_off = phase[snap_off]
            logger.info("Adding snap off pressures to queue")
            for T in net.throats():
                if not np.isnan(Pc_snap_off[T]):
                    hq.heappush(queue, [Pc_snap_off[T], T, "throat"])
        except KeyError:
            logger.warning(
                "Phase " + phase.name + " doesn't have " + "property " + snap_off
            )

    def set_residual(self, pores=[], overwrite=False):
        r"""
        Method to start invasion in a network w. residual saturation.
        Called after inlets are set.

        Parameters
        ----------
        pores : array_like
            The pores locations that are to be filled with invader at the
            beginning of the simulation.

        overwrite : boolean
            If ``True`` then all existing inlet locations will be removed and
            then the supplied locations will be added.  If ``False``, then
            supplied locations are added to any already existing locations.

        Notes
        -----
        Currently works for pores only and treats inner throats, i.e.
        those that connect two pores in the cluster as invaded and outer ones
        as uninvaded. Uninvaded throats are added to a new residual cluster
        queue but do not start invading independently if not connected to an
        inlet.

        Step 1. Identify clusters in the phase occupancy.
        Step 2. Look for clusters that are connected or contain an inlet
        Step 3. For those that are merge into inlet cluster. May be connected
        to more than one - run should sort this out
        Step 4. For those that are isolated set the queue to not invading.
        Step 5. (in run) When isolated cluster is met my invading cluster it
        merges in and starts invading


        """
        Ps = self._parse_indices(pores)
        if overwrite:
            self["pore.residual"] = False
        self["pore.residual"][Ps] = True
        residual = self["pore.residual"]
        net = self.project.network
        conns = net["throat.conns"]
        rclusters = site_percolation(conns, residual).sites
        rcluster_ids = np.unique(rclusters[rclusters > -1])
        initial_num = len(self.queue) - 1
        for rcluster_id in rcluster_ids:
            rPs = rclusters == rcluster_id
            existing = np.unique(self["pore.cluster"][rPs])
            existing = existing[existing > -1]
            if len(existing) > 0:
                # There was at least one inlet cluster connected to this
                # residual cluster, pick the first one.
                cluster_num = existing[0]
            else:
                # Make a new cluster queue
                cluster_num = len(self.queue)
                self.queue.append([])
            queue = self.queue[cluster_num]
            # Set the residual pores and inner throats as part of cluster
            self["pore.cluster"][rPs] = cluster_num
            Ts = net.find_neighbor_throats(pores=rPs, flatten=True, mode="xnor")
            self["throat.cluster"][Ts] = cluster_num
            self["pore.invasion_sequence"][rPs] = 0
            self["throat.invasion_sequence"][Ts] = 0
            self["pore.invasion_pressure"][rPs] = -np.inf
            self["throat.invasion_pressure"][Ts] = -np.inf
            # Add all the outer throats to the queue
            Ts = net.find_neighbor_throats(pores=rPs, flatten=True, mode="exclusive_or")
            for T in Ts:
                data = []
                # Pc
                data.append(self["throat.entry_pressure"][T])
                # Element Index
                data.append(T)
                # Element Type (Pore of Throat)
                data.append("throat")
                hq.heappush(queue, data)
        self.invasion_running = [True] * len(self.queue)
        # we have added new clusters that are currently isolated and we
        # need to stop them invading until they merge into an invading
        # cluster
        for c_num in range(len(self.queue)):
            if c_num > initial_num:
                self.invasion_running[c_num] = False

    def _invade_isolated_Ts(self):
        r"""
        Throats that are uninvaded connected to pores that are both invaded
        should be invaded too.
        """
        net = self.project.network
        Ts = net["throat.conns"].copy()
        invaded_Ps = self["pore.invasion_sequence"] > -1
        uninvaded_Ts = self["throat.invasion_sequence"] == -1
        isolated_Ts = np.logical_and(invaded_Ps[Ts[:, 0]], invaded_Ps[Ts[:, 1]])
        isolated_Ts = np.logical_and(isolated_Ts, uninvaded_Ts)
        inv_Pc = self["pore.invasion_pressure"]
        inv_seq = self["pore.invasion_sequence"]
        if np.any(isolated_Ts):
            max_array = Ts[:, 0]
            second_higher = inv_seq[Ts][:, 1] > inv_seq[Ts][:, 0]
            max_array[second_higher] = Ts[:, 1][second_higher]
            mPc = inv_Pc[max_array]
            mSeq = inv_seq[max_array]
            mClu = self["pore.cluster"][max_array]
            self["throat.invasion_pressure"][isolated_Ts] = mPc[isolated_Ts]
            self["throat.invasion_sequence"][isolated_Ts] = mSeq[isolated_Ts]
            self["throat.cluster"][isolated_Ts] = mClu[isolated_Ts]

    def _check_coop(self):
        r"""
        Not implemented in this class
        """
        pass
