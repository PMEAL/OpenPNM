import logging
import heapq as hq
import numpy as np
from tqdm import tqdm
from scipy.stats import rankdata
from collections import namedtuple
from openpnm.utils import Docorator
from openpnm.topotools import find_clusters
from openpnm.algorithms import GenericAlgorithm
from openpnm._skgraph.simulations import (
    bond_percolation,
    find_connected_clusters,
    find_trapped_sites,
)


__all__ = ['InvasionPercolation']


logger = logging.getLogger(__name__)
docstr = Docorator()


@docstr.get_sections(base='IPSettings',
                     sections=['Parameters', 'Other Parameters'])
@docstr.dedent
class IPSettings:
    r"""

    Parameters
    ----------
    %(GenericAlgorithmSettings.parameters)s
    pore_volume : str
        The dictionary key for the pore volume array
    throat_volume : str
        The dictionary key for the throat volume array
    entry_pressure : str
        The dictionary key for the throat capillary pressure

    """
    phase = ''
    pore_volume = 'pore.volume'
    throat_volume = 'throat.volume'
    entry_pressure = 'throat.entry_pressure'


class InvasionPercolation(GenericAlgorithm):
    r"""
    A classic invasion percolation algorithm optimized for speed with numba

    Parameters
    ----------
    network : GenericNetwork
        The Network upon which the invasion will occur

    Notes
    -----
    This algorithm uses a `binary heap <https://en.wikipedia.org/wiki/Binary_heap>`_
    to store a list of all accessible throats, sorted according to entry
    pressure.  This means that item [0] in the heap is the most easily invaded
    throat that is currently accessible by the invading fluid, so looking up
    which throat to invade next is computationally trivial. In order to keep
    the list sorted, adding new throats to the list takes more time; however,
    the heap data structure is very efficient at this.

    """

    def __init__(self, phase, name='ip_#', **kwargs):
        super().__init__(name=name, **kwargs)
        self.settings._update(IPSettings())
        self.settings['phase'] = phase.name
        self['pore.inlets'] = False
        self['pore.outlets'] = False
        self.reset()

    def reset(self):
        self['pore.invasion_sequence'] = -1
        self['throat.invasion_sequence'] = -1
        self['pore.trapped'] = False
        self['throat.trapped'] = False
        self['pore.residual'] = False
        self['throat.residual'] = False
        self.queue = []

    def set_inlets(self, pores=[], mode='add'):
        r"""
        Specifies from which pores the invasion process starts

        Parameters
        ----------
        pores : array_like
            The list of pores from which the invading phase can enter the network
        mode : str
            Controls how the given inlets are added.  Options are:

            ============ ======================================================
            mode         description
            ============ ======================================================
            'add'        Adds the given ``pores`` to list of inlets, while
                         keeping any already defined inlets
            'drop'       Removes given ``pores`` from list of inlets
            'clear'      Removes all currently specified inlets
            'overwrite'  Removes all existing inlets, then adds the given
                         ``pores``
            ============ =======================================================

        """
        if mode == 'add':
            self['pore.inlets'][pores] = True
        elif mode == 'drop':
            self['pore.inlets'][pores] = False
        elif mode == 'clear':
            self['pore.inlets'] = False
        elif mode == 'overwrite':
            self['pore.inlets'] = False
            self['pore.inlets'][pores] = True
        else:
            raise Exception(f'Unrecognized mode {mode}')

    def set_outlets(self, pores=[], mode='add'):
        r"""
        Specifies from which pores the defending fluid exits the domain

        Parameters
        ----------
        pores : array_like
            The list of pores from which the defending fluid exits
        mode : str
            Controls how the given inlets are added.  Options are:

            ============ ======================================================
            mode         description
            ============ ======================================================
            'add'        Adds the given ``pores`` to list of outlets, while
                         keeping any already defined outlets
            'drop'       Removes given ``pores`` from list of outlets
            'clear'      Removes all currently specified outlets
            'overwrite'  Removes all existing outlets, then adds the given
                         ``pores``
            ============ =======================================================

        """
        if mode == 'add':
            self['pore.outlets'][pores] = True
        elif mode == 'drop':
            self['pore.outlets'][pores] = False
        elif mode == 'clear':
            self['pore.outlets'] = False
        elif mode == 'overwrite':
            self['pore.outlets'] = False
            self['pore.outlets'][pores] = True
        else:
            raise Exception(f'Unrecognized mode {mode}')

    def _run_setup(self):
        self['pore.invasion_sequence'][self['pore.inlets']] = 0
        phase = self.project[self.settings['phase']]
        self['throat.entry_pressure'] = phase[self.settings['entry_pressure']]
        # Generated indices into t_entry giving a sorted list
        self['throat.sorted'] = np.argsort(self['throat.entry_pressure'], axis=0)
        self['throat.order'] = 0
        self['throat.order'][self['throat.sorted']] = np.arange(0, self.Nt)
        # Perform initial analysis on input pores
        pores = self['pore.invasion_sequence'] == 0
        Ts = self.project.network.find_neighbor_throats(pores=pores)
        for T in self['throat.order'][Ts]:
            hq.heappush(self.queue, T)

    def run(self, n_steps=None):
        r"""
        Performs the algorithm for the given number of steps

        Parameters
        ----------
        n_steps : int
            The number of throats to invade during the run.  This can be
            used to incrementally invading the network, allowing for
            simulations to occur between each call to ``run``. If ``None``
            (default) then the entire network is invaded.

        """

        # Setup arrays and info
        # This should be called conditionally so that it doesn't overwrite
        # existing data
        self._run_setup()

        if n_steps is None:
            n_steps = np.inf

        if len(self.queue) == 0:
            logger.warn('queue is empty, this network is fully invaded')
            return

        # Create incidence matrix to get neighbor throats later in _run method
        incidence_matrix = self.network.create_incidence_matrix(fmt='csr')
        t_inv, p_inv, p_inv_t = InvasionPercolation._run_accelerated(
            queue=self.queue,
            t_sorted=self['throat.sorted'],
            t_order=self['throat.order'],
            t_inv=self['throat.invasion_sequence'],
            p_inv=self['pore.invasion_sequence'],
            p_inv_t=np.zeros_like(self['pore.invasion_sequence']),
            conns=self.project.network['throat.conns'],
            idx=incidence_matrix.indices,
            indptr=incidence_matrix.indptr,
            n_steps=n_steps
        )

        self['throat.invasion_sequence'] = t_inv
        self['pore.invasion_sequence'] = p_inv
        self['throat.invasion_pressure'] = self['throat.entry_pressure']
        self['pore.invasion_pressure'] = self['throat.entry_pressure'][p_inv_t]
        self['pore.invasion_pressure'][self['pore.invasion_sequence'] == 0] = 0.0

    def apply_trapping(self, n_steps=100):
        r"""
        Analyze which pores and throats are trapped

        Parameters
        ----------
        n_steps : int
            The number of steps to divide the invasion sequence into when
            evaluating trapping. Ideally this number should be equal to the
            number of throats in the network, but this would result in very
            slow computational speeds.  The default is 100, which will incur
            some error but is a good compromise.
        """
        if not np.any(self['pore.outlets']):
            raise Exception('pore outlets must be specified first')
        # Firstly, any pores/throats with inv_seq > outlets are trapped
        N = self['pore.invasion_sequence'][self['pore.outlets']].max()
        self['pore.trapped'][self['pore.invasion_sequence'] > N] = True
        self['throat.trapped'][self['throat.invasion_sequence'] > N] = True
        # Now scan network and find pores/throats disconnected from outlets
        if n_steps is None:
            delta_n = 1
        else:
            delta_n = int(self.Nt/n_steps)
        msg = 'Evaluating trapping'
        for i in tqdm(range(delta_n, int(N), delta_n), msg):
            s, b = find_trapped_sites(conns=self.network.conns,
                                      occupied_sites=self['pore.invasion_sequence'] < i,
                                      outlets=self['pore.outlets'])
            P_trapped = s >= 0
            T_trapped = P_trapped[self.network.conns].any(axis=1)
            # ax = op.topotools.plot_connections(pn, throats=(b >= 0))
            # op.topotools.plot_coordinates(pn, color_by=(s + 10)*(s >= 0), ax=ax, s=200)
            self['pore.trapped'] += P_trapped
            self['throat.trapped'] += T_trapped
        # Set trapped pores and throats to uninvaded and adjust invasion sequence
        self['pore.invasion_sequence'][self['pore.trapped']] = -1
        self['throat.invasion_sequence'][self['throat.trapped']] = -1
        self['pore.invasion_sequence'] = \
            rankdata(self['pore.invasion_sequence'], method='dense') - 2
        self['throat.invasion_sequence'] = \
            rankdata(self['throat.invasion_sequence'], method='dense') - 2

    def apply_trapping_2(self, outlets):
        """
        Identify trapped pores and throats based on algorithm described by
        Y. Masson [1]. It is applied as a post-process and runs an invasion
        percolation algorithm in reverse.

        Parameters
        ----------
        outlets : list
            pore indices indicates where defending fluid can escape through

        Returns
        -------
        Creates a throat array called 'pore.clusters' in the Algorithm
        dictionary. Any positive number is a trapped cluster

        Also creates 2 boolean arrays Np and Nt long called '<element>.trapped'

        Reference
        ---------
        [1] Masson, Y., 2016. A fast two-step algorithm for invasion
        percolation with trapping. Computers & Geosciences, 90, pp.41-48

        Notes
        -----
        Consider the following scenario when running standard IP without
        trapping, three situations can happen after each invasion step:

            * The number of defending clusters stays the same and clusters
              can shrink
            * A cluster of size one is suppressed
            * A cluster is split into multiple clusters

        In reverse the following opposite situations can happen:

            * The number of defending clusters stays the same and clusters
              can grow
            * A cluster of size one is created
            * Mutliple clusters merge into one cluster

        With trapping the reversed rules are adjusted so that only clusters
        that do not connect to a sink can grow and merge. At the point that a
        neighbor connected to a sink is touched the trapped cluster stops
        growing as this is the point of trapping in forward invasion time.

        Logger info displays the invasion sequence and pore index and a message
        with condition number based on the modified trapping rules and the
        assignment of the pore to a given cluster.

        Initially all invaded pores are given cluster label -1
        Outlets / Sinks are given -2

        New clusters that grow into fully trapped clusters are either
        identified at the point of breakthrough or grow from nothing if the
        full invasion sequence is run, they are assigned numbers from 0 up.

        """
        # First see if network is fully invaded
        net = self.project.network
        invaded_ps = self['pore.invasion_sequence'] > -1
        if ~np.all(invaded_ps):
            # Put defending phase into clusters
            clusters = find_clusters(network=net, mask=~invaded_ps)
            # Identify clusters that are connected to an outlet and set to -2
            # -1 is the invaded fluid
            # -2 is the defender fluid able to escape
            # All others now trapped clusters which grow as invasion is reversed
            out_clusters = np.unique(clusters[outlets])
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
        all_neighbors = net.find_neighbor_pores(net.pores(), flatten=False,
                                                include_input=True)
        for un_seq, pore in inv_seq:
            if pore not in outlets and un_seq > 0:  # Skip inlets and outlets
                nc = clusters[all_neighbors[pore]]  # Neighboring clusters
                unique_ns = np.unique(nc[nc != -1])  # Unique Neighbors
                seq_pore = "S:"+str(un_seq)+" P:"+str(pore)
                if np.all(nc == -1):
                    # This is the start of a new trapped cluster
                    clusters[pore] = next_cluster_num
                    next_cluster_num += 1
                    msg = (seq_pore+" C:1 new cluster number: "
                           + str(clusters[pore]))
                    logger.info(msg)
                elif len(unique_ns) == 1:
                    # Grow the only connected neighboring cluster
                    if not stopped_clusters[unique_ns[0]]:
                        clusters[pore] = unique_ns[0]
                        msg = (seq_pore+" C:2 joins cluster number: "
                               + str(clusters[pore]))
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
                            msg = (seq_pore + " C:5 merge clusters: "
                                   + str(c) + " into "+str(new_num))
                            logger.info(msg)

        # And now return clusters
        self['pore.clusters'] = clusters
        logger.info("Number of trapped clusters"
                    + str(np.sum(np.unique(clusters) >= 0)))
        self['pore.trapped'] = self['pore.clusters'] > -1
        trapped_ts = net.find_neighbor_throats(self['pore.trapped'])
        self['throat.trapped'] = np.zeros([net.Nt], dtype=bool)
        self['throat.trapped'][trapped_ts] = True
        self['pore.invasion_sequence'][self['pore.trapped']] = -1
        self['throat.invasion_sequence'][self['throat.trapped']] = -1

    def get_intrusion_data(self):
        r"""
        Get the percolation data as the invader volume or number fraction vs
        the capillary pressure.

        """
        net = self.project.network
        pvols = net[self.settings['pore_volume']]
        tvols = net[self.settings['throat_volume']]
        tot_vol = np.sum(pvols) + np.sum(tvols)
        # Normalize
        pvols /= tot_vol
        tvols /= tot_vol
        # Remove trapped volume
        pvols[self['pore.invasion_sequence'] == -1] = 0.0
        tvols[self['throat.invasion_sequence'] == -1] = 0.0
        pseq = self['pore.invasion_sequence']
        tseq = self['throat.invasion_sequence']
        tPc = self['throat.invasion_pressure']
        pPc = self['pore.invasion_pressure']
        # Change the entry pressure for trapped pores and throats to be 0
        pPc[self['pore.invasion_sequence'] == -1] = 0.0
        tPc[self['throat.invasion_sequence'] == -1] = 0.0
        vols = np.concatenate((pvols, tvols))
        seqs = np.concatenate((pseq, tseq))
        Pcs = np.concatenate((pPc, tPc))
        data = np.rec.fromarrays([seqs, vols, Pcs], formats=['i', 'f', 'f'],
                                 names=['seq', 'vol', 'Pc'])
        data.sort(axis=0, order='seq')
        sat = np.cumsum(data.vol)
        pc_curve = namedtuple('pc_curve', ('Pcap', 'S_tot'))
        data = pc_curve(data.Pc, sat)
        return data

    def _run_accelerated(queue, t_sorted, t_order, t_inv, p_inv, p_inv_t,
                         conns, idx, indptr, n_steps):
        r"""
        Numba-jitted run method for InvasionPercolation class.

        Notes
        -----
        ``idx`` and ``indptr`` are properties are the network's incidence
        matrix, and are used to quickly find neighbor throats.

        Numba doesn't like foreign data types (i.e. GenericNetwork), and so
        ``find_neighbor_throats`` method cannot be called in a jitted method.

        Nested wrapper is for performance issues (reduced OpenPNM import)
        time due to local numba import

        """
        from numba import njit

        @njit
        def wrapper(queue, t_sorted, t_order, t_inv, p_inv, p_inv_t, conns,
                    idx, indptr, n_steps):
            count = 1
            while (len(queue) > 0) and (count < (n_steps + 1)):
                # Find throat at the top of the queue
                t = hq.heappop(queue)
                # Extract actual throat number
                t_next = t_sorted[t]
                t_inv[t_next] = count
                # If throat is duplicated
                while len(queue) > 0 and queue[0] == t:
                    _ = hq.heappop(queue)
                # Find pores connected to newly invaded throat
                Ps = conns[t_next]
                # Remove already invaded pores from Ps
                Ps = Ps[p_inv[Ps] < 0]
                if len(Ps) > 0:
                    p_inv[Ps] = count
                    p_inv_t[Ps] = t_next
                    for i in Ps:
                        Ts = idx[indptr[i]:indptr[i+1]]
                        Ts = Ts[t_inv[Ts] < 0]
                    for i in np.unique(Ts):  # Exclude repeated neighbor throats
                        hq.heappush(queue, t_order[i])
                count += 1
            return t_inv, p_inv, p_inv_t

        return wrapper(queue, t_sorted, t_order, t_inv, p_inv, p_inv_t, conns,
                       idx, indptr, n_steps)


if __name__ == '__main__':
    import openpnm as op
    import matplotlib.pyplot as plt

    pn = op.network.Cubic(shape=[250, 250, 1], spacing=1e-4)
    pn.add_model_collection(op.models.collections.geometry.spheres_and_cylinders)
    pn.regenerate_models()

    water = op.phase.Water(network=pn, name='h2o')
    water.add_model_collection(op.models.collections.physics.standard)
    water.regenerate_models()

    ip = InvasionPercolation(network=pn, phase=water)
    ip.set_inlets(pn.pores('left'))
    ip.set_outlets(pn.pores('right'))
    ip.run()
    # ip.apply_trapping(n_steps=100)

    # ip.plot_intrusion_curve()
    # %%
    fig, ax = plt.subplots(1, 1)
    ax.set_facecolor('grey')
    # ax = op.topotools.plot_coordinates(network=pn, c='w', ax=ax)
    ax = op.topotools.plot_connections(network=pn,
                                       throats=ip['throat.invasion_sequence'] >= 0,
                                       color_by=ip['throat.invasion_sequence'],
                                       ax=ax)























