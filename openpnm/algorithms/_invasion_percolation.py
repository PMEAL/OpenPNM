import logging
import heapq as hq
import numpy as np
from numba import njit
from numba.typed import List
from tqdm import tqdm
from scipy.stats import rankdata
from collections import namedtuple
from openpnm.utils import Docorator
from openpnm.topotools import find_clusters
from openpnm.algorithms import GenericAlgorithm
from openpnm._skgraph.simulations import (
    bond_percolation,
    site_percolation,
    mixed_percolation,
    find_connected_clusters,
    find_trapped_sites,
    find_trapped_bonds
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

    def set_residual(self, pores=None, throats=None, mode='add'):
        if mode == 'add':
            if pores is not None:
                self['pore.residual'][pores] = True
            if throats is not None:
                self['throat.residual'][throats] = True
        elif mode == 'drop':
            if pores is not None:
                self['pore.residual'][pores] = False
            if throats is not None:
                self['throat.residual'][throats] = False
        elif mode == 'clear':
            if pores is not None:
                self['pore.residual'] = False
            if throats is not None:
                self['throat.residual'] = False
        elif mode == 'overwrite':
            if pores is not None:
                self['pore.residual'] = False
                self['pore.residual'][pores] = True
            if throats is not None:
                self['throat.residual'] = False
                self['throat.residual'][throats] = True

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
        # TODO: This should be called conditionally so that it doesn't
        # overwrite existing data when doing a few steps at a time
        self._run_setup()

        # TODO: The following line is supposed to be numba's new list, but the
        # heap does not work with this
        # self.queue = List(self.queue)

        if n_steps is None:
            n_steps = np.inf

        # Create incidence matrix for use in _run_accelerated which is jit
        im = self.network.create_incidence_matrix(fmt='csr')

        # Perform initial analysis on input pores
        Ts = self.project.network.find_neighbor_throats(pores=self['pore.inlets'])
        t_start = self['throat.order'][Ts]
        t_inv, p_inv, p_inv_t = \
            _run_accelerated(
                t_start=t_start.astype(np.int32),
                t_sorted=self['throat.sorted'],
                t_order=self['throat.order'],
                t_inv=self['throat.invasion_sequence'],
                p_inv=self['pore.invasion_sequence'],
                p_inv_t=np.zeros_like(self['pore.invasion_sequence']),
                conns=self.project.network['throat.conns'],
                idx=im.indices,
                indptr=im.indptr,
                n_steps=n_steps)

        # Transfer results onto algorithm object
        self['throat.invasion_sequence'] = t_inv
        self['pore.invasion_sequence'] = p_inv
        self['throat.invasion_pressure'] = self['throat.entry_pressure']
        self['pore.invasion_pressure'] = self['throat.entry_pressure'][p_inv_t]
        # Set invasion pressure of inlets to 0
        self['pore.invasion_pressure'][self['pore.invasion_sequence'] == 0] = 0.0
        # Set invasion sequence and pressure of any residual pores/throats to 0
        self['throat.invasion_sequence'][self['throat.residual']] = 0
        self['pore.invasion_sequence'][self['pore.residual']] = 0
        self['pore.invasion_pressure'][self['pore.residual']] = 0

    def _run_setup(self):
        self['pore.invasion_sequence'][self['pore.inlets']] = 0
        # self['pore.invasion_sequence'][self['pore.residual']] = 0
        # self['throat.invasion_sequence'][self['throat.residual']] = 0
        # Set throats between inlets as trapped
        Ts = self.network.find_neighbor_throats(self['pore.inlets'], mode='xnor')
        self['throat.trapped'][Ts] = True
        # Get throat capillary pressures from phase and update
        phase = self.project[self.settings['phase']]
        self['throat.entry_pressure'] = phase[self.settings['entry_pressure']]
        self['throat.entry_pressure'][self['throat.residual']] = 0.0
        # Generated indices into t_entry giving a sorted list
        self['throat.sorted'] = np.argsort(self['throat.entry_pressure'], axis=0)
        self['throat.order'] = 0
        self['throat.order'][self['throat.sorted']] = np.arange(0, self.Nt)

    def pc_curve(self):
        r"""
        Get the percolation data as the invader volume vs the capillary
        pressure.

        """
        net = self.project.network
        pvols = net[self.settings['pore_volume']]
        tvols = net[self.settings['throat_volume']]
        tot_vol = np.sum(pvols) + np.sum(tvols)
        # Normalize
        pvols /= tot_vol
        tvols /= tot_vol
        # Remove trapped volume
        pmask = self['pore.invasion_sequence'] > -1
        tmask = self['throat.invasion_sequence'] > -1
        pvols = pvols[pmask]
        tvols = tvols[tmask]
        pseq = self['pore.invasion_sequence'][pmask]
        tseq = self['throat.invasion_sequence'][tmask]
        pPc = self['pore.invasion_pressure'][pmask]
        tPc = self['throat.invasion_pressure'][tmask]
        vols = np.concatenate((pvols, tvols))
        seqs = np.concatenate((pseq, tseq))
        Pcs = np.concatenate((pPc, tPc))
        data = np.rec.fromarrays([seqs, vols, Pcs],
                                 formats=['i', 'f', 'f'],
                                 names=['seq', 'vol', 'Pc'])
        data.sort(axis=0, order='seq')
        pc_curve = namedtuple('pc_curve', ('pc', 'snwp'))
        data = pc_curve(data.Pc, np.cumsum(data.vol))
        return data

    def apply_trapping2(self):
        """
        Apply trapping based on algorithm described by Y. Masson [1].

        References
        ----------
        [1] Masson, Y., 2016. A fast two-step algorithm for invasion
        percolation with trapping. Computers & Geosciences, 90, pp.41-48
        """
        # First see if network is fully invaded
        net = self.project.network
        outlets = self['pore.outlets']
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
        all_neighbors = net.find_neighbor_pores(net.pores(),
                                                flatten=False,
                                                include_input=True)
        clusters = reverse_ip(inv_seq,
                              clusters,
                              stopped_clusters,
                              next_cluster_num,
                              all_neighbors,
                              outlets)
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

    def apply_trapping(self, step_size=10, mode='mixed'):
        r"""
        Analyze which pores and throats are trapped using mixed percolation.

        Mixed percolation is crucial since it ensures that defending phase
        must escape via a connected path of throats *and* pores. Using bond
        percolation erroneously allows via two throats even if their shared
        pore has been invaded.

        Parameters
        ----------
        step_size: int
            The number of steps between between evaluations of trapping.
            A value of 1 will provide the "True" result, but this would
            require long computational time since a network clustering is
            performed for each step. The default is 10, which will incur
            some error (pores and throats are identified as invaded that
            are actually trapped), but is a good compromise.

        Returns
        -------
        This function does not return anything. It adjusts the
        ``'pore.invasion_sequence'`` and ``'throat.invasion_sequence'`` arrays
        on the object by setting trapped pores/throats to -1, and adjusting
        the sequence values to be contiguous.  It also puts ``True`` values
        into the ``'pore.trapped'`` and ``'throat.trapped'`` arrays

        Notes
        -----
        Outlet pores must be specified (using ``set_outlets`` or putting
        ``True`` values in ``alg['pore.outlets']``) or else an exception is
        raised.

        """
        if not np.any(self['pore.outlets']):
            raise Exception('pore outlets must be specified first')
        # TODO: This could be parallelized with dask since each loop is
        # independent of the others
        N = self['throat.invasion_sequence'].max()
        pseq = self['pore.invasion_sequence']
        tseq = self['throat.invasion_sequence']
        msg = 'Evaluating trapping'
        for i in tqdm(range(0, int(N), step_size), msg):
            # Performing cluster labeling
            if mode == 'bond':
                s, b = bond_percolation(conns=self.network.conns,
                                        occupied_bonds=tseq > i)
            elif mode == 'site':
                s, b = site_percolation(conns=self.network.conns,
                                        occupied_sites=pseq > i)
            elif mode == 'mixed':
                s, b = mixed_percolation(conns=self.network.conns,
                                         occupied_sites=pseq > i,
                                         occupied_bonds=tseq > i)
            # Note disconnected pores and throats as trapped
            s_def = np.unique(s[self['pore.outlets']])
            self['pore.trapped'] += np.isin(s, s_def, invert=True)*(s >= 0)
            self['throat.trapped'] += np.isin(b, s_def, invert=True)*(b >= 0)
        if 0:
            i += 1
            # %%
            s, b = mixed_percolation(conns=self.network.conns, occupied_sites=pseq > i, occupied_bonds=tseq > i)
            ax = op.topotools.plot_connections(pn, b>=0, color_by=b, linewidth=3)
            op.topotools.plot_connections(pn, tseq<=i, linestyle='--', c='r', linewidth=3, ax=ax)
            op.topotools.plot_coordinates(pn, s>=0, c='g', markersize=100, ax=ax)
            op.topotools.plot_coordinates(pn, pseq<=i, c='r', marker='x', markersize=100, ax=ax)
        # %%

        # Set trapped pores/throats to uninvaded and adjust invasion sequence
        self['pore.invasion_sequence'][self['pore.trapped']] = -1
        self['throat.invasion_sequence'][self['throat.trapped']] = -1
        # Set any residual pores within trapped clusters back to untrapped
        self['pore.trapped'][self['pore.residual']] = False
        self['throat.trapped'][self['throat.residual']] = False
        # The -2 is to shift uninvaded pores to -1 and initially invaded to 0
        # self['pore.invasion_sequence'] = \
        #     rankdata(self['pore.invasion_sequence'], method='dense') - 2
        # self['throat.invasion_sequence'] = \
        #     rankdata(self['throat.invasion_sequence'], method='dense') - 2


@njit
def reverse_ip(inv_seq,
               clusters,
               stopped_clusters,
               next_cluster_num,
               all_neighbors,
               outlets):
    for un_seq, pore in inv_seq:
        if pore not in outlets and un_seq > 0:  # Skip inlets and outlets
            nc = clusters[all_neighbors[pore]]  # Neighboring clusters
            unique_ns = np.unique(nc[nc != -1])  # Unique Neighbors
            if np.all(nc == -1):
                # This is the start of a new trapped cluster
                clusters[pore] = next_cluster_num
                next_cluster_num += 1
            elif len(unique_ns) == 1:
                # Grow the only connected neighboring cluster
                if not stopped_clusters[unique_ns[0]]:
                    clusters[pore] = unique_ns[0]
                else:
                    clusters[pore] = -2
            elif -2 in unique_ns:
                # We have reached a sink neighbor, stop growing cluster
                clusters[pore] = -2
                # Stop growth and merging
                stopped_clusters[unique_ns[unique_ns > -1]] = True
            else:  # We might be able to do some merging
                # Check if any stopped clusters are neighbors
                if np.any(stopped_clusters[unique_ns]):
                    clusters[pore] = -2
                    # Stop growing all neighboring clusters
                    stopped_clusters[unique_ns] = True
                else:
                    # Merge multiple un-stopped trapped clusters
                    new_num = unique_ns[0]
                    clusters[pore] = new_num
                    for c in unique_ns:
                        clusters[clusters == c] = new_num
    return clusters


@njit
def _run_accelerated(t_start, t_sorted, t_order, t_inv, p_inv, p_inv_t,
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
    queue = list(t_start)
    hq.heapify(queue)
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


if __name__ == '__main__':
    import openpnm as op
    import matplotlib.pyplot as plt

    np.random.seed(0)
    Nx, Ny, Nz = 10, 10, 1
    pn = op.network.Cubic(shape=[Nx, Ny, Nz], spacing=1e-4)
    pn.add_model_collection(op.models.collections.geometry.spheres_and_cylinders)
    pn.regenerate_models()
    pn['pore.volume@left'] = 0.0

    water = op.phase.Water(network=pn, name='h2o')
    water.add_model_collection(op.models.collections.physics.standard)
    water.regenerate_models()

    p_residual = np.random.randint(0, Nx*Ny, int(Nx*Ny/10))
    t_residual = pn.find_neighbor_throats(p_residual, mode='or')
    # %%
    ip = InvasionPercolation(network=pn, phase=water)
    ip.set_inlets(pn.pores('left'))
    # ip.set_residual(pores=p_residual, throats=t_residual)
    ip.run()
    ip.set_outlets(pn.pores('right'))
    # ip.apply_trapping(step_size=1, mode='mixed')
    # ip['pore.trapped_1'] = np.copy(ip['pore.trapped'])
    # ip['throat.trapped_1'] = np.copy(ip['throat.trapped'])
    # ip.reset()
    # ip.run()
    ip.apply_trapping2()

    # %%
    drn = op.algorithms.Drainage(network=pn, phase=water)
    drn.set_inlets(pn.pores('left'))
    # drn.set_residual(pores=p_residual, throats=t_residual)
    pressures = np.unique(ip['pore.invasion_pressure'])
    # pressures = np.logspace(np.log10(0.1e3), np.log10(2e4), 100)
    drn.run(pressures=pressures)
    drn.set_outlets(pn.pores('right'))
    drn.apply_trapping(pressures=None)

    # %%
    if 0:
        fig, ax = plt.subplots(1, 1)
        ax.step(*ip.pc_curve(), 'b', where='post')
        ax.step(*drn.pc_curve(), 'r--', where='post')
        ax.set_ylim([0, 1])

    # %%
    if 1:
        pseq = np.copy(ip['pore.invasion_sequence'])
        i = pseq.max()
        pseq[pseq == -1] = pseq.max() + 1
        tseq = np.copy(ip['throat.invasion_sequence'])
        tseq[tseq == -1] = tseq.max() + 1
        ax = op.topotools.plot_connections(pn, tseq<=i, linestyle='--', c='r', linewidth=3, ax=ax)
        op.topotools.plot_coordinates(pn, pseq<=i, c='r', marker='x', markersize=100, ax=ax)

    # %%
    if 0:
        from matplotlib import animation
        import openpnm as op
        pseq = ip['pore.invasion_sequence']
        pseq[pseq == -1] = pseq.max() + 1
        tseq = ip['throat.invasion_sequence']
        tseq[tseq == -1] = tseq.max() + 1
        for i in tqdm(np.unique(pseq)[:-1]):
            s, b = mixed_percolation(conns=ip.network.conns, occupied_sites=pseq > i, occupied_bonds=tseq > i)
            ax = op.topotools.plot_connections(pn, b>=0, color_by=b, linewidth=3)
            op.topotools.plot_connections(pn, tseq<=i, linestyle='--', c='r', linewidth=3, ax=ax)
            op.topotools.plot_coordinates(pn, s>=0, c='g', markersize=100, ax=ax)
            op.topotools.plot_coordinates(pn, pseq<=i, c='r', marker='x', markersize=100, ax=ax)
            plt.savefig(f"{str(i).zfill(3)}.png")
            plt.close()


















