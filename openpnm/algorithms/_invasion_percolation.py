import logging
import heapq as hq
import numpy as np
from numba import njit, jit
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
    find_trapped_bonds,
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

    def _apply_trapping2(self):
        """
        Apply trapping based on algorithm described by Y. Masson [1].

        References
        ----------
        [1] Masson, Y., 2016. A fast two-step algorithm for invasion
        percolation with trapping. Computers & Geosciences, 90, pp.41-48
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
            out_clusters = np.unique(clusters[self['pore.outlets']])
            for c in out_clusters:
                if c >= 0:
                    clusters[clusters == c] = -2
        else:
            # Go from end
            clusters = -np.ones(net.Np, dtype=int)
            clusters[self['pore.outlets']] = -2
        # Turn into a list for indexing
        inv_seq = np.vstack((self['pore.invasion_sequence'].astype(int),
                             np.arange(0, net.Np, dtype=int))).T
        # Reverse sort list
        inv_seq = inv_seq[inv_seq[:, 0].argsort()][::-1]
        # For all the steps after the inlets are set up to break-through
        # Reverse the sequence and assess the neighbors cluster state
        am = net.create_adjacency_matrix(fmt='csr')
        clusters = reverse_ip(inv_seq, clusters, am.indices, am.indptr)
        # And now return clusters
        self['pore.clusters'] = clusters
        self['pore.trapped'] = self['pore.clusters'] > -1
        trapped_ts = net.find_neighbor_throats(self['pore.trapped'])
        self['throat.trapped'] = False
        self['throat.trapped'][trapped_ts] = True
        self['pore.invasion_sequence'][self['pore.trapped']] = -1
        self['throat.invasion_sequence'][self['throat.trapped']] = -1
        # Apply a fix to remove incorrectly invaded throats
        # TODO: Fix the reverse_ip function to prevent the following problem
        # from occuring to start with
        # NOTE: The following fix is not perfect. It seems that pores on the
        # outlet are not treated as invaded by the clustering, so an escape
        # route that includes invaded outlet pores is considered valid
        tmask = np.stack((self['throat.invasion_sequence'],
                          self['throat.invasion_sequence'])).T
        pmask = self['pore.invasion_sequence'][net.conns]
        hits = ~np.any(pmask == tmask, axis=1)
        self['throat.trapped'][hits] = True
        self['throat.invasion_sequence'][hits] = -1

    def apply_trapping(self, step_size=1, mode='mixed'):
        r"""
        Analyze which pores and throats are trapped

        Parameters
        ----------
        step_size: int
            The number of steps between evaluations of trapping. A value
            of 1 will provide the "True" result, but this would require
            long computational time since a network clustering is performed
            for each step. The default is 1.
        mode : str
            Which algorithm to use when detecting trapped clusters. Options
            are:
            ========== ========================================================
            Mode       Description
            ========== ========================================================
            'reverse'  Uses the method of Masson to compute trapping by
                       scanning the invasion sequence in reverse and detecting
                       disconnected clusters using a disjoint data set. This
                       method is substantially faster than the others.
            'mixed'    Combines site and bond percolation to enforce that the
                       defending fluid can only escape via a path of consisting
                       of uninvaded sites and bonds.
            'bond'     The defending phase can escape via a path of connected
                       bonds regardless of whether any sites are invaded.
            'site'     The defending phase can escapse via path connected sites
                       regardless of whether any bonds are invaded.
            ========== ========================================================

        Returns
        -------
        This function does not return anything. It adjusts the
        ``'pore.invasion_sequence'`` and ``'throat.invasion_sequence'`` arrays
        on the object by setting trapped pores/throats to -1. It also puts
        ``True`` values into the ``'pore.trapped'`` and ``'throat.trapped'``
        arrays.

        Notes
        -----
        Outlet pores must be specified (using ``set_outlets`` or putting
        ``True`` values in ``alg['pore.outlets']``) or else an exception is
        raised.

        """
        if mode == 'reverse':
            self._apply_trapping2()
            return
        if mode == 'reverse2':
            outlets = np.where(self['pore.outlets'])[0]
            am = self.network.create_adjacency_matrix(fmt='csr')
            inv_seq = self['pore.invasion_sequence']
            indices = am.indices
            indptr = am.indptr
            self['pore.trapped'] = reverse2(inv_seq, am.indices, am.indptr, outlets)
            self['pore.invasion_sequence'][self['pore.trapped']] = -1
            trapped_ts = self.network.find_neighbor_throats(self['pore.trapped'])
            self['throat.trapped'] = False
            self['throat.trapped'][trapped_ts] = True
            self['throat.invasion_sequence'][self['throat.trapped']] = -1
            tmask = np.stack((self['throat.invasion_sequence'],
                              self['throat.invasion_sequence'])).T
            pmask = self['pore.invasion_sequence'][self.network.conns]
            hits = ~np.any(pmask == tmask, axis=1)
            self['throat.trapped'][hits] = True
            self['throat.invasion_sequence'][hits] = -1
            return
        # TODO: This could be parallelized with dask since each loop is
        # independent of the others
        N = self['throat.invasion_sequence'].max()
        pseq = self['pore.invasion_sequence']
        tseq = self['throat.invasion_sequence']
        msg = 'Evaluating trapping'
        for i in tqdm(range(0, int(N), step_size), msg):
            if mode == 'bond':
                i += 1
                s, b = bond_percolation(conns=self.network.conns,
                                        occupied_bonds=tseq > i)
            elif mode == 'site':
                s, b = site_percolation(conns=self.network.conns,
                                        occupied_sites=pseq > i)
            elif mode == 'mixed':
                s, b = mixed_percolation(conns=self.network.conns,
                                         occupied_sites=pseq > i,
                                         occupied_bonds=tseq > i)
            clusters = np.unique(s[self['pore.outlets']])
            self['pore.trapped'] += np.isin(s, clusters, invert=True)*(pseq > i)
            self['throat.trapped'] += np.isin(b, clusters, invert=True)*(tseq > i)
        # Set trapped pores/throats to uninvaded and adjust invasion sequence
        self['pore.invasion_sequence'][self['pore.trapped']] = -1
        self['throat.invasion_sequence'][self['throat.trapped']] = -1
        # Set any residual pores within trapped clusters back to untrapped
        self['pore.trapped'][self['pore.residual']] = False
        self['throat.trapped'][self['throat.residual']] = False


@njit
def reverse_ip(inv_seq,
               clusters,
               indices,
               indptr):
    stopped_clusters = np.zeros_like(clusters)
    next_cluster_num = np.max(clusters)+1
    i = -1
    for un_seq, pore in inv_seq:
        i += 1
        un_seq, pore = inv_seq[i]
        neighbors = indices[indptr[pore]:indptr[pore+1]]
        nc = clusters[neighbors]  # Get cluster numbers of all neighbor pores
        unique_ns = np.unique(nc[nc != -1])  # Isolate unique neighbor clusters
        if np.all(nc == -1):
            # The current pore and all of its neighbors have not been visited
            # so it has no cluster label, so let's start growing a new cluster
            # from this pore
            clusters[pore] = next_cluster_num
            next_cluster_num += 1
            print('case 1')
        elif (len(unique_ns) == 1) and (-2 not in unique_ns):
            # If all of the neighbors have the same cluster number, assign that
            # cluster number to the current pore...unless any of the neighbors
            # are 'stopped', in which case mark the current pore as stopped too
            if not stopped_clusters[unique_ns[0]]:
                clusters[pore] = unique_ns[0]
            else:
                clusters[pore] = -2
            print('case 2')
        elif -2 in unique_ns:
            # If any of the neighbors are an outlet, then set current pore to
            # have a cluster number of -2 also (seems like this may be the problem)
            clusters[pore] = -2
            # Label all neighboring pores as stopped
            stopped_clusters[unique_ns[unique_ns > -1]] = True
            print('case 3')
        else:
            print('case 4')
            # If any of the neighbors are set to stopped,
            if np.any(stopped_clusters[unique_ns]):
                clusters[pore] = -2
                # Stop growing all neighboring clusters
                stopped_clusters[unique_ns] = True
            else:
                # Merge multiple un-stopped trapped clusters
                new_num = unique_ns[0]
                clusters[pore] = new_num
                for c in unique_ns:
                    # TODO: This seems like it's a speed bottleneck as it
                    # does not use the disjoint dataset approach
                    clusters[clusters == c] = new_num
    return clusters


@njit
def reverse2(inv_seq, indices, indptr, outlets):
    Np = len(inv_seq)
    sorted_seq = np.vstack((inv_seq.astype(np.int_), np.arange(Np, dtype=np.int_))).T
    sorted_seq = sorted_seq[sorted_seq[:, 0].argsort()][::-1]
    cluster = -np.ones(Np, dtype=np.int_)
    trapped_pores = np.zeros(Np, dtype=np.bool_)
    trapped_clusters = np.zeros(Np, dtype=np.bool_)
    cluster_map = qupc_initialize(Np)
    next_cluster_num = 0
    i = -1
    for step, pore in sorted_seq:
        i += 1
        step, pore = sorted_seq[i, :]
        n = indices[indptr[pore]:indptr[pore+1]]
        nc = cluster_map[cluster[n]][inv_seq[n] > step]
        if len(nc) == 0:
            # Found an isolated pore, start a new cluster
            cluster[pore] = next_cluster_num
            # If pore is an outlet then note cluster as no longer trapped
            if pore in outlets:
                trapped_clusters[next_cluster_num] = False
            else:  # Otherwise note this cluster as being a trapped cluster
                trapped_clusters[next_cluster_num] = True
                # Note this pore as trapped as well
                trapped_pores[pore] = True
            # Increment cluster number for next time
            next_cluster_num += 1
        elif len(np.unique(nc)) == 1:
            c = np.unique(nc)[0]
            # Neighbor have one unique cluster number, so assign it to current pore
            cluster[pore] = c
            # If pore is an outlet then note cluster as no longer trapped
            if pore in outlets:
                trapped_clusters[c] = False
                # Also set all joined clusters to not trapped
                cluster_map = qupc_reduce(cluster_map)
                hits = np.where(cluster_map == cluster_map[c])[0]
                trapped_clusters[hits] = False
            # If this cluster number is part of a trapped cluster then
            # mark pore as trapped
            if trapped_clusters[c]:
                trapped_pores[pore] = True
        elif len(np.unique(nc)) > 1:
            cluster[pore] = min(np.unique(nc))
            # Merge all clusters into a single cluster
            for c in nc:
                qupc_update(cluster_map, c, min(np.unique(nc)))
            cluster_map = qupc_reduce(cluster_map)
            # If all neighboring clusters are trapped, then set current pore to
            # trapped as well
            if np.all(trapped_clusters[nc]):
                trapped_pores[pore] = True
            else:  # Otherwise set all neighbor clusters to untrapped!
                trapped_clusters[nc] = False
    return trapped_pores


# %%

@njit
def qupc_initialize(size):
    return np.arange(size, dtype=np.int_)


@njit
def qupc_update(arr, ind, val):
    if ind == val:
        arr[ind] = val
    else:
        # Update array and do path compression simultaneously
        while arr[ind] != val:
            arr[ind] = arr[val]
            ind = val
            val = arr[val]
    return arr


def qupc_compress(arr):
    temp = rankdata(arr, method='dense')
    arr[:] = temp
    arr -= 1
    return arr


@njit
def qupc_reduce(arr):
    for i in range(len(arr)-1, 0, -1):
        arr[i] = arr[arr[i]]
    return arr


if 0:
    a = qupc_initialize(10)
    qupc_update(a, 4, 2)
    qupc_update(a, 7, 4)
    qupc_update(a, 9, 6)
    qupc_update(a, 6, 2)
    qupc_update(a, 5, 9)
    assert np.all(a == [0, 1, 2, 3, 2, 6, 2, 2, 8, 6])
    qupc_reduce(a, compress=False)
    assert np.all(a == [0, 1, 2, 3, 2, 2, 2, 2, 8, 2])
    qupc_update(a, 9, 9)
    qupc_update(a, 0, 1)
    qupc_update(a, 8, 0)
    assert np.all(a == [1, 1, 2, 3, 2, 2, 2, 2, 1, 9])
    qupc_reduce(a, compress=True)
    assert np.all(a == [0, 0, 1, 2, 1, 1, 1, 1, 0, 3])
    print(a)

# %%
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
        # Find pores connected to newly invaded throat from am in coo format
        Ps = conns[t_next]
        # Remove already invaded pores from Ps
        Ps = Ps[p_inv[Ps] < 0]
        # If either of the neighboring pores are uninvaded (-1), set it to
        # invaded and add its neighboring throats to the queue
        if len(Ps) > 0:
            p_inv[Ps] = count
            p_inv_t[Ps] = t_next
            for i in Ps:
                # Get neighboring throat numbers from im in csr format
                Ts = idx[indptr[i]:indptr[i+1]]
                # Keep only throats which are uninvaded
                Ts = Ts[t_inv[Ts] < 0]
            for i in Ts:  # Add throat to the queue
                hq.heappush(queue, t_order[i])
        count += 1
    return t_inv, p_inv, p_inv_t


if __name__ == '__main__':
    import openpnm as op
    import matplotlib.pyplot as plt

    for seed in [2]:
        np.random.seed(seed)
        Nx, Ny, Nz = 25, 25, 1
        pn = op.network.Cubic(shape=[Nx, Ny, Nz], spacing=1e-4)
        pn.add_model_collection(op.models.collections.geometry.spheres_and_cylinders)
        pn.regenerate_models()
        pn['pore.volume@left'] = 0.0
        # op.topotools.trim(pn, pores=[380, 395])

        water = op.phase.Water(network=pn, name='h2o')
        water.add_model_collection(op.models.collections.physics.standard)
        water.regenerate_models()

        p_residual = np.random.randint(0, Nx*Ny, int(Nx*Ny/10))
        t_residual = pn.find_neighbor_throats(p_residual, mode='or')

        ip = InvasionPercolation(network=pn, phase=water)
        ip.set_inlets(pn.pores('left'))
        ip.run()
        ip.set_outlets(pn.pores('right'))
        ip.apply_trapping(step_size=1, mode='reverse2')

        ip2 = InvasionPercolation(network=pn, phase=water)
        ip2.set_inlets(pn.pores('left'))
        ip2.run()
        ip2.set_outlets(pn.pores('right'))
        ip2.apply_trapping(step_size=1, mode='mixed')

        assert np.all(ip['pore.trapped'] == ip2['pore.trapped'])

    # %%
    if 1:
        pseq = np.copy(ip['pore.invasion_sequence'])
        tseq = np.copy(ip['throat.invasion_sequence'])
        ax = op.topotools.plot_connections(pn, tseq>=0, linestyle='--', c='b', linewidth=3)
        op.topotools.plot_coordinates(pn, pseq>=0, c='r', marker='x', markersize=100, ax=ax)
        op.topotools.plot_coordinates(pn, pseq<0, c='c', marker='o', markersize=100, ax=ax)

    # %%
    if 0:
        drn = op.algorithms.Drainage(network=pn, phase=water)
        drn.set_inlets(pn.pores('left'))
        # drn.set_residual(pores=p_residual, throats=t_residual)
        pressures = np.unique(ip['pore.invasion_pressure'])
        # pressures = np.logspace(np.log10(0.1e3), np.log10(2e4), 100)
        drn.run(pressures=pressures)
        drn.set_outlets(pn.pores('right'))
        drn.apply_trapping(mode='mixed')

    # %%
    if 0:
        fig, ax = plt.subplots(1, 1)
        ax.step(*ip.pc_curve(), 'b', where='post')
        ax.step(*drn.pc_curve(), 'r--', where='post')
        ax.set_ylim([0, 1])


    # %%
    if 0:
        ip = ip2
        from matplotlib import animation
        import openpnm as op
        pseq = ip['pore.invasion_sequence']
        pseq[pseq == -1] = pseq.max() + 1
        tseq = ip['throat.invasion_sequence']
        tseq[tseq == -1] = tseq.max() + 1
        for j, i in enumerate(tqdm(np.unique(tseq)[:-1])):
            ax = op.topotools.plot_connections(pn, tseq<=i, linestyle='--', c='r', linewidth=3)
            op.topotools.plot_coordinates(pn, pseq<=i, c='b', marker='x', markersize=100, ax=ax)
            op.topotools.plot_coordinates(pn, ip['pore.trapped'], c='c', marker='o', markersize=100, ax=ax)
            plt.savefig(f"{str(j).zfill(3)}.png")
            plt.close()


















