import heapq as hq
import logging
from collections import namedtuple

import numpy as np
from numba import jit, njit
from tqdm.auto import tqdm

from openpnm._skgraph.queries import qupc_initialize, qupc_reduce, qupc_update
from openpnm._skgraph.simulations import bond_percolation, site_percolation
from openpnm.algorithms import Algorithm
from openpnm.utils import Docorator

__all__ = [
    'InvasionPercolation',
]


logger = logging.getLogger(__name__)
docstr = Docorator()


@docstr.get_sections(base='IPSettings',
                     sections=['Parameters', 'Other Parameters'])
@docstr.dedent
class IPSettings:
    r"""

    Parameters
    ----------
    %(AlgorithmSettings.parameters)s
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


class InvasionPercolation(Algorithm):
    r"""
    A classic invasion percolation algorithm optimized for speed with numba

    Parameters
    ----------
    network : Network
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

    def __init__(self, phase, name='ip_?', **kwargs):
        super().__init__(name=name, **kwargs)
        self.settings._update(IPSettings())
        self.settings['phase'] = phase.name
        self['pore.bc.inlet'] = False
        self['pore.bc.outlet'] = False
        self.reset()

    def reset(self):
        self['pore.invasion_sequence'] = -1
        self['throat.invasion_sequence'] = -1
        self['pore.trapped'] = False
        self['throat.trapped'] = False
        # self['pore.residual'] = False
        # self['throat.residual'] = False

    def _set_residual(self, pores=None, throats=None, mode='add'):  # pragma: no cover
        raise NotImplementedError("The ability to add residual nwp is not ready yet")
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

    def set_inlet_BC(self, pores=None, mode='add'):
        r"""
        Specifies which pores are treated as inlets for the invading phase

        Parameters
        ----------
        pores : ndarray
            The indices of the pores from which the invading fluid invasion
            should start
        mode : str or list of str, optional
            Controls how the boundary conditions are applied. Options are:

            ============ =====================================================
            mode         meaning
            ============ =====================================================
            'add'        (default) Adds the supplied boundary conditions to
                         the given locations. Raises an exception if values
                         of any type already exist in the given locations.
            'overwrite'  Adds supplied boundary conditions to the given
                         locations, including overwriting conditions of the
                         given type or any other type that may be present in
                         the given locations.
            'remove'     Removes boundary conditions of the specified type
                         from the specified locations. If ``bctype`` is not
                         specified then *all* types are removed. If no
                         locations are given then values are remvoed from
                         *all* locations.
            ============ =====================================================

            If a list of strings is provided, then each mode in the list is
            handled in order, so that ``['remove', 'add']`` will give the same
            results add ``'overwrite'``.

        """
        self.set_BC(pores=pores, bcvalues=True, bctype='inlet', mode=mode)
        self.reset()
        self['pore.invasion_sequence'][self['pore.bc.inlet']] = 0

    def set_outlet_BC(self, pores=None, mode='add'):
        r"""
        Specifies which pores are treated as outlets for the defending phase

        This must be specified if trapping is to be considered.

        Parameters
        ----------
        pores : ndarray
            The indices of the pores from which the defending fluid exits the
            domain
        mode : str or list of str, optional
            Controls how the boundary conditions are applied. Options are:

            ============ =====================================================
            mode         meaning
            ============ =====================================================
            'add'        (default) Adds the supplied boundary conditions to
                         the given locations. Raises an exception if values
                         of any type already exist in the given locations.
            'overwrite'  Adds supplied boundary conditions to the given
                         locations, including overwriting conditions of the
                         given type or any other type that may be present in
                         the given locations.
            'remove'     Removes boundary conditions of the specified type
                         from the specified locations. If ``bctype`` is not
                         specified then *all* types are removed. If no
                         locations are given then values are remvoed from
                         *all* locations.
            ============ =====================================================

            If a list of strings is provided, then each mode in the list is
            handled in order, so that ``['remove', 'add']`` will give the same
            results add ``'overwrite'``.

        """
        self.set_BC(pores=pores, bcvalues=True, bctype='outlet', mode=mode)

    def run(self):
        r"""
        Performs the algorithm for the given number of steps

        """
        # Setup arrays and info
        # TODO: This should be called conditionally so that it doesn't
        # overwrite existing data when doing a few steps at a time
        self._run_setup()
        n_steps = np.inf

        # Create incidence matrix for use in _run_accelerated which is jit
        im = self.network.create_incidence_matrix(fmt='csr')

        # Perform initial analysis on input pores
        Ts = self.project.network.find_neighbor_throats(pores=self['pore.bc.inlet'])
        t_start = self['throat.order'][Ts]
        t_inv, p_inv, p_inv_t = \
            _run_accelerated(
                t_start=t_start,
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
        # self['throat.invasion_sequence'][self['throat.residual']] = 0
        # self['pore.invasion_sequence'][self['pore.residual']] = 0

    def _run_setup(self):
        self['pore.invasion_sequence'][self['pore.bc.inlet']] = 0
        # self['pore.invasion_sequence'][self['pore.residual']] = 0
        # self['throat.invasion_sequence'][self['throat.residual']] = 0
        # Set throats between inlets as trapped
        Ts = self.network.find_neighbor_throats(self['pore.bc.inlet'], mode='xnor')
        self['throat.trapped'][Ts] = True
        # Get throat capillary pressures from phase and update
        phase = self.project[self.settings['phase']]
        self['throat.entry_pressure'] = phase[self.settings['entry_pressure']]
        # self['throat.entry_pressure'][self['throat.residual']] = 0.0
        # Generated indices into t_entry giving a sorted list
        self['throat.sorted'] = np.argsort(self['throat.entry_pressure'], axis=0)
        self['throat.order'] = 0
        self['throat.order'][self['throat.sorted']] = np.arange(0, self.Nt)

    def pc_curve(self):
        r"""
        Get the percolation data as the non-wetting phase saturation vs the
        capillary pressure.

        """
        net = self.project.network
        pvols = net[self.settings['pore_volume']]
        tvols = net[self.settings['throat_volume']]
        tot_vol = np.sum(pvols) + np.sum(tvols)
        # Normalize
        pvols /= tot_vol
        tvols /= tot_vol
        # Remove trapped volume
        pmask = self['pore.invasion_sequence'] >= 0
        tmask = self['throat.invasion_sequence'] >= 0
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

    def apply_trapping(self):
        r"""
        Adjusts the invasion sequence of pores and throats that are trapped.

        This method uses the reverse invasion percolation procedure outlined
        by Masson [1].

        Returns
        -------
        This function does not return anything. It adjusts the
        ``'pore.invasion_sequence'`` and ``'throat.invasion_sequence'`` arrays
        on the object by setting trapped pores/throats to ``ninf``. It also puts
        ``True`` values into the ``'pore.trapped'`` and ``'throat.trapped'``
        arrays.

        Notes
        -----
        Outlet pores must be specified (using ``set_outlet_BC`` or putting
        ``True`` values in ``alg['pore.bc.outlet']``) or else an exception is
        raised.

        References
        ----------
        [1] Masson, Y. https://doi.org/10.1016/j.cageo.2016.02.003

        """
        outlets = np.where(self['pore.bc.outlet'])[0]
        am = self.network.create_adjacency_matrix(fmt='csr')
        inv_seq = self['pore.invasion_sequence']
        self['pore.trapped'] = _find_trapped_pores(inv_seq, am.indices,
                                                   am.indptr, outlets)
        # Update invasion sequence
        self['pore.invasion_sequence'][self['pore.trapped']] = -1
        # Find which throats are trapped, including throats which were invaded
        # after both of it's pores were invaded (hence have a unique invasion
        # sequence number).
        pmask = self['pore.invasion_sequence'][self.network.conns]
        tmask = np.stack((self['throat.invasion_sequence'],
                          self['throat.invasion_sequence'])).T
        hits = ~np.any(pmask == tmask, axis=1)
        self['throat.trapped'] = hits
        self['throat.invasion_sequence'][hits] = -1
        # Make some adjustments
        Pmask = self['pore.invasion_sequence'] < 0
        Tmask = self['throat.invasion_sequence'] < 0
        self['pore.invasion_sequence'] = \
            self['pore.invasion_sequence'].astype(float)
        self['pore.invasion_sequence'][Pmask] = np.inf
        self['throat.invasion_sequence'] = \
            self['throat.invasion_sequence'].astype(float)
        self['throat.invasion_sequence'][Tmask] = np.inf

    def _apply_trapping_slow(self, step_size=1, mode='mixed'):  # pragma: no cover
        # TODO: Make sure this function actually works and keep it for debugging
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
            clusters = np.unique(s[self['pore.bc.outlet']])
            self['pore.trapped'] += np.isin(s, clusters, invert=True)*(pseq > i)
            self['throat.trapped'] += np.isin(b, clusters, invert=True)*(tseq > i)
        # Set trapped pores/throats to uninvaded and adjust invasion sequence
        self['pore.invasion_sequence'][self['pore.trapped']] = -1
        self['throat.invasion_sequence'][self['throat.trapped']] = -1
        # Set any residual pores within trapped clusters back to untrapped
        # self['pore.trapped'][self['pore.residual']] = False
        # self['throat.trapped'][self['throat.residual']] = False


@jit(forceobj=True)
def _find_trapped_pores(inv_seq, indices, indptr, outlets):
    Np = len(inv_seq)
    sorted_seq = np.vstack((inv_seq.astype(np.int_), np.arange(Np, dtype=np.int_))).T
    sorted_seq = sorted_seq[sorted_seq[:, 0].argsort()][::-1]
    cluster = -np.ones(Np, dtype=np.int_)
    trapped_pores = np.zeros(Np, dtype=bool)
    trapped_clusters = np.zeros(Np, dtype=bool)
    # cluster_map = qupc_initialize(Np)
    cluster_map = np.arange(Np, dtype=np.int_)
    next_cluster_num = 0
    i = -1
    for step, pore in sorted_seq:
        i += 1
        step, pore = sorted_seq[i, :]
        n = indices[indptr[pore]:indptr[pore+1]]
        nc = cluster_map[cluster[n]][inv_seq[n] > step]
        nc_uniq = np.unique(nc)
        if nc.size == 0:
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
        elif nc_uniq.size == 1:
            c = nc_uniq[0]
            # Neighbors have one unique cluster number, so assign it to current pore
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
        elif nc_uniq.size > 1:
            cluster[pore] = min(nc_uniq)
            # Merge all clusters into a single cluster
            for c in nc:
                qupc_update(cluster_map, c, min(nc_uniq))
            cluster_map = qupc_reduce(cluster_map)
            # If all neighboring clusters are trapped, then set current pore to
            # trapped as well
            if np.all(trapped_clusters[nc]):
                trapped_pores[pore] = True
            else:  # Otherwise set all neighbor clusters to untrapped!
                trapped_clusters[nc] = False
    return trapped_pores


@njit
def _run_accelerated(t_start, t_sorted, t_order, t_inv, p_inv, p_inv_t,
                     conns, idx, indptr, n_steps):  # pragma: no cover
    r"""
    Numba-jitted run method for InvasionPercolation class.

    Notes
    -----
    ``idx`` and ``indptr`` are properties are the network's incidence
    matrix, and are used to quickly find neighbor throats.

    Numba doesn't like foreign data types (i.e. Network), and so
    ``find_neighbor_throats`` method cannot be called in a jitted method.

    Nested wrapper is for performance issues (reduced OpenPNM import)
    time due to local numba import

    """
    # TODO: The following line is supposed to be numba's new list, but the
    # heap does not work with this
    # queue = List(t_start)
    queue = list(t_start)
    hq.heapify(queue)
    count = 1
    while count < (n_steps + 1):
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
        if len(queue) == 0:
            break
    return t_inv, p_inv, p_inv_t


# %%
if __name__ == '__main__':
    import matplotlib.pyplot as plt

    import openpnm as op

    np.random.seed(0)
    Nx, Ny, Nz = 25, 25, 1
    pn = op.network.Cubic(shape=[Nx, Ny, Nz], spacing=1e-4)
    pn.add_model_collection(op.models.collections.geometry.spheres_and_cylinders)
    pn.regenerate_models()
    pn['pore.volume@left'] = 0.0
    # op.topotools.trim(pn, pores=[380, 395])

    water = op.phase.Water(network=pn, name='h2o')
    water.add_model_collection(op.models.collections.physics.standard)
    water.regenerate_models()

    ip = InvasionPercolation(network=pn, phase=water)
    ip.set_inlet_BC(pn.pores('left'))
    ip.run()
    ip.set_outlet_BC(pn.pores('right'))
    # ip.apply_trapping()

    # %%
    if 0:
        pseq = np.copy(ip['pore.invasion_sequence'])
        tseq = np.copy(ip['throat.invasion_sequence'])
        ax = op.topotools.plot_connections(pn, tseq >= 0, linestyle='--', c='b',
                                           linewidth=3)
        op.topotools.plot_coordinates(pn, pseq >= 0, c='r', marker='x',
                                      markersize=100, ax=ax)
        op.topotools.plot_coordinates(pn, pseq < 0, c='c', marker='o',
                                      markersize=100, ax=ax)

    # %%
    if 1:
        drn = op.algorithms.Drainage(network=pn, phase=water)
        drn.set_inlet_BC(pn.pores('left'))
        pressures = np.unique(ip['pore.invasion_pressure'])
        # pressures = np.logspace(np.log10(0.1e3), np.log10(2e4), 100)
        drn.run(pressures=pressures)
        drn.set_outlet_BC(pn.pores('right'))
        # drn.apply_trapping()

        fig, ax = plt.subplots(1, 1)
        ax.step(*ip.pc_curve(), 'b', where='post')
        ax.step(*drn.pc_curve(), 'r--', where='post')
        ax.set_ylim([0, 1])

    # %%
    if 0:
        pseq = ip['pore.invasion_sequence']
        pseq[pseq == -1] = pseq.max() + 1
        tseq = ip['throat.invasion_sequence']
        tseq[tseq == -1] = tseq.max() + 1
        for j, i in enumerate(tqdm(np.unique(tseq)[:-1])):
            ax = op.topotools.plot_connections(pn, tseq <= i, linestyle='--',
                                               c='r', linewidth=3)
            op.topotools.plot_coordinates(pn, pseq <= i, c='b', marker='x',
                                          markersize=100, ax=ax)
            op.topotools.plot_coordinates(pn, ip['pore.trapped'], c='c',
                                          marker='o', markersize=100, ax=ax)
            plt.savefig(f"{str(j).zfill(3)}.png")
            plt.close()
