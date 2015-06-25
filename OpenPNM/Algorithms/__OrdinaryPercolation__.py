# -*- coding: utf-8 -*-
"""
===============================================================================
module __OrdinaryPercolation__: Ordinary Percolation Algorithm
===============================================================================

"""

import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from OpenPNM.Algorithms import GenericAlgorithm
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)


class OrdinaryPercolation(GenericAlgorithm):
    r"""
    Simulates a capillary drainage experiment by applying a list of increasing
    capillary pressures.

    Parameters
    ----------
    network : OpenPNM Network Object
        The network upon which the simulation will be run

    name : string, optional
        The name to assign to the Algorithm Object

    """

    def setup(self,
              inv_phase,
              p_inlets=[],
              p_outlets=[],
              t_entry='throat.capillary_pressure'):
        r"""
        inv_phase : OpenPNM Phase Object
            The invading phase to be injected into the Network

        p_inlets : array_like
            The injection points from which the invading phase accesses the
            Network.  If no inlets are specified then the algorithm assumes
            no access limitations apply to the invading phase, which is
            equivalent to performaing a standard bond ordinary percolation.


        Notes
        -----
        The 'inlet' pores are initially filled with invading fluid to start the
        simulation.  To avoid the capillary pressure curve showing a non-zero
        starting saturation at low pressures, it is necessary to apply boundary
        pores that have zero-volume, and set these as the inlets.
        """
        self['throat.entry_pressure'] = inv_phase[t_entry]
        self['pore.invaded'] = sp.inf
        self['throat.invaded'] = sp.inf
        self._inv_phase = inv_phase

    def set_inlets(self, pores):
        r"""
        Specify inlet locations

        Parameters
        ----------
        p_inlets : array_like
            The injection points from which the invading phase accesses the
            Network.  If no inlets are specified then the algorithm assumes
            no access limitations apply to the invading phase, which is
            equivalent to performaing a standard bond ordinary percolation.


        Notes
        -----
        The 'inlet' pores are initially filled with invading fluid to start the
        simulation.  To avoid the capillary pressure curve showing a non-zero
        starting saturation at low pressures, it is necessary to apply boundary
        pores that have zero-volume, and set these as the inlets.
        """
        Ps = sp.array(pores)
        if sp.size(Ps) > 0:
            if Ps.dtype == bool:
                Ps = self._net.Ps[Ps]
            self['pore.inlets'] = False
            self['pore.inlets'][Ps] = True

    def run(self, npts=25, inv_points=None, access_limited=True, **kwargs):
        r"""
        Parameters
        ----------
        npts : int (default = 25)
            The number of pressure points to apply.  The list of pressures
            is logarithmically spaced between the lowest and highest throat
            entry pressures in the network.

        inv_points : array_like, optional
            A list of specific pressure point(s) to apply.

        """
        self._AL = access_limited
        if inv_points is None:
            logger.info('Generating list of invasion pressures')
            min_p = sp.amin(self['throat.entry_pressure']) * 0.98  # nudge down
            max_p = sp.amax(self['throat.entry_pressure']) * 1.02  # bump up
            inv_points = sp.logspace(sp.log10(min_p),
                                     sp.log10(max_p),
                                     npts)

        self._npts = sp.size(inv_points)
        # Execute calculation
        self._do_outer_iteration_stage(inv_points)

    def _do_outer_iteration_stage(self, inv_points):
        # Generate curve from points
        for inv_val in inv_points:
            # Apply one applied pressure and determine invaded pores
            logger.info('Applying capillary pressure: ' + str(inv_val))
            self._do_one_inner_iteration(inv_val)

        # Find invasion sequence values (to correspond with IP algorithm)
        self['pore.inv_seq'] = sp.searchsorted(sp.unique(self['pore.invaded']),
                                               self['pore.invaded'])
        self['throat.inv_seq'] = sp.searchsorted(sp.unique(self['throat.invaded']),
                                                 self['throat.invaded'])

    def _do_one_inner_iteration(self, inv_val):
        r"""
        Determine which throats are invaded at a given applied capillary
        pressure.

        """
        # Generate a tlist containing boolean values for throat state
        Tinvaded = self['throat.entry_pressure'] <= inv_val
        # Find all pores that can be invaded at specified pressure
        [pclusters, tclusters] = self._net.find_clusters2(mask=Tinvaded,
                                                          t_labels=True)
        if self._AL:
            # Identify clusters connected to invasion sites
            inv_clusters = sp.unique(pclusters[self['pore.inlets']])
            inv_clusters = inv_clusters[inv_clusters >= 0]
        else:
            # All clusters are invasion sites
            inv_clusters = pclusters
        # Find pores on the invading clusters
        pmask = np.in1d(pclusters, inv_clusters)
        # Store current applied pressure in newly invaded pores
        inds = (self['pore.invaded'] == sp.inf) * (pmask)
        self['pore.invaded'][inds] = inv_val
        # Find throats on the invading clusters
        tmask = np.in1d(tclusters, inv_clusters)
        # Store current applied pressure in newly invaded throats
        inds = (self['throat.invaded'] == sp.inf) * (tmask)
        self['throat.invaded'][inds] = inv_val

    def evaluate_trapping(self, p_outlets):
        r"""
        Finds trapped pores and throats after a full ordinary
        percolation drainage has been run

        Parameters
        ----------
        p_outlets : array_like
            A list of pores that define the wetting phase outlets.
            Disconnection from these outlets results in trapping.

        """
        self['pore.trapped'] = sp.zeros([self.Np, ], dtype=float)
        self['throat.trapped'] = sp.zeros([self.Nt, ], dtype=float)
        try:
            # Get points used in OP
            inv_points = sp.unique(self['pore.invaded'])
        except:
            raise Exception('Orindary percolation has not been run!')
        tind = self._net.throats()
        conns = self._net.find_connected_pores(tind)
        for inv_val in inv_points[0:-1]:
            # Find clusters of defender pores
            Pinvaded = self['pore.invaded'] <= inv_val
            Cstate = sp.sum(Pinvaded[conns], axis=1)
            Tinvaded = self['throat.invaded'] <= inv_val
            # 0 = all open, 1=1 pore filled,
            # 2=2 pores filled 3=2 pores + 1 throat filled
            Cstate = Cstate + Tinvaded
            clusters = self._net.find_clusters(Cstate == 0)
            # Clean up clusters (invaded = -1, defended >=0)
            clusters = clusters * (~Pinvaded) - (Pinvaded)
            # Identify clusters connected to outlet sites
            out_clusters = sp.unique(clusters[p_outlets])
            trapped_pores = ~sp.in1d(clusters, out_clusters)
            trapped_pores[Pinvaded] = False
            if sum(trapped_pores) > 0:
                inds = (self['pore.trapped'] == 0) * trapped_pores
                self['pore.trapped'][inds] = inv_val
                trapped_throats = self._net.find_neighbor_throats(trapped_pores)
                trapped_throat_array = np.asarray([False] * len(Cstate))
                trapped_throat_array[trapped_throats] = True
                inds = (self['throat.trapped'] == 0) * trapped_throat_array
                self['throat.trapped'][inds] = inv_val
                inds = (self['throat.trapped'] == 0) * (Cstate == 2)
                self['throat.trapped'][inds] = inv_val
        self['pore.trapped'][self['pore.trapped'] > 0] = sp.inf
        self['throat.trapped'][self['throat.trapped'] > 0] = sp.inf

    def return_results(self, Pc=0, seq=None, sat=None, occupancy='occupancy'):
        r"""
        Updates the occupancy status of invading and defending phases
        as determined by the OP algorithm

        """
        p_inv = self['pore.invaded']
        self._inv_phase['pore.invaded'] = p_inv
        t_inv = self['throat.invaded']
        self._inv_phase['throat.invaded'] = t_inv
        # Apply invasion sequence values (to correspond with IP algorithm)
        p_seq = self['pore.inv_seq']
        self._inv_phase['pore.inv_seq'] = p_seq
        t_seq = self['throat.inv_seq']
        self._inv_phase['throat.inv_seq'] = t_seq
        # Apply saturation to pores and throats
        self._inv_phase['pore.inv_sat'] = self['pore.inv_sat']
        self._inv_phase['throat.inv_sat'] = self['throat.inv_sat']

        if sat is not None:
            p_inv = self['pore.inv_sat'] <= sat
            t_inv = self['throat.inv_sat'] <= sat
            # Apply occupancy to invading phase
            temp = sp.array(p_inv, dtype=sp.float_, ndmin=1)
            self._inv_phase['pore.' + occupancy] = temp
            temp = sp.array(t_inv, dtype=sp.float_, ndmin=1)
            self._inv_phase['throat.' + occupancy] = temp
            # Apply occupancy to defending phase
            if self._phase_def is not None:
                temp = sp.array(~p_inv, dtype=sp.float_, ndmin=1)
                self._phase_def['pore.' + occupancy] = temp
                temp = sp.array(~t_inv, dtype=sp.float_, ndmin=1)
                self._phase_def['throat.' + occupancy] = temp
        elif seq is not None:
            p_seq = self['pore.inv_seq'] <= seq
            t_seq = self['throat.inv_seq'] <= seq
            # Apply occupancy to invading phase
            temp = sp.array(p_seq, dtype=sp.float_, ndmin=1)
            self._inv_phase['pore.' + occupancy] = temp
            temp = sp.array(t_seq, dtype=sp.float_, ndmin=1)
            self._inv_phase['throat.' + occupancy] = temp
            # Apply occupancy to defending phase
            if self._phase_def is not None:
                temp = sp.array(~p_seq, dtype=sp.float_, ndmin=1)
                self._phase_def['pore.' + occupancy] = temp
                temp = sp.array(~t_seq, dtype=sp.float_, ndmin=1)
                self._phase_def['throat.' + occupancy] = temp
        else:
            p_inv = self['pore.inv_Pc'] <= Pc
            t_inv = self['throat.inv_Pc'] <= Pc
            # Apply occupancy to invading phase
            temp = sp.array(p_inv, dtype=sp.float_, ndmin=1)
            self._inv_phase['pore.' + occupancy] = temp
            temp = sp.array(t_inv, dtype=sp.float_, ndmin=1)
            self._inv_phase['throat.' + occupancy] = temp
            # Apply occupancy to defending phase
            if self._phase_def is not None:
                temp = sp.array(~p_inv, dtype=sp.float_, ndmin=1)
                self._phase_def['pore.' + occupancy] = temp
                temp = sp.array(~t_inv, dtype=sp.float_, ndmin=1)
                self._phase_def['throat.' + occupancy] = temp

    def plot_drainage_curve(self, pore_volume='volume', throat_volume='volume',
                            pore_label='all', throat_label='all'):
        r"""
        Plot drainage capillary pressure curve
        """
        try:
            PcPoints = sp.unique(self['pore.invaded'])
        except:
            raise Exception('Cannot print drainage curve: ordinary percolation \
                             simulation has not been run')
        pores = self._net.pores(labels=pore_label)
        throats = self._net.throats(labels=throat_label)
        Snwp_t = sp.zeros_like(PcPoints)
        Snwp_p = sp.zeros_like(PcPoints)
        Snwp_all = sp.zeros_like(PcPoints)
        Pvol = self._net['pore.' + pore_volume]
        Tvol = self._net['throat.' + throat_volume]
        Pvol_tot = sp.sum(Pvol)
        Tvol_tot = sp.sum(Tvol)
        vol_tot = Pvol_tot + Tvol_tot
        for i in range(0, sp.size(PcPoints)):
            Pc = PcPoints[i]
            Snwp_p[i] = sp.sum(Pvol[self['pore.invaded'][pores] <= Pc]) / vol_tot
            Snwp_t[i] = sp.sum(Tvol[self['throat.invaded'][throats] <= Pc]) / vol_tot
            Snwp_all[i] = (sp.sum(Tvol[self['throat.invaded'][throats] <= Pc]) +
                           sp.sum(Pvol[self['pore.invaded'][pores] <= Pc])) / vol_tot
        if sp.mean(self._inv_phase['pore.contact_angle']) < 90:
            Snwp_p = 1 - Snwp_p
            Snwp_t = 1 - Snwp_t
            PcPoints *= -1
        plt.plot(PcPoints, Snwp_all, 'g.-')
        plt.plot(PcPoints, Snwp_p, 'r.-')
        plt.plot(PcPoints, Snwp_t, 'b.-')
        r"""
        TODO: Add legend to distinguish the pore and throat curves
        """
        plt.show()

    def plot_primary_drainage_curve(self, pore_volume='volume',
                                    throat_volume='volume', pore_label='all',
                                    throat_label='all'):
        r"""
        Plot the primary drainage curve as the capillary pressure on ordinate
        and total saturation of the wetting phase on the abscissa.
        This is the preffered style in the petroleum engineering
        """
        try:
            PcPoints = sp.unique(self['pore.invaded'])
        except:
            raise Exception('Cannot print drainage curve: ordinary percolation \
                            simulation has not been run')
        pores = self._net.pores(labels=pore_label)
        throats = self._net.throats(labels=throat_label)
        p_inv = self['pore.invaded']
        t_inv = self['throat.invaded']
        Snwp_t = sp.zeros_like(PcPoints)
        Snwp_p = sp.zeros_like(PcPoints)
        Snwp_all = sp.zeros_like(PcPoints)
        Swp_all = sp.zeros_like(PcPoints)
        Pvol = self._net['pore.' + pore_volume]
        Tvol = self._net['throat.' + throat_volume]
        Pvol_tot = sp.sum(Pvol)
        Tvol_tot = sp.sum(Tvol)
        for i in range(0, sp.size(PcPoints)):
            Pc = PcPoints[i]
            Snwp_p[i] = sp.sum(Pvol[p_inv[pores] <= Pc]) / Pvol_tot
            Snwp_t[i] = sp.sum(Tvol[t_inv[throats] <= Pc]) / Tvol_tot
            Snwp_all[i] = (sp.sum(Tvol[t_inv[throats] <= Pc]) +
                           sp.sum(Pvol[p_inv[pores] <= Pc])) / \
                          (Tvol_tot + Pvol_tot)
            Swp_all[i] = 1 - Snwp_all[i]
        plt.plot(Swp_all, PcPoints, 'k.-')
        plt.xlim(xmin=0)
        plt.xlabel('Saturation of wetting phase')
        plt.ylabel('Capillary Pressure [Pa]')
        plt.title('Primay Drainage Curve')
        plt.grid(True)
        plt.show()
