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


class OrdinaryPercolation2(GenericAlgorithm):
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
              p_residual=[],
              t_entry='throat.capillary_pressure'):
        r"""
        inv_phase : OpenPNM Phase Object
            The invading phase to be injected into the Network

        p_inlets : array_like
            The injection points from which the invading phase accesses the
            Network.  If no inlets are specified then the algorithm assumes
            no access limitations apply to the invading phase, which is
            equivalent to performaing a standard bond ordinary percolation.

        p_outlets : array_like
            The points through which the defending phase exits the Network.  If
            outlets are given the algorithm calculates the trapping of the
            defending phase.  If no outlets are specified, then the algorithm
            assumes that no trapping occurs.

        Notes
        -----
        The 'inlet' pores are initially filled with invading fluid to start the
        simulation.  To avoid the capillary pressure curve showing a non-zero
        starting saturation at low pressures, it is necessary to apply boundary
        pores that have zero-volume, and set these as the inlets.
        """
        p_inlets = sp.array(p_inlets)
        if sp.size(p_inlets) > 0:
            if p_inlets.dtype == bool:
                p_inlets = self._net.Ps[p_inlets]
            self['pore.inlets'] = False
            self['pore.inlets'][p_inlets] = True

        p_outlets = sp.array(p_outlets)
        if sp.size(p_outlets) > 0:
            if p_outlets.dtype == bool:
                p_outlets = self._net.Ps[p_outlets]
            self['pore.outlets'] = p_outlets

        self['throat.entry_pressure'] = inv_phase[t_entry]

    def run(self, npts=25, inv_points=None, access_limited=True):
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
        # Execute calculation
        self._do_outer_iteration_stage(inv_points)

    def _do_outer_iteration_stage(self, inv_points):
        # Generate curve from points
        for inv_val in inv_points:
            # Apply one applied pressure and determine invaded pores
            logger.info('Applying capillary pressure: ' + str(inv_val))
            self._do_one_inner_iteration(inv_val)
        # Find invasion sequence values (to correspond with IP algorithm)
        self._p_seq = sp.searchsorted(sp.unique(self._p_inv), self._p_inv)
        self._t_seq = sp.searchsorted(sp.unique(self._t_inv), self._t_inv)
        self['pore.inv_seq'] = self._p_seq
        self['throat.inv_seq'] = self._t_seq
        # Calculate Saturations
        v_total = sp.sum(self._net['pore.volume']) + \
            sp.sum(self._net['throat.volume'])
        sat = 0.
        self['pore.inv_sat'] = 1.
        self['throat.inv_sat'] = 1.
        for i in range(self._npts):
            inv_pores = sp.where(self._p_seq == i)[0]
            inv_throats = sp.where(self._t_seq == i)[0]
            new_sat = (sum(self._net['pore.volume'][inv_pores]) +
                       sum(self._net['throat.volume'][inv_throats])) / v_total
            sat += new_sat
            self['pore.inv_sat'][inv_pores] = sat
            self['throat.inv_sat'][inv_throats] = sat

    def _do_one_inner_iteration(self, inv_val):
        r"""
        Determine which throats are invaded at a given applied capillary
        pressure.

        """
        # Generate a tlist containing boolean values for throat state
        Tinvaded = self['throat.entry_pressure'] <= inv_val
        # Find all pores that can be invaded at specified pressure
        clusters = self._net.find_clusters2(Tinvaded)
        temp = sp.unique(clusters)
        compressed_labels = {temp[i]: i for i in range(0, len(temp))}
        clusters = [compressed_labels[clusters[i]]
                    for i in range(0, len(clusters))]
        clusters = sp.array(clusters)
        if self._AL:
            # Identify clusters connected to invasion sites
            inv_clusters = sp.unique(clusters[self['pore.inlets']])
        else:
            # All clusters are invasion sites
            inv_clusters = clusters
        # Store invasion pressure in pores and throats
        pmask = np.in1d(clusters, inv_clusters)
        # Store result of invasion step
        inds = (self['pore.invaded'] == sp.inf) * (pmask)
        self['pore.invaded'][inds] = inv_val
        # Determine Pc_invaded for throats as well
        temp = self._net['throat.conns'].copy()
        tmask = (pmask[temp[:, 0]] + pmask[temp[:, 1]]) * Tinvaded
        self._t_inv[(self._t_inv == sp.inf) * (tmask)] = inv_val

    def evaluate_trapping(self, outlets):
        r"""
        Finds trapped pores and throats after a full ordinary
        percolation drainage has been run

        Parameters
        ----------
        outlets : array_like
            A list of pores that define the wetting phase outlets.
            Disconnection from these outlets results in trapping.

        """
        self._p_trap = sp.zeros_like(self._p_inv, dtype=float)
        self._t_trap = sp.zeros_like(self._t_inv, dtype=float)
        try:
            inv_points = sp.unique(self._p_inv)  # Get points used in OP
        except:
            logger.error('Orindary percolation has not been run!')
            raise Exception('Aborting algorithm')
        tind = self._net.throats()
        conns = self._net.find_connected_pores(tind)
        for inv_val in inv_points[0:-1]:
            # Find clusters of defender pores
            Pinvaded = self._p_inv <= inv_val
            Cstate = sp.sum(Pinvaded[conns], axis=1)
            Tinvaded = self._t_inv <= inv_val
            # 0 = all open, 1=1 pore filled,
            # 2=2 pores filled 3=2 pores + 1 throat filled
            Cstate = Cstate + Tinvaded
            clusters = self._net.find_clusters(Cstate == 0)
            # Clean up clusters (invaded = -1, defended >=0)
            clusters = clusters * (~Pinvaded) - (Pinvaded)
            # Identify clusters connected to outlet sites
            out_clusters = sp.unique(clusters[outlets])
            trapped_pores = ~sp.in1d(clusters, out_clusters)
            trapped_pores[Pinvaded] = False
            if sum(trapped_pores) > 0:
                self._p_trap[(self._p_trap == 0) * trapped_pores] = inv_val
                trapped_throats = self._net.find_neighbor_throats(trapped_pores)
                trapped_throat_array = np.asarray([False] * len(Cstate))
                trapped_throat_array[trapped_throats] = True
                self._t_trap[(self._t_trap == 0) * trapped_throat_array] = inv_val
                self._t_trap[(self._t_trap == 0) * (Cstate == 2)] = inv_val
        self._p_inv[self._p_trap > 0] = sp.inf
        self._t_inv[self._t_trap > 0] = sp.inf
        self['pore.inv_Pc'] = self._p_inv
        self['throat.inv_Pc'] = self._t_inv

    def return_results(self, Pc=0, seq=None, sat=None, occupancy='occupancy'):
        r"""
        Updates the occupancy status of invading and defending phases
        as determined by the OP algorithm

        """
        p_inv = self['pore.inv_Pc']
        self._phase_inv['pore.inv_Pc'] = p_inv
        t_inv = self['throat.inv_Pc']
        self._phase_inv['throat.inv_Pc'] = t_inv
        # Apply invasion sequence values (to correspond with IP algorithm)
        p_seq = self['pore.inv_seq']
        self._phase_inv['pore.inv_seq'] = p_seq
        t_seq = self['throat.inv_seq']
        self._phase_inv['throat.inv_seq'] = t_seq
        # Apply saturation to pores and throats
        self._phase_inv['pore.inv_sat'] = self['pore.inv_sat']
        self._phase_inv['throat.inv_sat'] = self['throat.inv_sat']

        if sat is not None:
            p_inv = self['pore.inv_sat'] <= sat
            t_inv = self['throat.inv_sat'] <= sat
            # Apply occupancy to invading phase
            temp = sp.array(p_inv, dtype=sp.float_, ndmin=1)
            self._phase_inv['pore.' + occupancy] = temp
            temp = sp.array(t_inv, dtype=sp.float_, ndmin=1)
            self._phase_inv['throat.' + occupancy] = temp
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
            self._phase_inv['pore.' + occupancy] = temp
            temp = sp.array(t_seq, dtype=sp.float_, ndmin=1)
            self._phase_inv['throat.' + occupancy] = temp
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
            self._phase_inv['pore.' + occupancy] = temp
            temp = sp.array(t_inv, dtype=sp.float_, ndmin=1)
            self._phase_inv['throat.' + occupancy] = temp
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
            PcPoints = sp.unique(self['pore.inv_Pc'])
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
            Snwp_p[i] = sp.sum(Pvol[self._p_inv[pores] <= Pc]) / vol_tot
            Snwp_t[i] = sp.sum(Tvol[self._t_inv[throats] <= Pc]) / vol_tot
            Snwp_all[i] = (sp.sum(Tvol[self._t_inv[throats] <= Pc]) +
                           sp.sum(Pvol[self._p_inv[pores] <= Pc])) / vol_tot
        if sp.mean(self._phase_inv['pore.contact_angle']) < 90:
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
            PcPoints = sp.unique(self['pore.inv_Pc'])
        except:
            raise Exception('Cannot print drainage curve: ordinary percolation \
                            simulation has not been run')
        pores = self._net.pores(labels=pore_label)
        throats = self._net.throats(labels=throat_label)
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
            Snwp_p[i] = sp.sum(Pvol[self._p_inv[pores] <= Pc]) / Pvol_tot
            Snwp_t[i] = sp.sum(Tvol[self._t_inv[throats] <= Pc]) / Tvol_tot
            Snwp_all[i] = (sp.sum(Tvol[self._t_inv[throats] <= Pc]) +
                           sp.sum(Pvol[self._p_inv[pores] <= Pc])) / \
                          (Tvol_tot + Pvol_tot)
            Swp_all[i] = 1 - Snwp_all[i]
        plt.plot(Swp_all, PcPoints, 'k.-')
        plt.xlim(xmin=0)
        plt.xlabel('Saturation of wetting phase')
        plt.ylabel('Capillary Pressure [Pa]')
        plt.title('Primay Drainage Curve')
        plt.grid(True)
        plt.show()
