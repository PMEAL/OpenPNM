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
    Simulates a capillary drainage experiment by looping through a list of
    capillary pressures.

    Parameters
    ----------
    network : OpenPNM Network Object
        The network upon which the simulation will be run

    invading_phase : OpenPNM Phase Object
        The phase to be forced into the network at increasingly high pressures

    defending_phase : OpenPNM Phase Object, optional
        The phase originally residing in the network prior to invasion.  This
        is only necessary so that the pressure at which the phase is drained
        can be stored on the phase.

    name : string, optional
        The name to assign to the Algorithm Object

    Examples
    --------
    >>> import OpenPNM
    >>> pn = OpenPNM.Network.TestNet()
    >>> geo = OpenPNM.Geometry.TestGeometry(network=pn, pores=pn.pores(),
    ...                                     throats=pn.throats())
    >>> phase1 = OpenPNM.Phases.TestPhase(network=pn)
    >>> phase2 = OpenPNM.Phases.TestPhase(network=pn)
    >>> phys1 = OpenPNM.Physics.TestPhysics(network=pn,
    ...                                     phase=phase1,
    ...                                     pores=pn.pores(),
    ...                                     throats=pn.throats())
    >>> phys2 = OpenPNM.Physics.TestPhysics(network=pn, phase=phase2,
    ...                                     pores=pn.pores(), throats=pn.throats())
    >>> OP = OpenPNM.Algorithms.OrdinaryPercolation(network=pn,
    ...                                             invading_phase=phase1,
    ...                                             defending_phase=phase2)
    >>> OP.run(inlets=pn.pores('top'))
    >>> med_Pc = sp.median(OP['pore.inv_Pc'])
    >>> OP.return_results(med_Pc)
    >>> print(len(phase1.pores('occupancy')))
    71

    To run this algorithm, use 'setup()' to provide the necessary simulation
    """

    def __init__(self, invading_phase=None, defending_phase=None,
                 residual_pores=None, residual_throats=None, **kwargs):
        super().__init__(**kwargs)
        self._phase_inv = invading_phase
        self._phase_def = defending_phase
        self._residual_pores = residual_pores
        self._residual_throats = residual_throats

        logger.debug('Create Drainage Percolation Algorithm Object')

    def run(self, inlets, outlets=[], npts=25, inv_points=None,
            capillary_pressure='capillary_pressure', access_limited=True,
            trapping=False, **kwargs):
        r"""
        Parameters
        ----------
        inlets : array_like
            The list of pores which are the injection sources

        npts : int, optional
            The number of pressure points to apply.  The list of pressures
            is logarithmically spaced between the lowest and highest throat
            entry pressures in the network.

        inv_points : array_like, optional
            A list of specific pressure points to apply.

        access_limited : boolean
            Only pores and throats connected to the inlet sites can be invaded

        trapping : boolean
            Wetting phase that is cut-off from the outlets becomes immobile.
            If outlet pores have not been provided then this argument is
            ignored.

        Notes
        -----
        The 'inlet' pores are initially filled with invading fluid to start the
        simulation.  To avoid the capillary pressure curve showing a non-zero
        starting saturation at low pressures, it is necessary to apply boundary
        pores that have 0 volume, and set these as the inlets.

        """
        self._inv_sites = inlets
        self._out_sites = outlets
        self._npts = npts
        self._p_cap = capillary_pressure  # Name of throat entry pressure prop
        self._AL = access_limited
        self._TR = trapping

        # Create pore and throat conditions lists to store inv_val at
        # which each is invaded
        self._p_inv = sp.zeros((self._net.num_pores(),), dtype=float)
        self._p_inv.fill(sp.inf)
        self._p_seq = sp.zeros_like(self._p_inv, dtype=int)
        self._t_inv = sp.zeros((self._net.num_throats(),), dtype=float)
        self._t_inv.fill(sp.inf)
        self._t_seq = sp.zeros_like(self._t_inv, dtype=int)
        # Determine the invasion pressures to apply
        try:
            self._t_cap = self._phase_inv['throat.' + self._p_cap]
        except:
            logger.error('Capillary pressure not assigned to invading phase ' +
                         self._phase_inv.name +
                         ', check for capillary pressure in defending phase ' +
                         self._phase_def.name + ' instead')
            try:
                self._t_cap = self._phase_def['throat.' + self._p_cap]
            except:
                pass
                logger.error('Capillary pressure neither assigned to defending \
                              phase ' + self._phase_def.name +
                             ' nor to invading phase ' + self._phase_inv.name)
        if inv_points is None:
            min_p = sp.amin(self._t_cap) * 0.98  # nudge min_p down slightly
            max_p = sp.amax(self._t_cap) * 1.02  # bump max_p up slightly
            logger.info('Generating list of invasion pressures')
            if min_p == 0:
                min_p = sp.linspace(min_p, max_p, self._npts)[1]
            self._inv_points = sp.logspace(sp.log10(min_p),
                                           sp.log10(max_p),
                                           self._npts)
        else:
            self._inv_points = inv_points
        self._do_outer_iteration_stage()

    def _do_outer_iteration_stage(self):
        # Generate curve from points
        for inv_val in self._inv_points:
            # Apply one applied pressure and determine invaded pores
            logger.info('Applying capillary pressure: ' + str(inv_val))
            self._do_one_inner_iteration(inv_val)
        # Store results using networks' get/set method
        self['pore.inv_Pc'] = self._p_inv
        self['throat.inv_Pc'] = self._t_inv
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
        if self._TR:
            logger.info('Evaluating trapping')
            self.evaluate_trapping(outlets=self._out_sites)

    def _do_one_inner_iteration(self, inv_val):
        r"""
        Determine which throats are invaded at a given applied capillary pressure.

        This function uses the scipy.csgraph module for the cluster labeling
        algorithm (connected_components).

        """
        # Generate a tlist containing boolean values for throat state
        Tinvaded = self._t_cap <= inv_val
        # Finding all pores that can be invaded at specified pressure
        clusters = self._net.find_clusters(Tinvaded)
        # Find all pores with at least 1 invaded throat (invaded)
        Pinvaded = sp.zeros_like(clusters, dtype=bool)
        Ts = self._net.throats()
        P12 = self._net.find_connected_pores(Ts)
        temp = P12[Tinvaded]
        temp = sp.hstack((temp[:, 0], temp[:, 1]))
        Pinvaded[temp] = True
        if self._AL:
            # Add injection sites to Pinvaded
            Pinvaded[self._inv_sites] = True
            # Clean up clusters (not invaded = -1, invaded >=0)
            clusters = clusters * (Pinvaded) - (~Pinvaded)
            # Identify clusters connected to invasion sites
            inv_clusters = sp.unique(clusters[self._inv_sites])
        else:
            # Clean up clusters (not invaded = -1, invaded >=0)
            clusters = clusters * (Pinvaded) - (~Pinvaded)
            # All clusters are invasion sites
            inv_clusters = sp.r_[0:self._net.num_pores()]
        # Store invasion pressure in pores and throats
        pmask = np.in1d(clusters, inv_clusters)
        # Store result of invasion step
        self._p_inv[(self._p_inv == sp.inf) * (pmask)] = inv_val
        # Determine Pc_invaded for throats as well
        temp = self._net['throat.conns']
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
        Pvol = self._net['pore.' + pore_volume]
        Tvol = self._net['throat.' + throat_volume]
        Pvol_tot = sum(Pvol)
        Tvol_tot = sum(Tvol)
        for i in range(0, sp.size(PcPoints)):
            Pc = PcPoints[i]
            Snwp_p[i] = sum(Pvol[self._p_inv[pores] <= Pc]) / Pvol_tot
            Snwp_t[i] = sum(Tvol[self._t_inv[throats] <= Pc]) / Tvol_tot
        if sp.mean(self._phase_inv['pore.contact_angle']) < 90:
            Snwp_p = 1 - Snwp_p
            Snwp_t = 1 - Snwp_t
            PcPoints *= -1
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
        Pvol_tot = sum(Pvol)
        Tvol_tot = sum(Tvol)
        for i in range(0, sp.size(PcPoints)):
            Pc = PcPoints[i]
            Snwp_p[i] = sum(Pvol[self._p_inv[pores] <= Pc]) / Pvol_tot
            Snwp_t[i] = sum(Tvol[self._t_inv[throats] <= Pc]) / Tvol_tot
            Snwp_all[i] = (sum(Tvol[self._t_inv[throats] <= Pc]) +
                           sum(Pvol[self._p_inv[pores] <= Pc])) / \
                          (Tvol_tot + Pvol_tot)
            Swp_all[i] = 1 - Snwp_all[i]
        plt.plot(Swp_all, PcPoints, 'k.-')
        plt.xlim(xmin=0)
        plt.xlabel('Saturation of wetting phase')
        plt.ylabel('Capillary Pressure [Pa]')
        plt.title('Primay Drainage Curve')
        plt.grid(True)
        plt.show()
