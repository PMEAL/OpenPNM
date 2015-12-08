# -*- coding: utf-8 -*-
"""
===============================================================================
module __Drainage__: Capillary Pressure Curve Simulation
===============================================================================

"""

import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from OpenPNM.Algorithms import GenericAlgorithm
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)


class Drainage(GenericAlgorithm):
    r"""
    Simulates a capillary drainage experiment by applying a list of increasing
    capillary pressures and invading all throat that are accessible and
    invadable at the given pressure

    Parameters
    ----------
    network : OpenPNM Network Object
        The network upon which the simulation will be run

    phase : OpenPNM Phase Object
        The phase which will

    name : string, optional
        The name to assign to the Algorithm Object

    throat_prop : string
        The Phase property containing the values to use for the invasion
        calculation.  The default is 'throat.capillary_pressure'.

    trapping : boolean
        Specifies whether trapping of the defending phase should be included.
        The default is False.  If you wish to apply trapping you can also
        specify outlet pores from which the defender can escape using
        ``set_outlets``. If no outlets are specified the it is assume the
        defender can escape anywhere that is not already filled by invading
        fluid.

    Notes
    -----

    """

    def __init__(self,
                 network,
                 name=None):
        super().__init__(network=network, name=name)

    def setup(self,
              invading_phase,
              throat_prop='throat.capillary_pressure',
              trapping=False):
        r"""
        This method is used to specify necessary arguments to the simulation.
        Note that all of these arguments can be sent to the initialization
        since ``setup`` is automatically called.  This method useful for
        resetting the Algorithm or apply more explicit control.

        Parameters
        ----------
        invading_phase : OpenPNM Phase object
            The Phase object containing the physical properties of the invading
            fluid.
        throat_prop : string
            The dictionary key on the Phase object where the throat entry
            pressure values can be found.  The default is
            'throat.capillary_pressure'.
        trapping : boolean
            Specifies whether wetting phase trapping should be included or not.
            The default is False.  Note that wetting phase outlets can be
            specified using the ``set_outlets`` method.  Otherwise it is
            assumed the wetting phase has no outlets.

        """
        self['throat.entry_pressure'] = invading_phase[throat_prop]
        self['pore.inv_Pc'] = sp.inf
        self['throat.inv_Pc'] = sp.inf
        self['pore.inlets'] = False
        self['pore.outlets'] = False
        self._inv_phase = invading_phase
        self._trapping = False

    def set_inlets(self, pores):
        r"""
        Set the locations from which the invading phase enters the network.

        Parameters
        ----------
        pores : array_like
            An array of pore numbers that are initially filled with invading
            phase, and from which clusters of invading phase grow and invade
            into the network.

        Notes
        -----
        The 'inlet' pores are initially filled with invading fluid to start the
        simulation.  To avoid the capillary pressure curve showing a non-zero
        starting saturation at low pressures, it is necessary to apply boundary
        pores that have zero-volume, and set these as the inlets.
        """
        Ps = self._parse_locations(pores)
        if 'pore.outlets' in self.keys():
            if sum(self['pore.outlets'][Ps]) > 0:
                raise Exception('Some inlets are already defined as outlets')
        self['pore.inlets'][Ps] = True

    def set_outlets(self, pores):
        r"""
        Set the locations through which defending phase exits the network.

        Parameters
        ----------
        pores : array_like
            An array of pore numbers where defending phase can exit.  Any
            defending phase that does not have access to these pores will be
            trapped.
        """
        Ps = self._parse_locations(pores)
        if 'pore.outlets' in self.keys():
            if sum(self['pore.inlets'][Ps]) > 0:
                raise Exception('Some outlets are already defined as inlets')
        self['pore.outlets'][Ps] = True

    def set_residual(self, pores=None, throats=None):
        r"""
        Specify locations of the residual wetting phase
        """
        if pores is not None:
            Ps = self._parse_locations(pores)
            self['pore.residual'] = False
            self['pore.residual'][Ps] = True
        if throats is not None:
            Ts = self._parse_locations(throats)
            self['throat.residual'] = False
            self['throat.residual'][Ts] = True

    def run(self, npts=25, inv_pressures=None):
        r"""
        Run the algorithm for specified number of points or at given capillary
        pressures.

        Parameters
        ----------
        npts : scalar
            The number of points to obtain on the curve.  The points are
            automatically selected to span the range of capillary pressures
            using a logarithmic spacing (more points are lower capillary
            pressure values).
        inv_pressures : array_like
            A list of capillary pressures to apply. List should contain
            increasing and unique values.
        """
        # If no invasion points are given then generate some
        if inv_pressures is None:
            logger.info('Generating list of invasion pressures')
            min_p = sp.amin(self['throat.entry_pressure']) * 0.98  # nudge down
            max_p = sp.amax(self['throat.entry_pressure']) * 1.02  # bump up
            inv_points = sp.logspace(sp.log10(min_p),
                                     sp.log10(max_p),
                                     npts)
        else:
            # Make sure the given invastion points are sensible
            inv_points = sp.unique(inv_pressures)

        # Ensure inlets are set
        if sp.sum(self['pore.inlets']) == 0:
            raise Exception('Inlet pores have not been specified')

        # Execute calculation
        self._do_outer_iteration_stage(inv_points)

    def _do_outer_iteration_stage(self, inv_points):
        r"""
        """
        # Generate curve from points
        for inv_val in inv_points:
            # Apply one applied pressure and determine invaded pores
            logger.info('Applying capillary pressure: ' + str(inv_val))
            self._do_one_inner_iteration(inv_val)

        # Find invasion sequence values (to correspond with IP algorithm)
        Pinv = self['pore.inv_Pc']
        self['pore.inv_seq'] = sp.searchsorted(sp.unique(Pinv), Pinv)
        Tinv = self['throat.inv_Pc']
        self['throat.inv_seq'] = sp.searchsorted(sp.unique(Tinv), Tinv)

    def _do_one_inner_iteration(self, inv_val):
        r"""
        Determine which throats are invaded at a given applied capillary
        pressure.
        """
        # Generate a tlist containing boolean values for throat state
        Tinvaded = self['throat.entry_pressure'] <= inv_val
        # Add residual throats, if any, to list of invaded throats
        if 'throat.residual' in self.keys():
            Tinvaded = Tinvaded + self['throat.residual']
        # Find all pores that can be invaded at specified pressure
        [pclusters, tclusters] = self._net.find_clusters2(mask=Tinvaded,
                                                          t_labels=True)
        # Identify clusters connected to inlet sites
        inv_clusters = sp.unique(pclusters[self['pore.inlets']])
        inv_clusters = inv_clusters[inv_clusters >= 0]
        # Find pores on the invading clusters
        pmask = np.in1d(pclusters, inv_clusters)
        # Store current applied pressure in newly invaded pores
        pinds = (self['pore.inv_Pc'] == sp.inf) * (pmask)
        self['pore.inv_Pc'][pinds] = inv_val
        # Find throats on the invading clusters
        tmask = np.in1d(tclusters, inv_clusters)
        # Store current applied pressure in newly invaded throats
        tinds = (self['throat.inv_Pc'] == sp.inf) * (tmask)
        self['throat.inv_Pc'][tinds] = inv_val
        # Set residual pores and throats, if any, to invaded
        if 'pore.residual' in self.keys():
            self['pore.inv_Pc'][self['pore.residual']] = 0
        if 'throat.residual' in self.keys():
            self['throat.inv_Pc'][self['throat.residual']] = 0

    def plot_primary_drainage_curve(self, pore_volume='volume',
                                    throat_volume='volume'):
        r"""
        Plot the drainage curve as the  non-wetting phase saturation vs the
        applied capillary pressure.

        Parameters
        ----------
        pore_volume and throat_volume : string
            The dictionary key  where the pore and throat volume arrays
            are stored.  The defaults are 'pore.volume' and 'throat.volume'.
        """
        # Infer list of applied capillary pressures
        PcPoints = sp.unique(self['pore.inv_Pc'])
        # Get pore and throat volumes
        Pvol = self._net['pore.' + pore_volume]
        Tvol = self._net['throat.' + throat_volume]
        Total_vol = sp.sum(Pvol) + sp.sum(Tvol)
        # Find cumulative filled volume at each applied capillary pressure
        Vnwp_t = []
        Vnwp_p = []
        for Pc in PcPoints:
            Vnwp_p.append(sp.sum(Pvol[self['pore.inv_Pc'] <= Pc]))
            Vnwp_t.append(sp.sum(Tvol[self['throat.inv_Pc'] <= Pc]))
        # Combine throat and pore volumes into single total phase volume
        Vnwp_all = [Vnwp_p[i] + Vnwp_t[i] for i in range(0, len(PcPoints))]
        # Convert volumes to saturations by normalizing with total pore volume
        Snwp_all = [V/Total_vol for V in Vnwp_all]
        # Begin creating nicely formatted plot
        fig = plt.figure()
        plt.plot(PcPoints, Snwp_all, 'ko-')
        plt.ylim(ymin=0)
        plt.ylabel('Non-Wetting Phase Saturation')
        plt.xlabel('Capillary Pressure [Pa]')
        plt.title('Primary Drainage Curve')
        plt.grid(True)
        return fig
