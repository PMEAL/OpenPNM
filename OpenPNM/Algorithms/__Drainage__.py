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
    capillary pressures and invading all throats that are accessible and
    invadable at the given pressure.

    Parameters
    ----------
    network : OpenPNM Network Object
        The network upon which the simulation will be run

    name : string (optional)
        The name to apply to the Algorithm for quick identification

    Notes
    -----
    This algorithm creates several arrays and adds them to its own dictionary.
    These include:

    **'pore(throat).inv_Pc'** : The applied capillary pressure at which each
    pore(throat) was invaded.

    **'pore(throat).inv_Seq'** : The sequence each pore(throat) was filled.

    Examples
    --------
    >>> import OpenPNM as op
    >>> pn = op.Network.Cubic(shape=[20, 20, 20], spacing=10)
    >>> pn.add_boundary_pores(pores=pn.pores('top'),
    ...                       offset=[0, 0, 10],
    ...                       apply_label='boundary_top')
    >>> geo = op.Geometry.Stick_and_Ball(network=pn, pores=pn.Ps, throats=pn.Ts)
    >>> water = op.Phases.Water(network=pn)
    >>> air = op.Phases.Air(network=pn)
    >>> phys = op.Physics.Standard(network=pn, phase=water, geometry=geo)

    Once the basic Core objects are setup, the Algorithm can be created and
    and run as follows:

    >>> alg = op.Algorithms.Drainage(network=pn)
    >>> alg.setup(invading_phase=water, defending_phase=air)
    >>> alg.set_inlets(pores=pn.pores('boundary_top'))
    >>> alg.run()
    >>> data = alg.get_drainage_data()

    The ``data`` variable is a dictionary containing the numerical values of
    the resultant capillary pressure curve.  The available values can be
    inspected by typing ``data.keys()`` at the command line. These values can
    of course be plotted with matplotlib or exported to a graphing program of
    your choice.

    The final step is to utilize these results elsewhere in your OpenPNM
    simulation.  All Algorithms possess a method called ``return_results``
    which as the name suggests send the results to the correct locations.  For
    the 'Drainage Algorithm' this works as follows:

    >>> alg.return_results(Pc=5000)

    This command determines which pores and throats were filled at the applied
    capillary pressure of 5000, and creates 'pore.occupancy' and
    'throat.occupancy' arrays on the Phase object that was specfied as the
    'invading_phase' in the ``setup_method``.

    """

    def __init__(self, network, name=None):
        super().__init__(network=network, name=name)

    def setup(self,
              invading_phase,
              defending_phase,
              entry_pressure='throat.capillary_pressure',
              trapping=False,
              pore_filling=None,
              throat_filling=None,
              pore_volume='pore.volume',
              throat_volume='throat.volume'):
        r"""
        Used to specify necessary arguments to the simulation.  This method is
        useful for resetting the Algorithm or applying more explicit control.

        Parameters
        ----------
        invading_phase : OpenPNM Phase object
            The Phase object containing the physical properties of the invading
            fluid.

        defending_phase : OpenPNM Phase object
            The Phase object containing the physical properties of the defending
            fluid.

        entry_pressure : string (optional)
            The dictionary key on the Phase object where the throat entry
            pressure values can be found.  The default is
            'throat.capillary_pressure'.

        trapping : boolean (optional)
            Specifies whether defending phase trapping should be included or
            not. The default is False.  Note that defending phase outlets can
            be specified using the ``set_outlets`` method.  Otherwise it is
            assumed the defending phase has no outlets.

        pore_filling and throat_filling: string (optional)
            The dictionary key on the Physics object where the late pore or
            throat filling model is located. The default is None, meaning that
            a pore or throat is completely filled upon penetration.

        pore_volume and throat_volume : string (optional)
            The dictionary key on the Geometry object where the pore or throat
            volume data is located.  The defaults is 'pore.volume' and
            'throat.volume'.

        """
        self['throat.entry_pressure'] = invading_phase[entry_pressure]
        self['pore.inv_Pc'] = sp.inf
        self['throat.inv_Pc'] = sp.inf
        self['pore.trapped'] = sp.inf
        self['throat.trapped'] = sp.inf
        self['pore.inlets'] = False
        self['pore.outlets'] = False
        self['pore.residual'] = False
        self['throat.residual'] = False
        self._inv_phase = invading_phase
        self._def_phase = defending_phase
        self._trapping = trapping
        self._pore_filling = pore_filling
        self._throat_filling = throat_filling
        self._throat_volume = 'throat.volume'
        self._pore_volume = 'pore.volume'

    def set_inlets(self, pores=None, mode='add'):
        r"""
        Set the locations from which the invading phase enters the network.

        Parameters
        ----------
        pores : array_like
            An array of pore numbers that are initially filled with invading
            phase, and from which clusters of invading phase grow and invade
            into the network.

        mode : string
            Controls how the new values are handled.  Options are:

            **'add'** : (default) Adds the newly recieved locations to any
            existing locations.  This is useful for incrementally adding
            inlets.

            **'overwrite'** : Deletes any present locations and adds new ones.
            This is useful for fixing mistaken inlets, or rerunning the
            algorithm with different inlet locations.

            **'remove'** : Removes the received locations from the list of
            inlet pores.

            **'clear'** : Removes the existing inlets and ignores any specified
            locations.  This is equivalent to calling the method with a mode
            of 'overwrite' and pores = [] or None.

        Notes
        -----
        The 'inlet' pores are initially filled with invading fluid to start the
        simulation.  To avoid the capillary pressure curve showing a non-zero
        starting saturation at low pressures, it is necessary to create
        boundary pores that have zero-volume, and set these as the inlets.
        """
        Ps = self._parse_locations(pores)
        if mode in ['clear', 'overwrite']:
            self['pore.inlets'] = False
        if sum(self['pore.outlets'][Ps]) > 0:
            raise Exception('Some inlets are already defined as outlets')
        bool_val = True
        if mode is 'remove':
            bool_val = False
        self['pore.inlets'][Ps] = bool_val

    def set_outlets(self, pores=None, mode='add'):
        r"""
        Set the locations through which defending phase exits the network.
        This is only necessary if 'trapping'was set to True when setup was
        called.

        Parameters
        ----------
        pores : array_like
            An array of pore numbers where defending phase can exit.  Any
            defending phase that does not have access to these pores will be
            trapped.

        mode : string
            Controls how the new values are handled.  Options are:

            **'add'** : (default) Adds the newly recieved locations to any
            existing locations.  This is useful for incrementally adding
            outlets.

            **'overwrite'** : Deletes any present locations and adds new ones.
            This is useful for fixing mistaken outlets, or rerunning the
            algorithm with different outlet locations.

            **'remove'** : Removes the received locations from the list of
            outlet pores.

            **'clear'** : Removes the existing outlets and ignores any
            specified locations. This is equivalent to calling the method with
            a mode of 'overwrite' and pores = [] or None.

        """
        if self._trapping is False:
            raise Exception('Setting outlets is meaningless unless trapping ' +
                            'was set to True during setup')
        Ps = self._parse_locations(pores)
        if mode in ['clear', 'overwrite']:
            self['pore.outlets'] = False
        if sum(self['pore.inlets'][Ps]) > 0:
            raise Exception('Some outlets are already defined as inlets')
        bool_val = True
        if mode is 'remove':
            bool_val = False
        self['pore.outlets'][Ps] = bool_val

    def set_residual(self, pores=None, throats=None, mode='add'):
        r"""
        Specify locations of the residual invading (nonwetting) phase

        Parameters
        ----------
        pores and throats : array_like
            The pore and throat locations that are to be filled with invading
            phase at the beginning of the simulation.

        mode : string
            Controls how the new values are handled.  Options are:

            **'add'** : (default) Adds the newly recieved locations to any
            existing locations.  This is useful for incrementally adding
            outlets.

            **'overwrite'** : Deletes any present locations and adds new ones.
            This is useful for fixing mistaken outlets, or rerunning the
            algorithm with different outlet locations.

            **'remove'** : Removes the received locations from the list of
            residual pores and/or throats.  Both can be specified
            simultaneously.

            **'clear'** : Removes all existing residual locations.  This
            ignores any specified locations and is equivalent to calling the
            method with a mode of 'overwrite' and pores = [] or None.

        Notes
        -----
        Setting pores as initially filled only affects the saturation at the
        start of the capillary pressure curve since pre-filled pores do not
        contribute to the invasion process.  Setting throats as filled,
        however, has a significant impact on the subsequent invasion since
        these throats act as bridges to the percolation process.  Of course,
        they also contribute to the starting saturation as well.

        """
        bool_val = True
        if mode is 'clear':
            self['pore.residual'] = False
            self['throat.residual'] = False
            return
        if pores is not None:
            Ps = self._parse_locations(pores)
            if mode in ['clear', 'overwrite']:
                self['pore.residual'] = False
            if mode is 'remove':
                bool_val = False
            self['pore.residual'][Ps] = bool_val
        if throats is not None:
            Ts = self._parse_locations(throats)
            if mode in ['clear', 'overwrite']:
                self['throat.residual'] = False
            if mode is 'remove':
                bool_val = False
            self['throat.residual'][Ts] = bool_val

    def set_boundary_conditions(self, bc_type=None, pores=None, throats=None,
                                mode='add'):
        r"""
        Parameters
        ----------
        bc_type : string
            The type of boundary condition to apply.  Options are:

            **'inlets'** : For specifying where invading phase enters the
            Network

            **'outlets'** : For specifying where defending phase exits the
            Network if trapping is to be considered.

            **'residual'** : For specifying the pore and throat locations of
            existing residual invading phase in the Network at the start of the
            drainage.

        pores and thorats: array_like
            The pore (and throat) locations where the specified boundary
            condition is to be applied.  Note that the 'throats' argument
            is only valid for setting 'residual' locations, since 'inlets' and
            'outlets' locations can only be pores.

        mode : string
            Controls how the new values are handled.  Options are:

            **'add'** : (default) Adds the newly recieved locations to any
            existing locations.  This is useful for incrementally adding
            outlets.

            **'overwrite'** : Deletes any present locations and adds new ones.
            This is useful for fixing mistaken outlets, or rerunning the
            algorithm with different outlet locations.

            **'remove'** : Removes boundary conditions of the specified type
            from the given locations.

            **'clear'** : Removes existing conditions of the specified type
            and ignores any specified locations. If 'bc_type' is not specified
            then all conditions are removed.

        """
        if bc_type is None:
            if mode == 'clear':
                self['pore.inlets'] = False
                self['pore.outlets'] = False
                self['pore.residual'] = False
                self['throat.residual'] = False
                return
            else:
                raise Exception('\'bc_type\' cannot be None unless mode \
                                is \'clear\'')
        if bc_type == 'residual':
            self.set_residual(pores=pores, throats=throats, mode=mode)
        elif bc_type == 'inlets':
            self.set_inlets(pores=pores, mode=mode)
        elif bc_type == 'outlets':
            self.set_outlets(pores=pores, mode=mode)
        else:
            raise Exception('Unrecognized \'bc_type\' specified')

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
        self._inv_points = inv_points

        # Ensure inlets are set
        if sp.sum(self['pore.inlets']) == 0:
            raise Exception('Inlet pores have not been specified')

        # Ensure outlet pores are set if trapping is enabled
        if self._trapping:
            if sp.sum(self['pore.outlets']) == 0:
                raise Exception('Outlet pores have not been specified')

        # Generate curve from points
        for inv_val in self._inv_points:
            # Apply one applied pressure and determine invaded pores
            logger.info('Applying capillary pressure: ' + str(inv_val))
            self._apply_percolation(inv_val)
            if self._trapping:
                logger.info('Checking for trapping')
                self._check_trapping(inv_val)

        # Find invasion sequence values (to correspond with IP algorithm)
        Pinv = self['pore.inv_Pc']
        self['pore.inv_seq'] = sp.searchsorted(sp.unique(Pinv), Pinv)
        Tinv = self['throat.inv_Pc']
        self['throat.inv_seq'] = sp.searchsorted(sp.unique(Tinv), Tinv)

    def _check_trapping(self, inv_val):
        r"""
        Determine which pores and throats are trapped by invading phase.  This
        method is called by ``run`` if 'trapping' is set to True.
        """
        # Generate a list containing boolean values for throat state
        Tinvaded = self['throat.inv_Pc'] < sp.inf
        # Add residual throats, if any, to list of invaded throats
        Tinvaded = Tinvaded + self['throat.residual']
        # Invert logic to find defending throats
        Tdefended = ~Tinvaded
        [pclusters, tclusters] = self._net.find_clusters2(mask=Tdefended,
                                                          t_labels=True)
        # See which outlet pores remain uninvaded
        outlets = self['pore.outlets']*(self['pore.inv_Pc'] == sp.inf)
        # Identify clusters connected to remaining outlet sites
        def_clusters = sp.unique(pclusters[outlets])
        temp = sp.in1d(sp.unique(pclusters), def_clusters, invert=True)
        trapped_clusters = sp.unique(pclusters)[temp]
        trapped_clusters = trapped_clusters[trapped_clusters >= 0]

        # Find defending clusters NOT connected to the outlet pores
        pmask = np.in1d(pclusters, trapped_clusters)
        # Store current applied pressure in newly trapped pores
        pinds = (self['pore.trapped'] == sp.inf) * (pmask)
        self['pore.trapped'][pinds] = inv_val

        # Find throats on the trapped defending clusters
        tinds = self._net.find_neighbor_throats(pores=pinds,
                                                mode='intersection')
        self['throat.trapped'][tinds] = inv_val
        self['throat.entry_pressure'][tinds] = 1000000

    def _apply_percolation(self, inv_val):
        r"""
        Determine which pores and throats are invaded at a given applied
        capillary pressure.  This method is called by ``run``.
        """
        # Generate a list containing boolean values for throat state
        Tinvaded = self['throat.entry_pressure'] <= inv_val
        # Add residual throats, if any, to list of invaded throats
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
        if sp.any(self['pore.residual']):
            self['pore.inv_Pc'][self['pore.residual']] = 0
        if sp.any(self['throat.residual']):
            self['throat.inv_Pc'][self['throat.residual']] = 0

    def get_drainage_data(self):
        r"""
        Obtain the numerical values of the resultant capillary pressure curve.

        Parameters
        ----------
        None

        Returns
        -------
        A dictionary containing arrays of applied capillary pressures and
        various phase saturations.  The dictionary keys explain the content of
        each array.

        Notes
        -----
        This method assumes that the pore and throat volumes are stored under
        the keys 'pore.volume' and 'throat.volume'.  This cannot be customized
        at this time.
        """
        # Infer list of applied capillary pressures
        PcPoints = self._inv_points
        if PcPoints[-1] == sp.inf:  # Remove infinity from PcPoints if present
            PcPoints = PcPoints[:-1]
        # Get pore and throat volumes
        Pvol = self._net[self._pore_volume]
        Tvol = self._net[self._throat_volume]
        Total_vol = sp.sum(Pvol) + sp.sum(Tvol)
        # Find cumulative filled volume at each applied capillary pressure
        Vnwp_t = []
        Vnwp_p = []
        Vnwp_all = []
        for Pc in PcPoints:
            # Calculate filled pore volumes
            p_inv = self['pore.inv_Pc'] <= Pc
            if self._pore_filling is None:
                Vp = sp.sum(Pvol[p_inv])
            else:
                Vp = self._calc_fractional_filling(pressure=Pc,
                                                   element='pore')
                Vp = sp.sum(Vp[p_inv])
            # Calculate filled throat volumes
            t_inv = self['throat.inv_Pc'] <= Pc
            if self._throat_filling is None:
                Vt = sp.sum(Tvol[t_inv])
            else:
                Vt = self._calc_fractional_filling(pressure=Pc,
                                                   element='throat')
                Vt = sp.sum(Vt[t_inv])
            Vnwp_p.append(Vp)
            Vnwp_t.append(Vt)
            Vnwp_all.append(Vp + Vt)
        # Convert volumes to saturations by normalizing with total pore volume
        Snwp_all = [V/Total_vol for V in Vnwp_all]
        data = {}
        data['capillary_pressure'] = PcPoints
        data['invading_phase_saturation'] = Snwp_all
        data['defending_phase_saturation'] = [1-s for s in Snwp_all]
        return data

    def _calc_fractional_filling(self, element, pressure):
        r"""
        Calculates the fractional filling of each pore or throat as the
        supplied capillary pressure

        Parameters
        ----------
        element : string
            Can either be 'pore' or 'throat' indicating which type of element
            to be calculated

        pressure : float
            The value of the capillary pressure to apply.  This value is sent
            the to Physics model stored in the ``_pore(throat)_filling``
            attribute of the object.  This attribute is set during the call
            to ``get_drainage_data``.

        Notes
        -----
        The 'pore(throat)_filling' model must accept the applied capillary
        pressure as 'Pc'.  This is not customizable at the moment.
        """
        if element == 'pore':
            key = self._pore_filling
            vol = self._pore_volume
        elif element == 'throat':
            key = self._throat_filling
            vol = self._throat_volume
        else:
            raise Exception('element must be either \'pore\' or \'throat\'')
        Snwp = sp.zeros((self._count(element),))
        for phys in self._inv_phase._physics:
            # Update Physics model with current Pc
            temp_Pc = phys.models[key]['Pc']  # Store old Pc
            phys.models[key]['Pc'] = pressure
            # Regenerate Physics model and capture output locally
            Snwp[phys.Pnet] = phys.models[key].run()
            # Re-populate the residual element with the non-wetting phase
            if sp.any(self[element+'.residual']):
                Snwp[self[element+'.residual']] = 1.0
            phys.models[key]['Pc'] = temp_Pc  # Return old Pc
        V = self._net[vol]*Snwp
        return V

    def plot_drainage_curve(self,
                            data=None,
                            x_values='capillary_pressure',
                            y_values='invading_phase_saturation'):
        r"""
        Plot the drainage curve as the non-wetting phase saturation vs the
        applied capillary pressure.

        Parameters
        ----------
        data : dictionary of arrays
            This dictionary should be obtained from the ``get_drainage_data``
            method.

        x_values and y_values : string
            The dictionary keys of the arrays containing the x-values and
            y-values

        """
        # Begin creating nicely formatted plot
        if data is None:
            data = self.get_drainage_data()
        xdata = data[x_values]
        ydata = data[y_values]
        fig = plt.figure()
        plt.plot(xdata, ydata, 'ko-')
        plt.ylabel(y_values)
        plt.xlabel(x_values)
        plt.grid(True)
        if sp.amax(xdata) <= 1:
            plt.xlim(xmin=0, xmax=1)
        if sp.amax(ydata) <= 1:
            plt.ylim(ymin=0, ymax=1)
        return fig

    def return_results(self, Pc):
        r"""
        This method determines which pores and throats are filled with non-
        wetting phase at the specified capillary pressure, and creates or
        updates 'pore.occupancy' and 'throat.occupancy' arrays on the
        associated Phase objects. Invasion pressure and sequence are also sent
        to the invading phase.

        Parameters
        ----------
        Pc : scalar
            The capillary pressure for which an invading phase configuration is
            required.

        Returns
        -------

        Notes
        -----
        The Phase object that receives the updated occupancy arrays is the
        one that was specified as the 'invading_phase' when ``setup`` was
        called. The defending phase has the opposite occupancy values and
        partial occupancy so that summing occupancy for both phases equals
        1.0 for every pore.
        """
        Psatn = self['pore.inv_Pc'] <= Pc
        Tsatn = self['throat.inv_Pc'] <= Pc
        self._inv_phase['pore.occupancy'] = sp.array(Psatn, dtype=float)
        self._inv_phase['throat.occupancy'] = sp.array(Tsatn, dtype=float)
        self._def_phase['pore.occupancy'] = sp.array(~Psatn, dtype=float)
        self._def_phase['throat.occupancy'] = sp.array(~Tsatn, dtype=float)
        self._inv_phase['pore.inv_Pc'] = self['pore.inv_Pc']
        self._inv_phase['pore.inv_seq'] = self['pore.inv_seq']
        self._inv_phase['throat.inv_Pc'] = self['throat.inv_Pc']
        self._inv_phase['throat.inv_seq'] = self['throat.inv_seq']
        if self._pore_filling:
            Vp = self._calc_fractional_filling(element='pore', pressure=Pc)
            Sp = Vp/self._net[self._pore_volume]
            # For pores that have zero volume (i.e. boundary pores in some cases)
            # fractional filling is meaningless so use the standard occupancy
            zero_ps = self._net[self._pore_volume] == 0.0
            Sp[zero_ps] = self._inv_phase['pore.occupancy'][zero_ps]
            self._inv_phase['pore.partial_occupancy'] = Sp
            self._def_phase['pore.partial_occupancy'] = 1 - Sp
        if self._throat_filling:
            Vt = self._calc_fractional_filling(element='throat', pressure=Pc)
            St = Vt/self._net[self._throat_volume]
            # For throats that have zero volume (i.e. boundary throats in some cases)
            # fractional filling is meaningless so use the standard occupancy
            zero_ts = self._net[self._throat_volume] == 0.0
            St[zero_ts] = self._inv_phase['throat.occupancy'][zero_ts]
            self._inv_phase['throat.partial_occupancy'] = St
            self._def_phase['throat.partial_occupancy'] = 1 - St
