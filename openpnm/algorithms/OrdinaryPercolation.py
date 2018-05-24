import scipy as sp
import numpy as np
from openpnm.algorithms import GenericAlgorithm
from openpnm.topotools import site_percolation, bond_percolation
from openpnm.topotools import remove_isolated_clusters, ispercolating
from openpnm.core import logging
logger = logging.getLogger(__name__)

default_settings = {'access_limited': True,
                    'mode': 'bond',
                    'pore_entry_threshold': 'pore.entry_pressure',
                    'throat_entry_threshold': 'throat.entry_pressure',
                    'pore_volume': '',
                    'throat_volume': ''}


class OrdinaryPercolation(GenericAlgorithm):
    r"""

    """

    def __init__(self, settings={}, **kwargs):
        super().__init__(**kwargs)
        self.settings.update(default_settings)
        # Apply user settings, if any
        self.settings.update(settings)
        # Use the reset method to initialize all arrays
        self.reset()

    def setup(self,
              phase=None,
              access_limited=None,
              mode='',
              throat_entry_pressure='',
              pore_entry_pressure='',
              pore_volume='',
              throat_volume=''):
        r"""
        Used to specify necessary arguments to the simulation.  This method is
        useful for resetting the algorithm or applying more explicit control.

        Parameters
        ----------
        phase : OpenPNM Phase object
            The Phase object containing the physical properties of the invading
            fluid.

        access_limited : boolean
            If ``True`` the invading phase can only enter the network from the
            invasion sites specified with ``set_inlets``.  Otherwise, invading
            clusters can appear anywhere in the network.  This second case is
            the normal *ordinary percolation* in the traditional sense.

        mode : string
            Specifies the type of percolation process to simulate.  Options
            are:

            **'bond'** - The percolation process is controlled by bond entry
            thresholds.

            **'site'** - The percolation process is controlled by site entry
            thresholds.

        pore_entry_pressure : string
            The dictionary key on the Phase object where the pore entry
            pressure values are stored.  The default is
            'pore.capillary_pressure'.  This is only accessed if the ``mode``
            is set to site percolation.

        throat_entry_pressure : string
            The dictionary key on the Phase object where the throat entry
            pressure values are stored.  The default is
            'throat.capillary_pressure'.  This is only accessed if the ``mode``
            is set to bond percolation.

        'pore_volume' : string
            The dictionary key containing the pore volume information.

        'throat_volume' : string
            The dictionary key containing the pore volume information.

        """
        if phase:
            self.settings['phase'] = phase.name
        if throat_entry_pressure:
            self.settings['throat_entry_pressure'] = throat_entry_pressure
            phase = self.project.find_phase(self)
            self['throat.entry_pressure'] = phase[throat_entry_pressure]
        if pore_entry_pressure:
            self.settings['pore_entry_pressure'] = pore_entry_pressure
            phase = self.project.find_phase(self)
            self['pore.entry_pressure'] = phase[pore_entry_pressure]
        if mode:
            self.settings['mode'] = mode
        if access_limited is not None:
            self.settings['access_limited'] = access_limited
        if pore_volume:
            self.settings['pore_volume'] = pore_volume
        if throat_volume:
            self.settings['throat_volume'] = throat_volume

    def reset(self):
        r"""
        Resets the various data arrays on the object back to their original
        state. This is useful for repeating a simulation at different inlet
        conditions, or invasion points for instance.
        """
        self['pore.invasion_pressure'] = np.inf
        self['throat.invasion_pressure'] = np.inf
        self['pore.invasion_sequence'] = -1
        self['throat.invasion_sequence'] = -1
        self['pore.inlets'] = False
        self['pore.outlets'] = False
        self['pore.residual'] = False
        self['throat.residual'] = False

    def set_inlets(self, pores=[], overwrite=False):
        r"""
        Set the locations from which the invader enters the network

        Parameters
        ----------
        pores : array_like
            Locations that are initially filled with invader, from which
            clusters grow and invade into the network

        overwrite : boolean
            If ``True`` then all existing inlet locations will be removed and
            then the supplied locations will be added.  If ``False`` (default),
            then supplied locations are added to any already existing inlet
            locations.


        """
        Ps = self._parse_indices(pores)
        if np.sum(self['pore.outlets'][Ps]) > 0:
            raise Exception('Some inlets are already defined as outlets')
        if overwrite:
            self['pore.inlets'] = False
        self['pore.inlets'][Ps] = True
        self['pore.invasion_pressure'][Ps] = 0
        self['pore.invasion_sequence'][Ps] = 0

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
        if self.settings['trapping'] is False:
            logger.warning('Setting outlets is meaningless unless trapping ' +
                           'was set to True during setup')
        Ps = self._parse_indices(pores)
        if np.sum(self['pore.inlets'][Ps]) > 0:
            raise Exception('Some outlets are already defined as inlets')
        if overwrite:
            self['pore.outlets'] = False
        self['pore.outlets'][Ps] = True

    def set_residual(self, pores=[], throats=[], overwrite=False):
        r"""
        Specify locations of any residual invader.  These locations are set
        to invaded at the start of the simulation.

        Parameters
        ----------
        pores : array_like
            The pores locations that are to be filled with invader at the
            beginning of the simulation.

        throats : array_like
            The throat locations that are to be filled with invader at the
            beginning of the simulation.

        overwrite : boolean
            If ``True`` then all existing inlet locations will be removed and
            then the supplied locations will be added.  If ``False``, then
            supplied locations are added to any already existing locations.

        """
        Ps = self._parse_indices(pores)
        if overwrite:
            self['pore.residual'] = False
        self['pore.residual'][Ps] = True
        Ts = self._parse_indices(throats)
        if overwrite:
            self['throat.residual'] = False
        self['throat.residual'][Ts] = True

    def get_percolation_threshold(self):
        r"""
        """
        if np.sum(self['pore.inlets']) == 0:
            raise Exception('Inlet pores must be specified first')
        if np.sum(self['pore.outlets']) == 0:
            raise Exception('Outlet pores must be specified first')
        else:
            Pout = self['pore.outlets']
        # Do a simple check of pressures on the outlet pores first...
        if self.settings['access_limited']:
            thresh = np.amin(self['pore.invasion_pressure'][Pout])
        else:
            raise Exception('This is currently only implemented for access ' +
                            'limited simulations')
        return thresh

    def is_percolating(self, applied_pressure):
        r"""
        Returns a True or False value to indicate if a percolating cluster
        spans between the inlet and outlet pores that were specified at the
        given applied pressure.

        Parameters
        ----------
        applied_pressure : scalar, float
            The pressure at which percolation should be checked

        Returns
        -------
        A simple boolean True or False if percolation has occured or not.

        """
        if np.sum(self['pore.inlets']) == 0:
            raise Exception('Inlet pores must be specified first')
        else:
            Pin = self['pore.inlets']
        if np.sum(self['pore.outlets']) == 0:
            raise Exception('Outlet pores must be specified first')
        else:
            Pout = self['pore.outlets']
        # Do a simple check of pressures on the outlet pores first...
        if np.amin(self['pore.invasion_pressure'][Pout]) > applied_pressure:
            val = False
        else:  # ... and do a rigorous check only if necessary
            mask = self['throat.invasion_pressure'] < applied_pressure
            am = self.project.network.create_adjacency_matrix(weights=mask,
                                                              fmt='coo')
            val = ispercolating(am=am, mode=self.settings['mode'],
                                inlets=Pin, outlets=Pout)
        return val

    def run(self, points=25, start=None, stop=None):
        r"""
        Runs the percolation algorithm to determine which pores and throats
        will be invaded at each given pressure point.

        Parameters
        ----------
        points: int or array_like
            An array containing the pressure points to apply.  If a scalar is
            given then an array will be generated with the given number of
            points spaced between the lowest and highest values of throat
            entry pressures using logarithmic spacing.  To specify low and
            high pressure points use the ``start`` and ``stop`` arguments.

        start : int
            The optional starting point to use when generating pressure points.

        stop : int
            The optional stopping point to use when generating pressure points.

        """
        phase = self.project.find_phase(self)
        # Parse inputs and generate list of invasion points if necessary
        if self.settings['mode'] == 'bond':
            self['throat.entry_pressure'] = \
                phase[self.settings['throat_entry_threshold']]
            if start is None:
                start = sp.amin(self['throat.entry_pressure'])*0.95
            if stop is None:
                stop = sp.amax(self['throat.entry_pressure'])*1.05

        elif self.settings['mode'] == 'site':
            self['pore.entry_pressure'] = \
                phase[self.settings['pore_entry_threshold']]
            if start is None:
                start = sp.amin(self['pore.entry_pressure'])*0.95
            if stop is None:
                stop = sp.amax(self['pore.entry_pressure'])*1.05
        else:
            raise Exception('Percolation type has not been set')
        if type(points) is int:
            points = sp.logspace(start=sp.log10(max(1, start)),
                                 stop=sp.log10(stop), num=points)

        # Ensure pore inlets have been set IF access limitations is True
        if self.settings['access_limited']:
            if sp.sum(self['pore.inlets']) == 0:
                raise Exception('Inlet pores must be specified first')
            else:
                Pin = self['pore.inlets']

        # Generate curve from points
        conns = self.project.network['throat.conns']
        for inv_val in points:
            if self.settings['mode'] == 'bond':
                t_invaded = self['throat.entry_pressure'] <= inv_val
                labels = bond_percolation(conns, t_invaded)
            elif self.settings['mode'] == 'site':
                p_invaded = self['pore.entry_pressure'] <= inv_val
                labels = site_percolation(conns, p_invaded)

            # Optionally remove clusters not connected to the inlets
            if self.settings['access_limited']:
                labels = remove_isolated_clusters(labels=labels,
                                                  inlets=Pin)

            # Store current applied pressure in newly invaded pores
            pinds = (self['pore.invasion_pressure'] == sp.inf) * \
                    (labels.sites >= 0)
            self['pore.invasion_pressure'][pinds] = inv_val
            # Store current applied pressure in newly invaded throats
            tinds = (self['throat.invasion_pressure'] == sp.inf) * \
                    (labels.bonds >= 0)
            self['throat.invasion_pressure'][tinds] = inv_val

        # Convert invasion pressures in sequence values
        Pinv = self['pore.invasion_pressure']
        Tinv = self['throat.invasion_pressure']
        Pseq = sp.searchsorted(sp.unique(Pinv), Pinv)
        Tseq = sp.searchsorted(sp.unique(Tinv), Tinv)
        self['pore.invasion_sequence'] = Pseq
        self['throat.invasion_sequence'] = Tseq

    def results(self, Pc):
        r"""
        This method determines which pores and throats are filled with invading
        phase at the specified capillary pressure, and creates several arrays
        indicating the occupancy status of each pore and throat for the given
        pressure.

        Parameters
        ----------
        Pc : scalar
            The capillary pressure for which an invading phase configuration
            is desired.

        Returns
        -------
        A dictionary containing an assortment of data about distribution
        of the invading phase at the specified capillary pressure.  The data
        include:

        **'pore.occupancy'** : A value between 0 and 1 indicating the
        fractional volume of each pore that is invaded.  If no late pore
        filling model was applied, then this will only be integer values
        (either filled or not).

        **'throat.occupancy'** : The same as 'pore.occupancy' but for throats.

        This dictionary can be passed directly to the ``update`` method of
        the *Phase* object. These values can then be accessed by models
        or algorithms.

        """
        Psatn = self['pore.invasion_pressure'] <= Pc
        Tsatn = self['throat.invasion_pressure'] <= Pc
        inv_phase = {}
        inv_phase['pore.occupancy'] = sp.array(Psatn, dtype=float)
        inv_phase['throat.occupancy'] = sp.array(Tsatn, dtype=float)
        return inv_phase
