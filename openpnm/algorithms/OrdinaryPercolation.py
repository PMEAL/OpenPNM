import scipy as sp
import porespy as ps
import numpy as np
import scipy.sparse as sprs
import scipy.sparse.csgraph as csg
from collections import namedtuple
from openpnm.algorithms import GenericAlgorithm
from openpnm.core import logging
logger = logging.getLogger(__name__)


class OrdinaryPercolation(GenericAlgorithm):

    def __init__(self, settings={}, **kwargs):
        r"""
        """
        super().__init__(**kwargs)
        defaults = {'access_limited': True,
                    'mode': 'bond',
                    'pore_entry_pressure': 'pore.capillary_pressure',
                    'throat_entry_pressure': 'throat.capillary_pressure',
                    'pore_volume': '',
                    'throat_volume': ''}
        self.settings.update(defaults)
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
            the supplied locations are added to any already existing locations.

        """
        Ps = self._parse_indices(pores)
        if sum(self['pore.outlets'][Ps]):
            raise Exception('Some given indices are already set as outlets')
        if overwrite:
            self['pore.inlets'] = False
        self['pore.inlets'][Ps] = True

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
        Ps = self._parse_indices(pores)
        if sum(self['pore.inlets'][Ps]):
            raise Exception('Some given indices are already set as inlets')
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
            If ``True`` then all existing residual locations will be removed
            and then the supplied locations will be added.  If ``False``, the
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

    def run(self, points=25, start=None, stop=None):
        r"""
        Runs the percolation algorithm to determine which pores and throats
        will be invaded at each given

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
                phase[self.settings['throat_entry_pressure']]
            if type(points) is int:
                if start is None:
                    start = sp.amin(self['throat.entry_pressure'])*0.95
                if stop is None:
                    stop = sp.amax(self['throat.entry_pressure'])*1.05
                points = sp.logspace(start=sp.log10(max(1, start)),
                                     stop=sp.log10(stop),
                                     num=points)
        elif self.settings['mode'] == 'site':
            self['pore.entry_pressure'] = \
                phase[self.settings['pore_entry_pressure']]
            if type(points) is int:
                if start is None:
                    start = sp.amin(self['pore.entry_pressure'])*0.95
                if stop is None:
                    stop = sp.amax(self['pore.entry_pressure'])*1.05
                points = sp.logspace(start=sp.log10(max(1, start)),
                                     stop=sp.log10(stop),
                                     num=points)
        else:
            raise Exception('Percolation type has not been set')
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
                labels = self.bond_percolation(conns, t_invaded)
            elif self.settings['mode'] == 'site':
                p_invaded = self['pore.entry_pressure'] <= inv_val
                labels = self.site_percolation(conns, p_invaded)

            # Optionally remove clusters not connected to the inlets
            if self.settings['access_limited']:
                labels = self.remove_isolated_clusters(labels=labels,
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

    def get_percolation_threshold(self):
        r"""
        """
        if sp.sum(self['pore.inlets']) == 0:
            raise Exception('Inlet pores must be specified first')
        else:
            Pin = self['pore.inlets']
        if sp.sum(self['pore.outlets']) == 0:
            raise Exception('Outlet pores must be specified first')
        else:
            Pout = self['pore.outlets']
        # Do a simple check of pressures on the outlet pores first...
        if self.settings['access_limited']:
            thresh = sp.amin(self['pore.invasion_pressure'][Pout])
        else:
            raise Exception('This is currently only implemented for access ' +
                            'limited simulations')
        return thresh

    def is_percolating(self, applied_pressure):
        r"""
        Returns a True of False value to indicate if a percolating cluster
        spans the between the inlet and outlet pores that were specified.

        Parameters
        ----------
        applied_pressure : scalar, float
            The pressure for which percolation should be checked.

        Returns
        -------
        A simple boolean True or False if percolation has occured or not.

        """
        if sp.sum(self['pore.inlets']) == 0:
            raise Exception('Inlet pores must be specified first')
        else:
            Pin = self['pore.inlets']
        if sp.sum(self['pore.outlets']) == 0:
            raise Exception('Outlet pores must be specified first')
        else:
            Pout = self['pore.outlets']
        # Do a simple check of pressures on the outlet pores first...
        if sp.amin(self['pore.invasion_pressure'][Pout]) > applied_pressure:
            val = False
        else:  # ... and do a rigorous check only if necessary
            mask = self['throat.invasion_pressure'] < applied_pressure
            am = self.project.network.create_adjacency_matrix(weights=mask,
                                                              fmt='coo')
            val = self._is_percolating(am=am, mode=self.settings['mode'],
                                       inlets=Pin, outlets=Pout)
        return val

    def _is_percolating(self, am, inlets, outlets, mode='site'):
        r"""
        Determines if a percolating clusters exists in the network spanning
        the given inlet and outlet sites

        Parameters
        ----------
        am : adjacency_matrix
            The adjacency matrix with the ``data`` attribute indicating
            if a bond is occupied or not

        inlets : array_like
            An array of indices indicating which sites are part of the inlets

        outlets : array_like
            An array of indices indicating which sites are part of the outlets

        mode : string
            Indicates which type of percolation to apply, either `'site'` or
            `'bond'`

        """
        if am.format is not 'coo':
            am = am.to_coo()
        ij = sp.vstack((am.col, am.row)).T
        if mode.startswith('site'):
            occupied_sites = sp.zeros(shape=am.shape[0], dtype=bool)
            occupied_sites[ij[am.data].flatten()] = True
            clusters = self.site_percolation(ij, occupied_sites)
        elif mode.startswith('bond'):
            occupied_bonds = am.data
            clusters = self.bond_percolation(ij, occupied_bonds)
        ins = sp.unique(clusters.sites[inlets])
        if ins[0] == -1:
            ins = ins[1:]
        outs = sp.unique(clusters.sites[outlets])
        if outs[0] == -1:
            outs = outs[1:]
        hits = sp.in1d(ins, outs)
        return sp.any(hits)

    def remove_isolated_clusters(self, labels, inlets):
        r"""
        Finds cluster labels not attached to the inlets, and sets them to
        unoccupied (-1)

        Parameters
        ----------
        labels : tuple of site and bond labels
            This information is provided by the ``site_percolation`` or
            ``bond_percolation`` functions

        inlets : array_like
            A list of which sites are inlets.  Can be a boolean mask or an
            array of indices.

        Returns
        -------
        A tuple containing a list of site and bond labels, with all clusters
        not connected to the inlet sites set to not occupied.

        """
        # Identify clusters of invasion sites
        inv_clusters = sp.unique(labels.sites[inlets])
        # Remove cluster numbers == -1, if any
        inv_clusters = inv_clusters[inv_clusters >= 0]
        # Find all pores in invading clusters
        p_invading = np.in1d(labels.sites, inv_clusters)
        labels.sites[~p_invading] = -1
        t_invading = np.in1d(labels.bonds, inv_clusters)
        labels.bonds[~t_invading] = -1
        return labels

    def site_percolation(self, ij, occupied_sites):
        r"""
        Calculates the site and bond occupancy status for a site percolation
        process given a list of occupied sites.

        Parameters
        ----------
        ij : array_like
            An N x 2 array of [site_A, site_B] connections.  If two connected
            sites are both occupied they are part of the same cluster, as it
            the bond connecting them.

        occupied_sites : boolean
            A list indicating whether sites are occupied or not

        Returns
        -------
        A tuple containing a list of site and bond labels, indicating which
        cluster each belongs to.  A value of -1 indicates unoccupied.

        Notes
        -----
        The ``connected_components`` function of scipy.csgraph will give ALL
        sites a cluster number whether they are occupied or not, so this
        function essentially adjusts the cluster numbers to represent a
        percolation process.

        """
        from collections import namedtuple
        Np = sp.size(occupied_sites)
        occupied_bonds = sp.all(occupied_sites[ij], axis=1)
        adj_mat = sprs.csr_matrix((occupied_bonds, (ij[:, 0], ij[:, 1])),
                                  shape=(Np, Np))
        adj_mat.eliminate_zeros()
        clusters = csg.connected_components(csgraph=adj_mat, directed=False)[1]
        clusters[~occupied_sites] = -1
        s_labels = ps.tools.make_contiguous(clusters + 1) - 1
        b_labels = sp.amin(s_labels[ij], axis=1)
        tup = namedtuple('cluster_labels', ('sites', 'bonds'))
        return tup(s_labels, b_labels)

    def bond_percolation(self, ij, occupied_bonds):
        r"""
        Calculates the site and bond occupancy status for a bond percolation
        process given a list of occupied bonds.

        Parameters
        ----------
        ij : array_like
            An N x 2 array of [site_A, site_B] connections.  A site is
            considered occupied if any of it's connecting bonds are occupied.

        occupied_bonds: boolean
            A list indicating whether a bond is occupied or not

        Returns
        -------
        A tuple contain a list of site and bond labels, indicating which
        cluster each belongs to.  A value of -1 indicates uninvaded.

        Notes
        -----
        The ``connected_components`` function of scipy.csgraph will give ALL
        sites a cluster number whether they are occupied or not, so this
        function essentially adjusts the cluster numbers to represent a
        percolation process.
        """
        from collections import namedtuple
        Np = sp.amax(ij) + 1
        adj_mat = sprs.csr_matrix((occupied_bonds, (ij[:, 0], ij[:, 1])),
                                  shape=(Np, Np))
        adj_mat.eliminate_zeros()
        clusters = csg.connected_components(csgraph=adj_mat, directed=False)[1]
        valid_clusters = sp.bincount(clusters) > 1
        mapping = -sp.ones(shape=(clusters.max()+1, ), dtype=int)
        mapping[valid_clusters] = sp.arange(0, valid_clusters.sum())
        s_labels = mapping[clusters]
        b_labels = sp.amin(s_labels[ij], axis=1)
        tup = namedtuple('cluster_labels', ('sites', 'bonds'))
        return tup(s_labels, b_labels)

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

    def get_percolation_data(self):
        r"""
        Obtain the numerical values of the calculated percolation curve

        Returns
        -------
        A named-tuple containing arrays of applied capillary pressures and
        invading phase saturation.

        """
        net = self.project.network
        # Infer list of applied capillary pressures
        points = np.unique(self['throat.invasion_pressure'])
        # Add a low pressure point to the list to improve graph
        points = np.concatenate(([0], points))
        if points[-1] == np.inf:  # Remove infinity from PcPoints if present
            points = points[:-1]
        # Get pore and throat volumes
        if self.settings['pore_volume']:
            Pvol = net[self.settings['pore_volume']]
        else:
            Pvol = sp.ones(shape=(self.Np, ), dtype=int)
        if self.settings['throat_volume']:
            Tvol = net[self.settings['throat_volume']]
        else:
            Tvol = sp.zeros(shape=(self.Nt, ), dtype=int)
        Total_vol = np.sum(Pvol) + np.sum(Tvol)
        # Find cumulative filled volume at each applied capillary pressure
        Vnwp_t = []
        Vnwp_p = []
        Vnwp_all = []
        for p in points:
            # Calculate filled pore volumes
            p_inv = self['pore.invasion_pressure'] <= p
            Vp = np.sum(Pvol[p_inv])
            # Calculate filled throat volumes
            t_inv = self['throat.invasion_pressure'] <= p
            Vt = np.sum(Tvol[t_inv])
            Vnwp_p.append(Vp)
            Vnwp_t.append(Vt)
            Vnwp_all.append(Vp + Vt)
        # Convert volumes to saturations by normalizing with total pore volume
        Snwp_all = [V/Total_vol for V in Vnwp_all]
        pc_curve = namedtuple('pc_curve', ('Pcap', 'Snwp'))
        data = pc_curve(points, Snwp_all)
        return data

    def plot_percolation_curve(self):
        r"""
        Plot the percolation curve as the invader volume or number fraction vs
        the applied capillary pressure.

        """
        # Begin creating nicely formatted plot
        data = self.get_percolation_data()
        xdata = data.Pcap
        ydata = data.Snwp
        fig = plt.figure()
        plt.semilogx(xdata, ydata, 'ko-')
        plt.ylabel('Invading Phase Saturation')
        plt.xlabel('Capillary Pressure')
        plt.grid(True)
        if np.amax(xdata) <= 1:
            plt.xlim(xmin=0, xmax=1)
        if np.amax(ydata) <= 1:
            plt.ylim(ymin=0, ymax=1)
        return fig