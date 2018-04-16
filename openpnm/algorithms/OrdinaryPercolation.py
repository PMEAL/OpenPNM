import scipy as sp
import numpy as np
import scipy.sparse as sprs
import scipy.sparse.csgraph as csg
from collections import namedtuple
from openpnm.algorithms import GenericPercolation
from openpnm.core import logging
logger = logging.getLogger(__name__)
tup = namedtuple('cluster_labels', ('sites', 'bonds'))


class OrdinaryPercolation(GenericPercolation):

    def __init__(self, settings={}, **kwargs):
        r"""
        """
        super().__init__(**kwargs)
        self.settings.update({'trapping': False,
                              'access_limited': True,
                              'mode': 'bond'})
        self.settings.update(settings)
        self.reset()

    def setup(self,
              phase=None,
              trapping=None,
              access_limited=None,
              mode='',
              entry_pressure=''):
        r"""
        Used to specify necessary arguments to the simulation.  This method is
        useful for resetting the algorithm or applying more explicit control.

        Parameters
        ----------
        phase : OpenPNM Phase object
            The Phase object containing the physical properties of the invading
            fluid.

        entry_pressure : string
            The dictionary key on the Phase object where the throat entry
            pressure values are stored.  The default is
            'throat.capillary_pressure'.

        trapping : boolean
            Specifies whether defending phase trapping should be included or
            not. The default is False.  Note that defending phase outlets can
            be specified using the ``set_outlets`` method.  Otherwise it is
            assumed the defending phase has no outlets.

        access_limited : boolean
            If ``True`` the invading phase can only enter the network from the
            invasion sites specified with ``set_inlets``.  Otherwise, invading
            clusters can appear anywhere in the network.  This second case is
            the normal *ordinary percolation*.

        mode : string
            Specifies they type of percolation process to simulate.  Options
            are:

            **'bond'** - The percolation process is controlled by bond entry
            thresholds.

            **'site'** - The percolation process is controlled by site entry
            thresholds.

        """
        if phase:
            self.settings['phase'] = phase.name
        if entry_pressure:
            self.settings['entry_pressure'] = entry_pressure
            phase = self.project.find_phase(self)
            self['throat.entry_pressure'] = phase[entry_pressure]
        if trapping is not None:
            self.settings['trapping'] = trapping
        if access_limited is not None:
            self.settings['access_limited'] = access_limited

    def run(self, points=25):
        r"""
        Runs the percolation algorithm to determine which pores and throats
        will be invaded at each given

        Parameters
        ----------
        points: int or array_like
            An array containing the pressure points to apply.  If a scalar is
            given then an array will be generated with the given number of
            points spaced between the lowest and highest values of throat
            entry pressures.

        """
        phase = self.project.find_phase(self)
        self['throat.entry_pressure'] = phase['throat.capillary_pressure']
        # Parse inputs and generate list of invasion points if necessary
        if type(points) is int:
            logger.info('Generating list of invasion pressures')
            # Nudge values up and down
            if self.settings['mode'] == 'bond':
                min_p = sp.amin(self['throat.entry_pressure'])*0.95
                max_p = sp.amax(self['throat.entry_pressure'])*1.05
            elif self.settings['mode'] == 'site':
                min_p = sp.amin(self['pore.entry_pressure'])*0.95
                max_p = sp.amax(self['pore.entry_pressure'])*1.05
            else:
                raise Exception('Percolation type has not been set')
            points = sp.logspace(sp.log10(min_p), sp.log10(max_p), points)

        # Generate curve from points
        for inv_val in points:

            # Apply one applied pressure and determine invaded pores
            conns = self.project.network['throat.conns']
            labels = self._apply_percolation(inv_val)

            # Optionally remove clusters not connected to the inlets
            if self.settings['access_limited']:
                labels = self.remove_isolated_clusters

            if self.settings['trapping']:
                t_defended = labels.bonds == -1
                temp = self.bond_percolation(conns, t_defended)

            # Store current applied pressure in newly invaded pores
            pinds = (self['pore.invasion_pressure'] == sp.inf)*(labels.sites >= 0)
            self['pore.invasion_pressure'][pinds] = inv_val
            # Store current applied pressure in newly invaded throats
            tinds = (self['throat.invasion_pressure'] == sp.inf)*(labels.bonds >= 0)
            self['throat.invasion_pressure'][tinds] = inv_val

        # Convert invasion pressures in sequence values
        Pinv = self['pore.invasion_pressure']
        Tinv = self['throat.invasion_pressure']
        Pseq = sp.searchsorted(sp.unique(Pinv), Pinv)
        Tseq = sp.searchsorted(sp.unique(Tinv), Tinv)
        self['pore.invasion_sequence'] = Pseq
        self['throat.invasion_sequence'] = Tseq

    def _apply_percolation(self, inv_val):
        if self.settings['mode'] == 'bond':
            t_invaded = self['throat.entry_pressure'] <= inv_val
            labels = self.bond_percolation(conns, t_invaded)
        elif self.settings['mode'] == 'site':
            p_invaded = self['pore.entry_pressure'] <= inv_val
            labels = self.site_percolation(conns, p_invaded)
        elif self.settings['mode'] == 'mixed':
            raise Exception('Mixed percolation is not implemented yet')
        else:
            raise Exception('Unrecognized percolation process specified')

    def remove_isolated_clusters(self, labels):
        # Identify clusters on invasion sites
        inv_clusters = sp.unique(labels.sites[self['pore.inlets']])
        # Remove cluster numbers == -1, if any
        inv_clusters = inv_clusters[inv_clusters >= 0]
        # Find all pores in invading clusters
        p_invading = np.in1d(labels.sites, inv_clusters)
        labels.sites[~p_invading] = -1
        t_invading = np.in1d(labels.bonds, inv_clusters)
        labels.bonds[~t_invading] = -1
        return labels

    def site_percolation(self, ij, p_invaded):
        r"""
        """
        Np = p_invaded.size
        t_invaded = sp.all(p_invaded[ij], axis=1)
        adj_mat = sprs.csr_matrix((t_invaded, (ij[:, 0], ij[:, 1])),
                                  shape=(Np, Np))
        adj_mat.eliminate_zeros()
        clusters = csg.connected_components(csgraph=adj_mat, directed=False)[1]
        valid_clusters = sp.bincount(clusters) > 1
        mapping = -sp.ones(shape=(clusters.max()+1, ), dtype=int)
        mapping[valid_clusters] = sp.arange(0, valid_clusters.sum())
        p_labels = mapping[clusters]
        t_labels = sp.amin(p_labels[ij], axis=1)
        result = tup(p_labels, t_labels)
        return result

    def bond_percolation(self, ij, t_invaded):
        r"""
        """
        Np = sp.amax(ij) + 1
        adj_mat = sprs.csr_matrix((t_invaded, (ij[:, 0], ij[:, 1])),
                                  shape=(Np, Np))
        adj_mat.eliminate_zeros()
        clusters = csg.connected_components(csgraph=adj_mat, directed=False)[1]
        valid_clusters = sp.bincount(clusters) > 1
        mapping = -sp.ones(shape=(clusters.max()+1, ), dtype=int)
        mapping[valid_clusters] = sp.arange(0, valid_clusters.sum())
        p_labels = mapping[clusters]
        t_labels = sp.amin(p_labels[ij], axis=1)
        result = tup(p_labels, t_labels)
        return result

    def results(self, Pc):
        r"""
        This method determines which pores and throats are filled with non-
        wetting phase at the specified capillary pressure, and creates or
        updates 'pore.occupancy' and 'throat.occupancy' arrays on the
        associated Phase objects. Invasion pressure and sequence are also sent
        to the invading phase.

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
        proj = self.project
        net = proj.network
        Psatn = self['pore.invasion_pressure'] <= Pc
        Tsatn = self['throat.invasion_pressure'] <= Pc
        inv_phase = {}
        inv_phase['pore.occupancy'] = sp.array(Psatn, dtype=float)
        inv_phase['throat.occupancy'] = sp.array(Tsatn, dtype=float)
        if self.settings['pore_filling']:
            Vp = self._calc_fractional_filling(element='pore', pressure=Pc)
            Sp = Vp/net[self.settings['pore_volume']]
            zero_ps = net[self.settings['pore_volume']] == 0.0
            Sp[zero_ps] = inv_phase['pore.occupancy'][zero_ps]
            inv_phase['pore.occupancy'] = Sp
        if self.settings['throat_filling']:
            Vt = self._calc_fractional_filling(element='throat', pressure=Pc)
            St = Vt/net[self.settings['throat_volume']]
            zero_ts = net[self.settings['throat_volume']] == 0.0
            St[zero_ts] = inv_phase['throat.occupancy'][zero_ts]
            inv_phase['throat.occupancy'] = St
        return inv_phase
