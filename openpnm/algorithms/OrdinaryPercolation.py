import scipy as sp
import numpy as np
import scipy.sparse as sprs
import scipy.sparse.csgraph as csg
from collections import namedtuple
from openpnm.algorithms import GenericAlgorithm
from openpnm.core import logging
logger = logging.getLogger(__name__)
tup = namedtuple('cluster_labels', ('sites', 'bonds'))


class OrdinaryPercolation(GenericAlgorithm):

    def __init__(self, settings={}, **kwargs):
        r"""
        """
        super().__init__(**kwargs)
        settings.update({'trapping': False,
                         'access_limited': True,
                         'kind': 'bond'})

    def setup(self,
              invading_phase,
              defending_phase=None,
              trapping=False,
              access_limited=True,
              kind='bond',
              entry_pressure='throat.capillary_pressure',
              pore_filling=None,
              throat_filling=None,
              pore_volume='pore.volume',
              throat_volume='throat.volume'):
        r"""
        Used to specify necessary arguments to the simulation.  This method is
        useful for resetting the algorithm or applying more explicit control.

        Parameters
        ----------
        invading_phase : OpenPNM Phase object
            The Phase object containing the physical properties of the invading
            fluid.

        defending_phase : OpenPNM Phase object
            The Phase object containing the physical properties of the
            defending fluid.

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
        self['pore.invasion_pressure'] = sp.inf
        self['throat.invasion_pressure'] = sp.inf
        self['pore.trapped'] = sp.inf
        self['throat.trapped'] = sp.inf
        self['pore.inlets'] = False
        self['pore.outlets'] = False
        self['pore.residual'] = False
        self['throat.residual'] = False
        self.settings['defending_phase'] = None
        self.settings.update({'invading_phase': invading_phase.name,
                              'trapping': trapping,
                              'access_limited': access_limited,
                              'pore_filling': pore_filling,
                              'throat_filling': throat_filling,
                              'throat_volume': throat_volume,
                              'pore_volume': pore_volume})
        if defending_phase:
            self.settings['defending_phase'] = defending_phase.name

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
        # Parse inputs and generate list of invasion points if necessary
        if type(points) is int:
            logger.info('Generating list of invasion pressures')
            # Nudge values up and down
            if self.settings['percolation_type'] == 'bond':
                min_p = sp.amin(self['throat.entry_pressure'])*0.99
                max_p = sp.amax(self['throat.entry_pressure'])*1.01
            elif self.settings['percolation_type'] == 'sites':
                min_p = sp.amin(self['pore.entry_pressure'])*0.99
                max_p = sp.amax(self['pore.entry_pressure'])*1.01
            else:
                raise Exception('Percolation type has not been set')
            points = sp.logspace(sp.log10(min_p), sp.log10(max_p), points)

        # Generate curve from points
        for inv_val in points:

            conns = net['throat.conns']
            # Apply one applied pressure and determine invaded pores
            if self.settings['kind'] == 'bond':
                t_invaded = self['throat.invasion_pressure'] < inv_val
                labels = self.bond_percolation(ij=conns, t_invaded)
            elif self.settings['kind'] == 'site':
                p_invaded = self['pore.invasion_pressure'] < inv_val
                labels = self.site_percolation(ij=conns, p_invaded)
            elif self.settings['kind'] == 'mixed':
                raise Exception('Mixed percolation is not implemented yet')
            else:
                raise Exception('Unrecognized percolation process specified')

            # Optionally remove clusters not connected to the inlets
            if self.settings['access_limited']:
                # Identify clusters on invasion sites
                inv_clusters = sp.unique(labels.sites[self['pore.inlets']])
                # Remove cluster numbers == -1, if any
                inv_clusters = inv_clusters[inv_clusters >= 0]
                # Find all pores in invading clusters
                p_invaded = np.in1d(labels.sites, inv_clusters)
                t_invaded = np.in1d(labels.bonds, inv_clusters)
                labels = tup(p_labels, t_labels)

            # Store current applied pressure in newly invaded pores
            pinds = (self['pore.invasion_pressure'] == sp.inf)*(labels.sites >= 0)
            self['pore.inv_Pc'][pinds] = inv_val
            # Store current applied pressure in newly invaded throats
            tinds = (self['throat.inv_Pc'] == sp.inf)*(labels.bonds >= 0)
            self['throat.inv_Pc'][tinds] = inv_val

        # Convert invasion pressures in sequence values
        Pinv = self['pore.invasion_pressure']
        Tinv = self['throat.invasion_pressure']
        Pseq = sp.searchsorted(sp.unique(Pinv), Pinv)
        Tseq = sp.searchsorted(sp.unique(Tinv), Tinv)
        self['pore.invasion_sequence'] = Pseq
        self['throat.invasion_sequence'] = Tseq

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

    def make_contiguous(self, array):
        valid_clusters = sp.bincount(clusters) > 1
        mapping = -sp.ones(shape=(clusters.max()+1, ), dtype=int)
        mapping[valid_clusters] = sp.arange(0, valid_clusters.sum())
        p_labels = mapping[clusters]

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

    def evaluate_trapping(self):
        r"""
        Finds trapped pores and throats after a full ordinary percolation
        simulation has been run.

        Returns
        -------
        It creates arrays called ``pore.trapped`` and ``throat.trapped``, but
        also adjusts the ``pore.invasion_pressure`` and
        ``throat.invasion_pressure`` arrays to set trapped locations to have
        infinite invasion pressure (i.e. not invaded).

        """
        network = self.project.network
        p_outlets = self['pore.outlets']
        self['pore.trapped'] = sp.zeros([self.Np, ], dtype=float)
        self['throat.trapped'] = sp.zeros([self.Nt, ], dtype=float)
        try:
            # Get points used in OP
            inv_points = sp.unique(self['pore.invasion_pressure'])
        except KeyError:
            raise Exception('Orindary percolation has not been run!')
        tind = network.throats()
        conns = network.find_connected_pores(tind)
        for inv_val in inv_points[0:-1]:
            # Find clusters of defender pores
            Pinvaded = self['pore.inv_Pc'] <= inv_val
            Cstate = sp.sum(Pinvaded[conns], axis=1)
            Tinvaded = self['throat.inv_Pc'] <= inv_val
            # 0 = all open, 1=1 pore filled,
            # 2=2 pores filled 3=2 pores + 1 throat filled
            Cstate = Cstate + Tinvaded
            clusters = network.find_clusters(Cstate == 0)
            # Clean up clusters (invaded = -1, defended >=0)
            clusters = clusters * (~Pinvaded) - (Pinvaded)
            # Identify clusters connected to outlet sites
            out_clusters = sp.unique(clusters[p_outlets])
            trapped_pores = ~sp.in1d(clusters, out_clusters)
            trapped_pores[Pinvaded] = False
            if sum(trapped_pores) > 0:
                inds = (self['pore.trapped'] == 0) * trapped_pores
                self['pore.trapped'][inds] = inv_val
                trapped_throats = network.find_neighbor_throats(trapped_pores)
                trapped_throat_array = np.asarray([False] * len(Cstate))
                trapped_throat_array[trapped_throats] = True
                inds = (self['throat.trapped'] == 0) * trapped_throat_array
                self['throat.trapped'][inds] = inv_val
                inds = (self['throat.trapped'] == 0) * (Cstate == 2)
                self['throat.trapped'][inds] = inv_val
        self['pore.inv_Pc'][self['pore.trapped'] > 0] = sp.inf
        self['throat.inv_Pc'][self['throat.trapped'] > 0] = sp.inf

    def _check_trapping(self, inv_val):
        r"""
        Determine which pores and throats are trapped by invading phase.  This
        method is called by ``run`` if 'trapping' is set to True.
        """
        net = self.project.network
        # Generate a list containing boolean values for throat state
        Tinvaded = self['throat.invasion_pressure'] < sp.inf
        # Add residual throats, if any, to list of invaded throats
        Tinvaded = Tinvaded + self['throat.residual']
        # Invert logic to find defending throats
        Tdefended = ~Tinvaded
        [pclusters, tclusters] = net.find_clusters2(mask=Tdefended,
                                                    t_labels=True)
        # See which outlet pores remain uninvaded
        outlets = self['pore.outlets']*(self['pore.invasion_pressure'] == sp.inf)
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
        tinds = net.find_neighbor_throats(pores=pinds, mode='intersection')
        self['throat.trapped'][tinds] = inv_val
        self['throat.entry_pressure'][tinds] = 1000000

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
            The capillary pressure for which an invading phase configuration is
            required.

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
