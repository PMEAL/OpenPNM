# -*- coding: utf-8 -*-
"""
===============================================================================
MatFile: Subclass to import Networks from a Matlab '.mat' file
===============================================================================

"""
import scipy as sp
import scipy.io as spio
import os
from OpenPNM.Network import GenericNetwork
import OpenPNM.Geometry
import OpenPNM.Utilities.misc as misc
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)


class MatFile(GenericNetwork):
    r"""
    MatFile - constructs a pore network from a perfectly formatted .mat file (MATLAB)
    Create network from Matlab file. Returns OpenPNM.Network.GenericNetwork()
    object. The extra data of 'type' will trigger internal and boundary pores.

    Parameters
    ----------
    filename : string
        filename = 'standard_cubic_5x5x5.mat' (default)
        Name of mat file
    path : string
        path='' (default)
        the full path to the mat file on your computer
        leaving blank searches for the file in the local directory
    xtra_pore_data : list of strings
        xtra_pore_data = ['type', 'shape', 'material']
        any additional props to look for in the dictionary
    xtra_throat_data : list of strings
        xtra_throat_data = ['type', 'shape', 'material']
        any additional props to look for in the dictionary

    Notes
    ------
    Matfiles should include the following variables

    +----------------+------------+----------------------------------+
    | Variable Name  | Value      | Description                      |
    +================+============+==================================+
    | pcoords        | <Npx3>     | physical coordinates, in meters, |
    |                | float      | of pores to be imported          |
    +----------------+------------+----------------------------------+
    | pdiameter      | <Npx1>     | pore diamters, in meters         |
    |                | float      |                                  |
    +----------------+------------+----------------------------------+
    | pvolume        | <Npx1>     | pore volumes, in cubic meters    |
    |                | float      |                                  |
    +----------------+------------+----------------------------------+
    | pnumbering     | <Npx1>     | = 0:1:Np-1                       |
    |                | int        |                                  |
    +----------------+------------+----------------------------------+
    | ptype          | <Npx1>     | (optional) designates surfaces   |
    |                | int        | of pores in network.             |
    +----------------+------------+----------------------------------+
    | tconnections   | <Ntx2>     | pore numbers of the two pores    |
    |                | int        | that each throat connects        |
    +----------------+------------+----------------------------------+
    | tdiameter      | <Ntx1>     | throat diameters, in meters      |
    |                | float      |                                  |
    +----------------+------------+----------------------------------+
    | tnumbering     | <Ntx1>     | = 0:1:Nt-1                       |
    |                | int        |                                  |
    +----------------+------------+----------------------------------+
    | ttype          | <Ntx1>     | (optional) designates surface    |
    |                | int        | throats in network.              |
    +----------------+------------+----------------------------------+
    """

    def __init__(self, filename='', path='', xtra_pore_data=None,
                 xtra_throat_data=None, **kwargs):

        super().__init__(**kwargs)
        if filename == '':
            return
        if path == '':
            path = os.path.abspath('.')
        self._path = path
        filepath = os.path.join(self._path, filename)
        self._xtra_pore_data = xtra_pore_data
        self._xtra_throat_data = xtra_throat_data
        self._dictionary = spio.loadmat(filepath)

        self._Np = sp.size(self._dictionary['pnumbering'])
        self._Nt = sp.size(self._dictionary['tnumbering'])

        # Run through generation steps
        self._add_pores()
        self._add_throats()
        self._remove_disconnected_clusters()
        self._add_xtra_pore_data()
        self._add_xtra_throat_data()
        self._add_geometry()

    def _add_pores(self):
        Pind = sp.arange(0, self._Np)
        self['pore.all'] = sp.ones_like(Pind, dtype=bool)
        logger.info('Writing pore data')
        self['pore.coords'] = sp.array(self._dictionary['pcoords'], float)

    def _add_throats(self):
        Tind = sp.arange(0, self._Nt)
        self['throat.all'] = sp.ones_like(Tind, dtype=bool)
        logger.info('Writing throat data')
        self['throat.conns'] = sp.array(self._dictionary['tconnections'], int)

    def _remove_disconnected_clusters(self):
        bad_pores = sp.array([], dtype=int)
        self._pore_map = self.pores()
        self._throat_map = self.throats()
        health = self.check_network_health()
        if health['disconnected_clusters'] == []:
            self._throat_map = self.throats()
            self._pore_map = self.pores()
        else:
            Np = self.num_pores()
            Nt = self.num_throats()
            cluster_sizes = [sp.shape(x)[0] for x in health['disconnected_clusters']]
            # 50 or less, if it's a really small network.
            acceptable_size = min([min([50, Np/2]), max(cluster_sizes)])
            # Step through each cluster of pores. If its a small cluster,
            # add it to the list
            for cluster in health['disconnected_clusters']:
                if sp.shape(cluster)[0] < acceptable_size:
                    bad_pores = sp.append(bad_pores, sp.ravel(cluster))
            bad_throats = sp.unique(self.find_neighbor_throats(bad_pores))
            # Create map for pores
            if sp.shape(bad_pores)[0] > 0:
                i = 0
                self._pore_map = sp.zeros((Np-sp.shape(bad_pores)[0],), dtype=int)
                for pore in self.pores():
                    if pore not in bad_pores:
                        self._pore_map[i] = pore
                        i += 1
            # Create map for throats
            if sp.shape(bad_throats)[0] > 0:
                i = 0
                self._throat_map = sp.zeros((Nt - sp.shape(bad_throats)[0],),
                                            dtype=int)
                for throat in self.throats():
                    if throat not in bad_throats:
                        self._throat_map[i] = throat
                        i += 1
            # Fix the pore transformer
            try:
                if sp.shape(bad_pores)[0] > 0:
                    i = 0
                    old_transform = self._dictionary['pname_transform']
                    self._dictionary['pname_transform'] = \
                        sp.zeros((Np-sp.shape(bad_pores)[0],), dtype=int)
                    for pore in self.pores():
                        if pore not in bad_pores:
                            self._dictionary['pname_transform'][i] = \
                                old_transform[pore]
                            i += 1
            except:
                logger.info('Could not update pname_transform. Imported network \
                             may not have had it.')
                pass
            self.trim(pores=bad_pores)

    def _add_geometry(self):
        try:
            boundary_pores = sp.where(self['pore.type'] != 0)[0]
            boundary_throats = sp.where(self['throat.type'] != 0)[0]
            self['throat.top'] = sp.ravel(self['throat.type'] == 1)
            self['throat.bottom'] = sp.ravel(self['throat.type'] == 6)
            self['throat.left'] = sp.ravel(self['throat.type'] == 2)
            self['throat.right'] = sp.ravel(self['throat.type'] == 5)
            self['throat.front'] = sp.ravel(self['throat.type'] == 3)
            self['throat.back'] = sp.ravel(self['throat.type'] == 4)
            self['pore.top'] = self.tomask(pores=sp.ravel(
                self.find_connected_pores(self.throats('top'))))
            self['pore.bottom'] = self.tomask(sp.ravel(
                self.find_connected_pores(self.throats('bottom'))))
            self['pore.left'] = self.tomask(sp.ravel(
                self.find_connected_pores(self.throats('left'))))
            self['pore.right'] = self.tomask(sp.ravel(
                self.find_connected_pores(self.throats('right'))))
            self['pore.front'] = self.tomask(sp.ravel(
                self.find_connected_pores(self.throats('front'))))
            self['pore.back'] = self.tomask(sp.ravel(
                self.find_connected_pores(self.throats('back'))))
            add_boundaries = True
        except:
            boundary_pores = sp.array([])
            boundary_throats = sp.array([])
            logger.info('No boundary pores added.')
            add_boundaries = False
        Ps = sp.where([pore not in boundary_pores for pore in self.pores()])[0]
        Ts = sp.where(
            [throat not in boundary_throats for throat in self.throats()])[0]
        geom = OpenPNM.Geometry.GenericGeometry(network=self, pores=Ps, throats=Ts)
        geom['pore.volume'] = sp.ravel(sp.array(
            self._dictionary['pvolume'][self._pore_map[Ps]], float))
        geom['pore.diameter'] = sp.ravel(sp.array(
            self._dictionary['pdiameter'][self._pore_map[Ps]], float))
        geom['throat.diameter'] = sp.ravel(sp.array(
            self._dictionary['tdiameter'][self._throat_map[Ts]], float))
        geom.add_model(propname='pore.area',
                       model=OpenPNM.Geometry.models.pore_area.spherical)
        geom.add_model(propname='throat.area',
                       model=OpenPNM.Geometry.models.throat_area.cylinder)

        if add_boundaries:
            boun = OpenPNM.Geometry.Boundary(network=self,
                                             pores=boundary_pores,
                                             throats=boundary_throats,
                                             name='boundary')
            self['pore.top_boundary'] = \
                self.tomask(pores=self.pores(['top', 'boundary'],
                            mode='intersection'))
            self['pore.bottom_boundary'] = \
                self.tomask(pores=self.pores(['bottom', 'boundary'],
                            mode='intersection'))
            self['pore.left_boundary'] = \
                self.tomask(pores=self.pores(['left', 'boundary'],
                            mode='intersection'))
            self['pore.right_boundary'] = \
                self.tomask(pores=self.pores(['right', 'boundary'],
                            mode='intersection'))
            self['pore.front_boundary'] = \
                self.tomask(pores=self.pores(['front', 'boundary'],
                            mode='intersection'))
            self['pore.back_boundary'] = \
                self.tomask(pores=self.pores(['back', 'boundary'],
                            mode='intersection'))

            self['throat.top_boundary'] = \
                self.tomask(throats=self.throats(['top', 'boundary'],
                            mode='intersection'))
            self['throat.bottom_boundary'] = \
                self.tomask(throats=self.throats(['bottom', 'boundary'],
                            mode='intersection'))
            self['throat.left_boundary'] = \
                self.tomask(throats=self.throats(['left', 'boundary'],
                            mode='intersection'))
            self['throat.right_boundary'] = \
                self.tomask(throats=self.throats(['right', 'boundary'],
                            mode='intersection'))
            self['throat.front_boundary'] = \
                self.tomask(throats=self.throats(['front', 'boundary'],
                            mode='intersection'))
            self['throat.back_boundary'] = \
                self.tomask(throats=self.throats(['back', 'boundary'],
                            mode='intersection'))

    def _add_xtra_pore_data(self):
        xpdata = self._xtra_pore_data
        if xpdata is not None:
            if isinstance(xpdata, type([])):
                for pdata in xpdata:
                    try:
                        self['pore.' + pdata] = \
                            self._dictionary['p' + pdata][self._pore_map]
                    except:
                        logger.warning('Could not add pore data: ' + pdata +
                                       ' to network')
                        pass
            else:
                try:
                    self['pore.' + xpdata] = \
                        self._dictionary['p'+xpdata][self._pore_map]
                except:
                    logger.warning('Could not add pore data: '+xpdata+' to network')
                    pass

    def _add_xtra_throat_data(self):
        xtdata = self._xtra_throat_data
        if xtdata is not None:
            if isinstance(xtdata, type([])):
                for tdata in xtdata:
                    try:
                        self['throat.' + tdata] = \
                            self._dictionary['t' + tdata][self._throat_map]
                    except:
                        logger.warning('Could not add throat data: ' + tdata +
                                       ' to network')
                        pass
            else:
                try:
                    self['throat.' + xtdata] = \
                        self._dictionary['t' + xtdata][self._throat_map]
                except:
                    logger.warning('Could not add throat data: ' + xtdata +
                                   ' to network')
                    pass
