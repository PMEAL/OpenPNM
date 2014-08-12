# -*- coding: utf-8 -*-
# Author: CEF PNM Team
# License: TBD
# Copyright (c) 2012

#from __future__ import print_function

"""
module __GenericGeometry__: Base class to construct pore networks
==================================================================

.. warning:: The classes of this module should be loaded through the 'Geometry.__init__.py' file.

"""

import OpenPNM
import scipy as sp
import scipy.io as spio
import os
from .__GenericNetwork__ import GenericNetwork

class MatFile(GenericNetwork):
    r"""
    MatFile - constructs a pore network from a perfectly formatted .mat file (MATLAB)
    
    This class contains the interface definition for the construction
    of networks
    
    Parameters
    ----------

    loglevel : int
        Level of the logger (10=Debug, 20=Info, 30=Warning, 40=Error, 50=Critical)
    all other parameters are in the generate() command

    
    """
    def __init__(self, **kwargs):
        
        r"""
        Initialize
        """
        super(MatFile,self).__init__(**kwargs)
        
    def generate(self,filename='standard_cubic_5x5x5.mat', path='LocalFiles', xtra_pore_data=None, xtra_throat_data=None):
        '''
        Create network from Matlab file. Returns OpenPNM.Network.GenericNetwork() object.

        Parameters
        ----------

        Critical\n
        filename : string
            filename = 'standard_cubic_5x5x5.mat' (default)\n
            Name of mat file\n
        path : string
            path='LocalFiles' (default)\n
            the location of the mat file on your computer \n
        xtra_pore_data : list of strings
            xtra_pore_data = ['type','shape','material']
            any additional props to look for in the dictionary
        xtra_throat_data : list of strings
            xtra_throat_data = ['type','shape','material']
            any additional props to look for in the dictionary

        Examples:
        ---------

        generate network using example mat file

        >>> import OpenPNM as PNM
        >>> pn=PNM.Geometry.MatFile(filename='standard_cubic_5x5x5.mat', path='LocalFiles')
        '''
        if path == 'LocalFiles':
            long_path = os.path.abspath(__file__)
            short_path, fname = os.path.split(long_path)
            short_path, foldername = os.path.split(short_path)  
            path, foldername = os.path.split(short_path)  
            path = os.path.join(path,'LocalFiles')
        self._path = path
        filepath = os.path.join(self._path,filename)
        self._xtra_pore_data=xtra_pore_data
        self._xtra_throat_data=xtra_throat_data
        self._dictionary=spio.loadmat(filepath)
        
        self._Np=sp.size(self._dictionary['pnumbering'])
        self._Nt=sp.size(self._dictionary['tnumbering'])
        
        #Run through generation steps
        self._add_pores()
        self._add_throats()
        self._add_geometry()
        self._add_xtra_pore_data()
        self._add_xtra_throat_data()
        
    def _add_pores(self):
        Pind = sp.arange(0,self._Np)
        self['pore.all'] = sp.ones_like(Pind,dtype=bool)
        self._logger.info('Writing pore data')
        self['pore.coords']=self._dictionary['pcoords']
        
    def _add_throats(self):
        Tind = sp.arange(0,self._Nt)
        self.set_throat_info(label='all',locations=sp.ones_like(Tind))
        self._logger.info('Writing throat data')
        self['throat.conns']=self._dictionary['tconnections']
        
    def _add_geometry(self):
        
        geom = OpenPNM.Geometry.GenericGeometry(network=self,name='imported')
        geom.set_locations(pores='all',throats='all')
        
        data = self._dictionary['pvolume']
        geom.add_property(prop='pore_volume',model='constant',value=data)
        data = self._dictionary['pdiameter']
        geom.add_method(prop='pore_diameter',model='constant',value=data)
        
        data = self._dictionary['tdiameter']
        geom.add_property(prop='throat_diameter',model='constant',value=data)
        
        geom.add_property(prop='throat_length',model='straight')
        geom.regenerate()
    
    def _add_xtra_pore_data(self):
        xpdata = self._xtra_pore_data
        if xpdata is not None:
            if type(xpdata) is type([]):
                for pdata in xpdata:
                    try:
                        self['pore.'+pdata]=self._dictionary['p'+pdata])
                    except:
                        self._logger.warning('Could not add pore data: '+pdata+' to network')
            else:
                try:
                    self['pore.'+xpdata]=self._dictionary['p'+xpdata]
                except:
                    self._logger.warning('Could not add pore data: '+xpdata+' to network')

    def _add_xtra_throat_data(self):
        xtdata = self._xtra_throat_data
        if xtdata is not None:
            if type(xtdata) is type([]):
                for tdata in xtdata:
                    try:
                        self['throat.'+tdata]=self._dictionary['t'+tdata]
                    except:
                        self._logger.warning('Could not add throat data: '+tdata+' to network')
            else:
                try:
                    self['throat.'+xtdata]=self._dictionary['t'+xtdata]
                except:
                    self._logger.warning('Could not add throat data: '+xtdata+' to network')