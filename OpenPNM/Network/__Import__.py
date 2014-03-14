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
        
    def generate(self,filename='standard_cubic_5x5x5.mat', path='LocalFiles'):
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
        self._dictionary=spio.loadmat(filepath)
        
        self._Np=sp.size(self._dictionary['pnumbering'])
        self._Nt=sp.size(self._dictionary['tnumbering'])
        
        #Run through generation steps
        self._add_pores()
        self._add_throats()
        self._add_geometry()
 
        return self
        
    def _add_pores(self):
        Pind = sp.arange(0,self._Np)
        self.set_pore_info(label='all',locations=sp.ones_like(Pind))
        self._logger.info('Writing pore data')
        self.set_pore_data(prop='coords',data=self._dictionary['pcoords'])
        
    def _add_throats(self):
        Tind = sp.arange(0,self._Nt)
        self.set_throat_info(label='all',locations=sp.ones_like(Tind))
        self._logger.info('Writing throat data')
        self.set_throat_data(prop='connections',data=self._dictionary['tconnections'])
        
    def _add_geometry(self):
        
        geom = OpenPNM.Geometry.GenericGeometry(network=self,name='imported')
        
        data = self._dictionary['pvolume']
        geom.add_method(prop='pore_volume',model='constant',value=data)
        data = self._dictionary['pdiameter']
        geom.add_method(prop='pore_diameter',model='constant',value=data)
        
        data = self._dictionary['tdiameter']
        geom.add_method(prop='throat_diameter',model='constant',value=data)
    
