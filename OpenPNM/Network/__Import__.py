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
import numpy as np
import scipy.sparse as sprs
import scipy.stats as spst
from .__GenericNetwork__ import GenericNetwork
import os

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
        self._mat=OpenPNM.Network.ImportMat(filename=filename,path=path)
        self._Np=np.size(self._mat.getvar('pnumbering'))
        self._Nt=np.size(self._mat.getvar('tnumbering'))
        self._net=OpenPNM.Network.GenericNetwork(name=self.name,num_pores=self._Np, num_throats=self._Nt)
        Pind = np.arange(0,self._Np)
        self._net.set_pore_info(label='all',locations=np.ones_like(Pind))
        self._net.set_pore_data(prop='numbering',phase='',data=Pind)
        Tind = np.arange(0,self._Nt)
        self._net.set_throat_info(label='all',locations=np.ones_like(Tind))
        self._net.set_throat_data(prop='numbering',data=Tind)        
        self._logger.info('Writing pore properties')
        self._logger.debug('Writing pore volumes')
        self._net.set_pore_data(prop='volume',data=np.reshape(self._mat.getvar('pvolume'),(self._Np)))        
        self._logger.debug('Writing pore seeds')
        self._net.set_pore_data(prop='seed',data=np.zeros((self._Np),dtype=np.float))
        self._logger.debug('Writing pore diameters')
        self._net.set_pore_data(prop='diameter',data=np.reshape(self._mat.getvar('pdiameter'),(self._Np)))
        self._logger.debug('Writing pore numbering')
        self._net.set_pore_data(prop='numbering',data=np.reshape(self._mat.getvar('pnumbering'),(self._Np)))
        self._logger.debug('Writing pore type')
        self._net.set_pore_data(prop='type',data=np.reshape(self._mat.getvar('ptype'),(self._Np)))
        self._logger.debug('Writing pore coordinates')
        self._net.set_pore_data(prop='coords',data=self._mat.getvar('pcoords'))
        
        self._logger.info('Writing throat properties')
        self._logger.debug('Writing throat lengths')
        self._net.set_throat_data(prop='length',data=np.zeros((self._Nt),dtype=np.float))
        self._logger.debug('Writing throat seeds')
        self._net.set_throat_data(prop='seed',data=np.zeros((self._Nt),dtype=np.float))
        self._logger.debug('Writing throat diameters')
        self._net.set_throat_data(prop='diameter',data=np.reshape(self._mat.getvar('tdiameter'),(self._Nt)))
        self._logger.debug('Writing throat numbering')
        self._net.set_throat_data(prop='numbering',data=np.reshape(self._mat.getvar('tnumbering'),(self._Nt)))
        self._logger.debug('Writing throat type')
        self._net.set_throat_data(prop='type',data=np.reshape(self._mat.getvar('ttype'),(self._Nt)))
        self._logger.debug('Writing throat volumes')
        self._net.set_throat_data(prop='volume',data=np.zeros((self._Nt),dtype=np.float))
        self._logger.debug('Writing throat connections')
        self._net.set_throat_data(prop='connections',data=self._mat.getvar('tconnections'))
 
        return self._net  