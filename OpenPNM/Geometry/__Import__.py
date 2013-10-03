# -*- coding: utf-8 -*-
# Author: CEF PNM Team
# License: TBD
# Copyright (c) 2012

#from __future__ import print_function

"""
module __GenericGenerator__: Base class to construct pore networks
==================================================================

.. warning:: The classes of this module should be loaded through the 'Geometry.__init__.py' file.

"""

import OpenPNM
import scipy as sp
import numpy as np
import scipy.sparse as sprs
import scipy.stats as spst
from __GenericGenerator__ import GenericGenerator

class MatFile(GenericGenerator):
    r"""
    MatFile - constructs a pore network from a perfectly formatted .mat file (MATLAB)
    
    This class contains the interface definition for the construction
    of networks
    
    Parameters
    ----------
    filename: str
        name of input file
    path: str
        path of the input file
    loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)
    
    
        
    Attributes
    ----------
        net :   OpenPNM.Network.GenericNetwork
            Initialized network
        
            
    Examples
    --------
    
    To import the example_network.mat file in your LocalFiles folder
    
    >>> import OpenPNM as PNM
    >>> net=PNM.Geometry.MatFile(filename='example_network', path='D:\\AFCC code\\GitHub projects\\OpenPNM\\LocalFiles').generate()
    
    """
    def __init__(self, filename='example_network', path='D:\\AFCC code\\GitHub projects\\OpenPNM\\LocalFiles', **kwargs):
        
        r"""
        Initialize
        """
        super(MatFile,self).__init__(**kwargs)
        self._mat=OpenPNM.Utilities.ImportMat(filename=filename,path=path)
        self._Np=np.size(self._mat.getvar('pnumbering'))
        self._Nt=np.size(self._mat.getvar('tnumbering'))
        self._net=OpenPNM.Network.GenericNetwork(num_pores=self._Np, num_throats=self._Nt)
        
    def generate(self):
        r"""
        Generate the network
        """
        self._logger.info('Writing pore properties')
        self._logger.debug('Writing pore volumes')
        self._net.pore_properties['volume']=np.reshape(self._mat.getvar('pvolume'),(self._Np))
        self._logger.debug('Writing pore seeds')
        self._net.pore_properties['seed']=np.zeros((self._Np),dtype=np.float)
        self._logger.debug('Writing pore diameters')
        self._net.pore_properties['diameter']=np.reshape(self._mat.getvar('pdiameter'),(self._Np))
        self._logger.debug('Writing pore numbering')
        self._net.pore_properties['numbering']=np.reshape(self._mat.getvar('pnumbering'),(self._Np))
        self._logger.debug('Writing pore type')
        self._net.pore_properties['type']=np.reshape(self._mat.getvar('ptype'),(self._Np))
        self._logger.debug('Writing pore coordinates')
        self._net.pore_properties['coords']=self._mat.getvar('pcoords')
        
        self._logger.info('Writing throat properties')
        self._logger.debug('Writing throat lengths')
        self._net.throat_properties['length']=np.zeros((self._Nt),dtype=np.float)
        self._logger.debug('Writing throat seeds')
        self._net.throat_properties['seed']=np.zeros((self._Nt),dtype=np.float)
        self._logger.debug('Writing throat diameters')
        self._net.throat_properties['diameter']=np.reshape(self._mat.getvar('tdiameter'),(self._Nt))
        self._logger.debug('Writing throat numbering')
        self._net.throat_properties['numbering']=np.reshape(self._mat.getvar('tnumbering'),(self._Nt))
        self._logger.debug('Writing throat type')
        self._net.throat_properties['type']=np.reshape(self._mat.getvar('ttype'),(self._Nt))
        self._logger.debug('Writing throat volumes')
        self._net.throat_properties['volume']=np.zeros((self._Nt),dtype=np.float)
        self._logger.debug('Writing throat connections')
        self._net.throat_properties['connections']=self._mat.getvar('tconnections')
        
        
        #self._logger.debug("self.generate()")
        #self.generate_pores()
        #self.generate_throats()
        #self._net.update()
        #self.add_boundaries()
        #self.generate_pore_seeds()
        #self.generate_throat_seeds()
        #self.generate_pore_diameters()
        #self.generate_throat_diameters()
        #self.calc_pore_volumes()
        #self.calc_throat_lengths()
        #self.calc_throat_volumes()
        #self._net.update()
        #self._logger.debug("\t end of self.generate()")
        return self._net  

if __name__ == '__main__':
    self=MatFile(filename='example_network', path='D:\\AFCC code\\GitHub projects\\OpenPNM\\LocalFiles',loglevel=10)
    pn=self.generate()
    inlets = np.nonzero(pn.pore_properties['type']==1)[0]
    outlets = np.nonzero(pn.pore_properties['type']==6)[0]
    OpenPNM.Algorithms.InvasionPercolation(net=pn,inlets=inlets,outlets=outlets).run()
    OpenPNM.Visualization.NetToVtp(net=pn)