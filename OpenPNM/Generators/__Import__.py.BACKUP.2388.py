# -*- coding: utf-8 -*-
<<<<<<< HEAD

"""
module __Import__: Import networks from file
==========================================================
=======
# Author: CEF PNM Team
# License: TBD
# Copyright (c) 2012

#from __future__ import print_function

"""
module __GenericGenerator__: Base class to construct pore networks
==================================================================
>>>>>>> develop

.. warning:: The classes of this module should be loaded through the 'Generators.__init__.py' file.

"""

import OpenPNM
import scipy as sp
import numpy as np
<<<<<<< HEAD
from time import clock

=======
import scipy.sparse as sprs
import scipy.stats as spst
>>>>>>> develop
from __GenericGenerator__ import GenericGenerator

class MatFile(GenericGenerator):
    r"""
<<<<<<< HEAD
    Cubic - Class to create a basic cubic network
    
    Parameters
    ----------
    
    domain_size : list with 3 float elements 
        Shape of the cube [Lx,Ly,Lz]
    lattice_spacing : list of three floats
        Spacing between pore centers in each spatial directions
    loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)
        
    Examples
    --------
    
    .. plot::
        
       import pylab as pl
       import OpenPNM
       gen = OpenPNM.Generators.Cubic()
       net = gen.generate()
       pl.spy(net._adjmatrix)
       pl.show()
    
    TODO:
        - Check for 3D shape
        - Check for correct divisions
    """
    
    def __init__(self,  filename = [],path = [],**kwargs):

        super(MatFile,self).__init__(**kwargs)
        self._logger.debug("Execute constructor")
        self._logger.info("Reading .mat file")
        self._mat = OpenPNM.IO.ImportMat(filename=filename,path=path)
        self._Np = np.size(self._mat.getvar('pnumbering'))
        self._Nt = np.size(self._mat.getvar('tnumbering'))
        
        self._logger.info("Generate network object")        
        #Instantiate object(correct terminology?)
        self._net=OpenPNM.Network.GenericNetwork(num_pores=self._Np, num_throats=self._Nt)
        
        
=======
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
    >>> net=PNM.Generators.MatFile(filename='example_network', path='D:\\AFCC code\\GitHub projects\\OpenPNM\\LocalFiles').generate()
    
    """
    def __init__(self, filename='example_network', path='D:\\AFCC code\\GitHub projects\\OpenPNM\\LocalFiles', **kwargs):
        
        r"""
        Initialize
        """
        super(MatFile,self).__init__(**kwargs)
        self._mat=OpenPNM.IO.ImportMat(filename=filename,path=path)
        self._Np=np.size(self._mat.getvar('pnumbering'))
        self._Nt=np.size(self._mat.getvar('tnumbering'))
        self._net=OpenPNM.Network.GenericNetwork(num_pores=self._Np, num_throats=self._Nt)
        
>>>>>>> develop
    def generate(self):
        r"""
        Generate the network
        """
<<<<<<< HEAD
        self._logger.info("Importing Pores")
        self._net.pore_properties['coords'] = np.float64(self._mat.getvar('pcoords'))
        self._net.pore_properties['numbering'] = np.reshape(np.int32(self._mat.getvar('pnumbering')),self._Np)
        self._net.pore_properties['type']= np.reshape(np.int8(self._mat.getvar('ptype')),self._Np)
        self._net.pore_properties['volume']= np.reshape(np.float64(self._mat.getvar('pvolume')),self._Np)
        self._net.pore_properties['diameter']= np.reshape(np.float64(self._mat.getvar('pdiameter')),self._Np)
        
        self._logger.info("Importing Throats")
        self._net.throat_properties['connections'] = np.int32(self._mat.getvar('tconnections'))
        self._net.throat_properties['numbering'] = np.reshape(np.int32(self._mat.getvar('tnumbering')),self._Nt)
        self._net.throat_properties['type'] = np.reshape(np.int8(self._mat.getvar('ttype')),self._Nt)
        self._net.throat_properties['diameter'] = np.reshape(np.float64(self._mat.getvar('tdiameter')),self._Nt)
        self._net.throat_properties['volume'] = np.zeros(self._Nt,dtype=np.float64)
        self._net.throat_properties['length'] = np.zeros(self._Nt,dtype=np.float64)
        
        #self.generate_pore_seeds()
        #self.generate_throat_seeds()
        #self.calc_throat_lengths()
        #self.calc_throat_volumes()
        #self._net.update()
        self._logger.debug("\t end of self.generate()")
        return self._net  
        
    def generate_pores(self):
        r"""
        Generate the pores (numbering, types and coordinates)
        """
        
    def generate_throats(self):        
        r"""
        Generate the throats (numbering and types)
        """
        
    def generate_throat_type(self):
        self._logger.info("update_throat_type: Start of method")
        self._net.throat_properties["type"] = self._net._incmatrix.transpose() * self._net.pore_properties["type"]
        self._logger.debug("update_throat_type: End of method")

    def generate_pore_seeds(self):
        r"""
        Assigns random seed to pores
        """
        self._logger.info("generate_pore_seeds: Assign each pore a random seed")
        Np = self._net.get_num_pores()
        self._net.pore_properties['seed'] = np.random.rand(Np)
        #Set boundary pore to 1 (max possible value) so throats adopt seed from interal pore
        self._net.pore_properties['seed'][self._net.pore_properties['type']>0] = 1
        self._logger.debug("generate_pore_seeds: End of method")

    def generate_throat_seeds(self):
        r"""
        Assigns random seeds to throats based on smaller of neighboring pores
        """
        self._logger.info("generate_throat_seeds: Assign each throat its smaller neighboring pore seed")
        self._net.throat_properties['seed'] = np.amin(self._net.pore_properties['seed'][self._net.throat_properties['connections']],1)
        self._logger.debug("generate_throat_seeds: End of method")
        
    def generate_pore_diameters(self):
        r"""
        Calculates pore diameter from given statisical distribution
        """
        self._logger.info("generate_pore_diameters: Generate pore diameter from "+self._psd_dist+" distribution")
        prob_fn = getattr(spst,self._psd_dist)
        P = prob_fn(self._psd_shape,loc=self._psd_loc,scale=self._psd_scale)
        self._net.pore_properties['diameter'] = P.ppf(self._net.pore_properties['seed'])
        #Set boundadry pores to size 0
        self._net.pore_properties['diameter'][self._net.pore_properties['type']>0] = 0
        self._logger.debug("generate_pore_diameters: End of method")

    def generate_throat_diameters(self):
        r"""
        Calculates throat diameter from given statisical distribution
        """
        self._logger.info("generate_throat_diameters: Generate throat diameter from "+self._tsd_dist+" distribution")
        prob_fn = getattr(spst,self._tsd_dist)
        P = prob_fn(self._tsd_shape,loc=self._tsd_loc,scale=self._tsd_scale)
        self._net.throat_properties['diameter'] = P.ppf(self._net.throat_properties['seed'])
        self._logger.debug("generate_throat_diameters: End of method")
        
    def calc_pore_volumes(self):
        r"""
        Calculates pore volume
        """
        self._logger.info("calc_pore_volumes: Setting pore volumes assuming cubic bodies")
        #Set internal pore volumes to 1
        self._net.pore_properties['volume'] = self._net.pore_properties['diameter']**3
        #Set boundary pore volumes to 0
        self._net.pore_properties['volume'][self._net.pore_properties['type']>0] = 0
        self._logger.debug("calc_pore_volumes: End of method")
        
    def calc_throat_volumes(self):
        r"""
        Calculates throat volume from diameter and length
        """
        self._logger.info("calc_throat_volumes: Setting throat volumes assuming square cross-section")
        #Set internal pore volumes to 1
        self._net.throat_properties['volume'] = self._net.throat_properties['length']*self._net.throat_properties['diameter']**2
        self._logger.debug("calc_throat_volumes: End of method")
        
    def calc_throat_lengths(self):
        self._logger.info("calc_throat_lengths: Determine throat length from distance between pores")
        #Initialize throat_property['length']
        self._net.throat_properties['length'] = np.zeros_like(self._net.throat_properties['type'])
        #Find length of throats
        Lx = 0.001
        Ly = 0.001
        Lz = 0.0004
        poffset = np.array([[ 0.,  0.,  0.],
                            [ 0.,  0.,  Lz],
                            [ Lx,  0.,  0.],
                            [ 0.,  Ly,  0.],
                            [ 0., -Ly,  0.],
                            [-Lx,  0.,  0.],
                            [ 0.,  0., -Lz]])
        T1 = self._net.throat_properties['type']
        C1 = self._net.pore_properties['coords'][self._net.throat_properties['connections'][:,0]]
        C2 = self._net.pore_properties['coords'][self._net.throat_properties['connections'][:,1]]
        E = np.sqrt(np.sum((C1-C2)**2,axis=1))  #Euclidean distance between pores
        D1 = self._net.pore_properties['diameter'][self._net.throat_properties['connections'][:,0]]
        D2 = self._net.pore_properties['diameter'][self._net.throat_properties['connections'][:,1]]
        self._net.throat_properties['length'] = E - (D1 + D2)/2
        #Perform check for unphysical throat lengths
        if np.sum(self._net.throat_properties['length']<0):
            self._logger.warning("calc_throat_lengths: Some negative throat lengths exist, some pores overlap!")
        self._logger.debug("calc_throat_lengths: End of method")
        
    def add_boundaries(self):
        self._logger.error("add_boundaries: not implemented")
        
    def get_net(self):
        r"""
        Return the produced network
        """
        return self._net

if __name__ == '__main__':
    filename = 'example_network'
    path = 'D:\\AFCC code\\GitHub projects\\OpenPNM\\LocalFiles'
    self=MatFile(filename=filename,path=path,loggername="TestGenerator")
    pn=self.generate()        
    IP = OpenPNM.Algorithms.InvasionPercolation(net=pn,inlets=np.nonzero(pn.pore_properties['type']==1)[0],outlets=np.nonzero(pn.pore_properties['type']==6)[0],loglevel=10)
    IP.run()
    OpenPNM.IO.NetToVtp(net=pn,filename="testimport.vtp",scaling_factor=1000000)
    
    

=======
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
    OpenPNM.IO.NetToVtp(net=pn)
>>>>>>> develop
