#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: CEF PNM Team
# License: TBD
# Copyright (c) 2012

#from __future__ import print_function

"""
module __GenericGenerator__: Base class to construct pore networks
==================================================================

.. warning:: The classes of this module should be loaded through the 'Generators.__init__.py' file.

"""

import OpenPNM
import scipy as sp
import numpy as np
import scipy.sparse as sprs
import scipy.stats as spst

class GenericGenerator(OpenPNM.Base.OpenPNMbase):
    r"""
    GenericGenerator - Base class to construct pore networks
    
    This class contains the interface definition for the construction
    of networks
    
    Parameters
    ----------
    loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)
    
    
        
    Attributes
    ----------
        net :   OpenPNM.Network.GenericNetwork
            Initialized network
        
            
    Examples
    --------
    
    To reserve space for a network with the default number of pores
    and throats execute
    
    >>> import OpenPNM as PNM
    >>> net=PNM.Generators.GenericGenerator().generate()
    >>> OpenPNM.IO.NetToVtp(net)   
    This generates the following network:
    
    """
    pore_properties   = {}
    throat_properties = {}
    _net              = None
    def __init__(self,  psd_dist = 'gamma',
                        psd_shape = 2,
                        psd_loc = 0,
                        psd_scale = 1,
                        tsd_dist = 'gamma',
                        tsd_shape = 2,
                        tsd_loc = 0,
                        tsd_scale = 1,
                        btype = [0,0,0],
                        **kwargs):
        
        r"""
        Initialize
        """
        super(GenericGenerator,self).__init__(**kwargs)
        self._logger.debug("Method: Constructor")
        
        #Set statistical distribution info
        self._psd_dist = psd_dist
        self._psd_shape = psd_shape
        self._psd_loc = psd_loc
        self._psd_scale = psd_scale
        self._tsd_dist = tsd_dist
        self._tsd_shape = tsd_shape
        self._tsd_loc = tsd_loc
        self._tsd_scale = tsd_scale
        self._btype = btype
        
    def generate(self):
        r"""
        Generate the network
        """
        self._logger.debug("self.generate()")
        self.generate_pores()
        self.generate_throats()
        self._net.update()
        self.generate_boundary()
        #self.add_boundaries()
        self.generate_pore_seeds()
        self.generate_throat_seeds()
        self.generate_pore_diameters()
        self.generate_throat_diameters()
        self.calc_pore_volumes()
        self.calc_throat_lengths()
        self.calc_throat_volumes()
        self._net.update()
        self._logger.debug("\t end of self.generate()")
        return self._net  
        
    def generate_pores(self):
        r"""
        Generate the pores (numbering, types and coordinates)
        """
        self._logger.error("generate_pores: not implemented")
        
    def generate_throats(self):        
        r"""
        Generate the throats (numbering and types)
        """
        self._logger.error("generate_throats: not implemented")
        
    def generate_boundary(self):
        r"""
        This calls generate_throats and pores once again to create Nx*Ny*1...Nx*1*Nz....1*Ny*Nz sized layers.
        A script for translating this boundary in space will modify all of the parameters that are calculated with the pores and throats are generated. 
        An additional script for stitching will append the larger pore network to include these extra layers after they have been translated.
        """
        self._logger.error("generate_bounaries: not implemented")
        
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
    test=GenericGenerator(loggername="TestGenerator")        

