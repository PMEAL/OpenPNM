#! /usr/bin/env python
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
import scipy.sparse as sprs
import scipy.stats as spst

class GenericGeometry(OpenPNM.Utilities.OpenPNMbase):
    r"""
    GenericGeometry - Base class to construct pore networks
    
    This class contains the interface definition for the construction of networks
    
    Parameters
    ----------
    loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)
    
        
    """

    def __init__(self, **kwargs):
        
        r"""
        Initialize
        """
        super(GenericGeometry,self).__init__(**kwargs)
        self._logger.debug("Method: Constructor")
        
    def generate(self, **params):
        r"""
        Generate the network
        """
        self._logger.debug("self.generate()")
        self._generate_setup(**params)
        self._generate_pores()
        self._generate_throats()
#        self.add_boundaries()
        self._generate_pore_seeds()
        self._generate_throat_seeds()
        self._generate_pore_diameters(params['psd_info'])
        self._generate_throat_diameters(params['tsd_info'])
        self._calc_pore_volumes()
        self._calc_throat_lengths()
        self._calc_throat_volumes()
        self._logger.debug("\t end of self.generate()")
        return self._net  

    def _generate_setup(self,**params):
        r"""
        Perform applicable preliminary checks and calculations required for generation
        
        Notes
        -----
        This method is not implemented in the GenericGeometry and must be sub-classed to produce desired network topology        
        """
        self._logger.error("generation_setup: not implemented")

    def _generate_pores(self):
        r"""
        Generate the pores (numbering, types and coordinates)
        
        Notes
        -----
        This method is not implemented in the GenericGeometry and must be sub-classed to produce desired network topology        
        """
        self._logger.error("generate_pores: not implemented")
        
    def _generate_throats(self):        
        r"""
        Generate the throats (numbering and types)
        
        Notes
        -----
        This method is not implemented in the GenericGeometry and must be sub-classed to produce desired network topology.
        """
        self._logger.error("generate_throats: not implemented")
        
    def _add_boundaries(self):
        r"""
        Add boundary pores around network (numbering and types)
        
        Notes
        -----
        This method is not implemented in the GenericGeometry and must be sub-classed to produce desired network topology.
        """
        self._logger.error("add_boundaries: not implemented")

    def _generate_pore_seeds(self):
        r"""
        Assign random seed to pores
        
        Notes
        -----
        To reproduce an identical network it is necessary to set the seed of the random number generator
        """
        self._logger.info("generate_pore_seeds: Assign each pore a random seed")
        Np = self._net.get_num_pores()
        self._net.pore_properties['seed'] = sp.random.rand(Np)
        #Set boundary pore to 1 (max possible value) so throats adopt seed from interal pore
        self._net.pore_properties['seed'][self._net.pore_properties['type']>0] = 1
        self._logger.debug("generate_pore_seeds: End of method")

    def _generate_throat_seeds(self):
        r"""
        Assigns random seeds to throats based on smaller of neighboring pores
        
        Notes
        -----
        
        """
        self._logger.info("generate_throat_seeds: Assign each throat its smaller neighboring pore seed")
        self._net.throat_properties['seed'] = sp.amin(self._net.pore_properties['seed'][self._net.throat_properties['connections']],1)
        self._logger.debug("generate_throat_seeds: End of method")
        
    def _generate_pore_diameters(self,psd_info):
        r"""
        Calculates pore diameter from given statisical distribution using the random seeds provided by generate_pore_seeds()
        
        Parameters
        ----------
        psd_info : dictionary
            Contains 'key : value' pairs required by chosen distribution
        
        Notes
        -----
                
        """
        self._logger.info("generate_pore_diameters: Generate pore diameter from "+psd_info['name']+" distribution")
        prob_fn = getattr(spst,psd_info['name'])
        P = prob_fn(psd_info['shape'],loc=psd_info['loc'],scale=psd_info['scale'])
        self._net.pore_properties['diameter'] = P.ppf(self._net.pore_properties['seed'])
        #Set boundadry pores to size 0
        self._net.pore_properties['diameter'][self._net.pore_properties['type']>0] = 0
        self._logger.debug("generate_pore_diameters: End of method")

    def _generate_throat_diameters(self,tsd_info):
        r"""
        Calculates throat diameter from given statisical distribution using the random seeds provided by generate_throat_seeds()
        
        Parameters
        ----------
        tsd_info : dictionary
            Contains 'key : value' pairs required by chosen distribution
        
        Notes
        -----
        
        """
        self._logger.info("generate_throat_diameters: Generate throat diameter from "+tsd_info['name']+" distribution")
        prob_fn = getattr(spst,tsd_info['name'])
        P = prob_fn(tsd_info['shape'],loc=tsd_info['loc'],scale=tsd_info['scale'])
        self._net.throat_properties['diameter'] = P.ppf(self._net.throat_properties['seed'])
        self._logger.debug("generate_throat_diameters: End of method")
        
    def _calc_pore_volumes(self):
        r"""
        Calculates pore volume from diameter assuming a spherical pore
        """
        self._logger.info("calc_pore_volumes: Setting pore volumes assuming cubic bodies")
        #Set internal pore volumes to 1
        self._net.pore_properties['volume'] = self._net.pore_properties['diameter']**3
        #Set boundary pore volumes to 0
        self._net.pore_properties['volume'][self._net.pore_properties['type']>0] = 0
        self._logger.debug("calc_pore_volumes: End of method")
        
    def _calc_throat_volumes(self):
        r"""
        Calculates throat volume from diameter and length assuming a cylindrical pore of constant cross-section
        
        Notes
        -----
        
        """
        self._logger.info("calc_throat_volumes: Setting throat volumes assuming square cross-section")
        #Set internal pore volumes to 1
        self._net.throat_properties['volume'] = self._net.throat_properties['length']*self._net.throat_properties['diameter']**2
        self._logger.debug("calc_throat_volumes: End of method")
        
    def _calc_throat_lengths(self):
        r"""
        Determine throat length from distance between pores
        
        Notes
        -----
        
        """
        self._logger.info("calc_throat_lengths: Determine throat length from distance between pores")
        #Initialize throat_property['length']
        self._net.throat_properties['length'] = sp.zeros_like(self._net.throat_properties['type'])
        C1 = self._net.pore_properties['coords'][self._net.throat_properties['connections'][:,0]]
        C2 = self._net.pore_properties['coords'][self._net.throat_properties['connections'][:,1]]
        E = sp.sqrt(sp.sum((C1-C2)**2,axis=1))  #Euclidean distance between pores
        D1 = self._net.pore_properties['diameter'][self._net.throat_properties['connections'][:,0]]
        D2 = self._net.pore_properties['diameter'][self._net.throat_properties['connections'][:,1]]
        self._net.throat_properties['length'] = E - (D1 + D2)/2
        #Perform check for unphysical throat lengths
        if sp.sum(self._net.throat_properties['length']<0):
            self._logger.warning("calc_throat_lengths: Some negative throat lengths exist, some pores overlap!")
        self._logger.debug("calc_throat_lengths: End of method")     
        
    def translate_coordinates(self,net,displacement=[0,0,0]):
        r"""
        Translate pore network coordinates by specified amount
        
        Parameters
        ----------
        net : OpenPNM Network Object
            The network to which translation should be applied
            
        displacement : array_like
            A vector containing the amount to translate in each dimension. [0,0,0] yeilds no translation.
        
        """
        Lc = net.pore_properties['domain_size'][3]
        N = net.pore_properties['divisions']
        ind = N[0]*N[1]*N[2]
        net.pore_properties['coords'] = Lc*(0.5 + displacement + sp.array(sp.unravel_index(sp.arange(0,ind),dims=net.pore_properties['divisions'],order='F')).T)
        
    def scale_coordinates(self,net,scale=[1,1,1]):
        r"""
        Scale pore network coordinates by specified amount
        
        Parameters
        ----------
        net : OpenPNM Network Object
            The network to which translation should be applied
            
        scale : array_like
            A vector containing the amount to scale in each dimension.  [1,1,1] yeilds no scaling.
            
        """
        net.pore_properties['coords'] = net.pore_properties['coords']*scale 

    def stitch(self,net1,net2,edge = 0):
        r"""
        Stitch two networks together
        
        Parameters
        ----------
        net1 : OpenPNM Network Object
            The network that is stiched to
        
        net2 : OpenPNM Network Object
            The network that is stitched
        
        """
        
        #manipulate all pore_properties. This mainly includes checks on the validity of coords, addition of numbering, and concatening conserved properties. 
        net2.pore_properties['numbering']   = len(net1.pore_properties['numbering']) + net2.pore_properties['numbering']
        net1.pore_properties['numbering']   = sp.concatenate((net1.pore_properties['numbering'],net2.pore_properties['numbering']),axis=0)        
        #print net2.pore_properties['numbering']
        
        #The coordinates should have been translated by this point already using the "translate_coordinates" method. 
        # Checks to ensure that the coordinates are properly translated: ...
        
        net1.pore_properties['coords']      = sp.concatenate((net1.pore_properties['coords'],net2.pore_properties['coords']),axis = 0)
        
        # All of these properties should be conserved when we stitch 2 networks. 
        net1.pore_properties['volume']      = sp.concatenate((net1.pore_properties['volume'],net2.pore_properties['volume']),axis = 0)
        net1.pore_properties['seed']        = sp.concatenate((net1.pore_properties['seed'],net2.pore_properties['seed']),axis = 0)
        net1.pore_properties['type']        = sp.concatenate((net1.pore_properties['type'],net2.pore_properties['type']),axis = 0)
        net1.pore_properties['diameter']    = sp.concatenate((net1.pore_properties['diameter'],net2.pore_properties['diameter']),axis = 0)
        
        net1.throat_properties['numbering'] = len(net1.throat_properties['numbering']) + net2.throat_properties['numbering']
        net1.throat_properties['numbering'] = sp.concatenate((net1.throat_properties['numbering'],net2.throat_properties['numbering']),axis=0)
        
        net1.throat_properties['volume']    = sp.concatenate((net1.throat_properties['volume'],net2.throat_properties['volume']),axis=0)
        net1.throat_properties['diameter']  = sp.concatenate((net1.throat_properties['diameter'],net2.throat_properties['diameter']),axis=0)
        net1.throat_properties['length']    = sp.concatenate((net1.throat_properties['length'],net2.throat_properties['length']),axis=0)
        net1.throat_properties['seed']      = sp.concatenate((net1.throat_properties['seed'],net2.throat_properties['seed']),axis=0)
        
        # We can define the type of the stitched network (net2) to have a value corresponding to an edge type.
        # This will have functionality if we stitch boundaries rather than defining them after the fact. 
        net2.throat_properties['type']      = sp.repeat(edge,len(net2.throat_properties['type']))
        net1.throat_properties['type']      = sp.concatenate((net1.throat_properties['type'],net2.throat_properties['type']),axis=0)
        

        N1 = net1.pore_properties['divisions']
        N2 = net2.pore_properties['divisions']
        N = sp.vstack((N1,N2))
        
        Nx = sp.sum(N,0)[0]
        Ny = sp.sum(N,0)[1]
        Nz = sp.sum(N,0)[2]
        
        Nx = 5
        Ny = 10
        Nz = 15
        Np = Nx*Ny*Nz
        ind = sp.arange(0,Np)
        
        #Generate throats based on pattern of the adjacency matrix, This is the exact code from the cubic script.
        
        # There are problems when we have the stitched network at a different size than the original network .
        # My goal is to find all connecting neighbors in the unraveled index, but I want to concatenate the stack 
            # to include the upper layer. 
        
        tpore1_1 = ind[(ind%Nx)<(Nx-1)]
        tpore2_1 = tpore1_1 + 1
        tpore1_2 = ind[(ind%(Nx*Ny))<(Nx*(Ny-1))]
        print tpore1_2
        tpore2_2 = tpore1_2 + Nx
        tpore1_3 = ind[(ind%Np)<(Nx*Ny*(Nz-1))]
        tpore2_3 = tpore1_3 + Nx*Ny
        tpore1 = sp.hstack((tpore1_1,tpore1_2,tpore1_3))
        tpore2 = sp.hstack((tpore2_1,tpore2_2,tpore2_3))
        connections = sp.vstack((tpore1,tpore2)).T
        connections = connections[sp.lexsort((connections[:, 1], connections[:, 0]))]
        net1.throat_properties['connections'] = connections
        
        I = net1.throat_properties['connections'][:,0]
        J = net1.throat_properties['connections'][:,1]
        V = sp.ones(len(net1.throat_properties['connections']))
        A = sp.sparse.coo_matrix((V,(I,J))).todense()
        #print A
        
        
        print "Stitch: Method Incomplete"

if __name__ == '__main__':
    test=GenericGeometry(loggername="TestGenerator")        

