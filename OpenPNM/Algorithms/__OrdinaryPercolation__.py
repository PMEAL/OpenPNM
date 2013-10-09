#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jeff Gostick (jeff@gostick.ca)
# License: TBD
# Copyright (c) 2013


"""
module __OrdinaryPercolation__: Ordinary Percolation Algorithm
========================================================================

.. warning:: The classes of this module should be loaded through the 'Algorithms.__init__.py' file.

"""

import OpenPNM
import scipy as sp
import numpy as np
import scipy.sparse as sprs
import matplotlib.pyplot as plt

from __GenericAlgorithm__ import GenericAlgorithm

class OrdinaryPercolation(GenericAlgorithm):
    r"""
    Calculates a capillary pressure curve by looping through a list
    of capillary pressures
    
    This function produces a pore_properties array called 'Pc_invaded' which contains the
    pressure at which a given pore was invaded. This list can be useful for
    reproducing the simulation for plotting or determining late pore filling.

    Parameters
    ----------
    
    npts: int
        number of simulation steps (pressures); default 25

    inv_sites: array_like
        invasion pores i.e. [1,2,3,4,10]; default = [0]
    

    """
    
    def __init__(self, net ,npts=25,inv_sites=[0],**kwargs):
        r"""
        
        """
        super(OrdinaryPercolation,self).__init__(net = net,**kwargs)
        self._logger.debug("Create Drainage Percolation Algorithm Object")
        self._npts = npts
        self._inv_src = inv_sites
        self._setup()
        
    def _setup(self):
        #Create a pore and throat list to store inv_val at which each is invaded
        self._net.pore_properties['Pc_invaded'] = np.zeros(self._net.get_num_pores(),np.float)
        self._net.throat_properties['Pc_invaded'] = np.zeros(self._net.get_num_throats(),np.float)
        #Create a throat list to temporarily store the invasion state of throats
        self._net.throat_properties['invaded'] = np.zeros(self._net.get_num_throats())
        #Determine the invasion pressures to apply
        min_p = np.min(self._net.throat_properties['Pc_entry'])
        max_p = np.max(self._net.throat_properties['Pc_entry'])
        self._inv_points = np.logspace(np.log10(min_p),np.log10(max_p),self._npts)
        print self._inv_points
        
    def _do_outer_iteration_stage(self):
        #Generate curve from points
        for inv_val in self._inv_points:
            self._logger.info("Applying Pc = "+str(int(inv_val)))
            #Apply one applied pressure and determine invaded pores
            inv_clusters = self._do_one_inner_iteration(inv_val)
            #Store result of invasion step
            self._net.pore_properties['Pc_invaded'][(self._net.pore_properties['Pc_invaded']==0)*(inv_clusters>0)] = inv_val
            r"""
            TODO:
            Tracking the pressure at which each throat is invaded has not been
            implimented yet.  This means that calculation of invaded volume is
            based only on volume of the pores.
            """
#            invaded_pores = self._net.pore_properties['numbering'][inv_clusters>0]
#            connected_throats = self._net.get_neighbor_throats(invaded_pores)
#            self._net.throat_properties['Pc_invaded'][self._net.throat_properties['Pc_invaded']==0] = inv_val
        del self._net.throat_properties['invaded']
            
    def _do_one_inner_iteration(self,inv_val):
        r"""
        Determines which throats are invaded at a given applied capillary pressure
        
        This function uses the scipy.csgraph module for the graph theory cluster
        labeling algorithm (connected_components)
        
        """
        #Generate a tlist containing boolean values for throat state
        self._net.throat_properties['invaded'] = self._net.throat_properties['Pc_entry']<inv_val
        #Fill adjacency matrix with invasion state info
        self._net.create_adjacency_matrix('invaded',sprsfmt='csr',dropzeros=True)
        #This step seems to be very slow!
        clusters = sprs.csgraph.connected_components(self._net.adjacency_matrix['csr']['invaded'])[1]
        #Clean up (not invaded = -2, invaded >0)
        clusters = (clusters[0:]>=0)*(clusters[0:]+1)
        #Identify clusters connected to invasion sites
        inj_clusters = np.in1d(self._inv_src,self._net.pore_properties['numbering'][clusters>0])
        #Trim non-connected clusters
        inv_clusters2 = sp.zeros([np.size(clusters,0)],np.int32)
        inv_clusters2[np.in1d(clusters,inj_clusters)] = 1
        temp = sp.unique(inj_clusters[sp.nonzero(inj_clusters)])
        inv_clusters = sp.zeros([np.size(clusters,0)],np.int32)
        for i in range(0,np.shape(temp)[0]):
            pores=sp.where(clusters==temp[i])[0]
            inv_clusters[pores] = temp[i]
        return(inv_clusters)
            
    def _plot_results(self):
        r"""
        Plot the numerical output of the OP algorithm
        """
        PcPoints = np.unique(self._net.pore_properties['Pc_invaded'])
        Snwp = np.zeros_like(PcPoints)
        Ps = np.where(self._net.pore_properties['type']==0)
        for i in range(1,np.size(PcPoints)):
            Pc = PcPoints[i]
            Snwp[i] = sum((self._net.pore_properties['Pc_invaded'][Ps]<Pc)*(self._net.pore_properties['volume'][Ps]))/sum(self._net.pore_properties['volume'][Ps])
        plt.plot(PcPoints,Snwp,'r.-')
        plt.show(block = False)
        self._results = np.vstack((PcPoints,Snwp)).T
        
        
if __name__ == '__main__':
    print "Create a test"