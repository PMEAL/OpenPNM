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
from time import clock

from __GenericAlgorithm__ import GenericAlgorithm

class OrdinaryPercolation(GenericAlgorithm):
    r"""
    Simulates a capillary drainage experiment by looping through a list of capillary pressures

    This function produces a pore_properties array called 'Pc_invaded' which contains the pressure at which a given pore was invaded. This list can be useful for reproducing the simulation for plotting or determining late pore filling.

    Parameters
    ----------

    npts: int
        number of simulation steps (pressures); default 25

    inv_sites: array_like
        invasion pores i.e. [1,2,3,4,10]; default = [0]

    Note
    ----
    Some info about OP?
    """

    def __init__(self, **kwargs):
        r"""

        """
        super(OrdinaryPercolation,self).__init__(**kwargs)
        self._logger.debug("Create Drainage Percolation Algorithm Object")

    def _setup(self, npts=25, inv_sites=[0]):
        self._npts = npts
        self._inv_sites = inv_sites
        #Create a pore and throat list to store inv_val at which each is invaded
        self._net.pore_conditions['Pc_invaded'] = np.zeros(self._net.get_num_pores(),np.float)
        self._net.throat_conditions['Pc_invaded'] = np.zeros(self._net.get_num_throats(),np.float)
        #Create a throat list to temporarily store the invasion state of throats
        self._net.throat_properties['invaded'] = np.zeros(self._net.get_num_throats())
        #Determine the invasion pressures to apply
        min_p = np.min(self._net.throat_conditions['Pc_entry'])
        max_p = np.max(self._net.throat_conditions['Pc_entry'])
        self._inv_points = np.logspace(np.log10(min_p),np.log10(max_p),self._npts)

    def _do_outer_iteration_stage(self):
        #Generate curve from points
        for inv_val in self._inv_points:
            self._logger.info("Applying Pc = "+str(int(inv_val)))
            #Apply one applied pressure and determine invaded pores
            pmask = self._do_one_inner_iteration(inv_val)
            #Store result of invasion step
            self._net.pore_conditions['Pc_invaded'][(self._net.pore_conditions['Pc_invaded']==0)*(pmask)] = inv_val
            #Determine Pc_invaded for throats as well
            temp = self._net.throat_properties['connections']
            tmask = (pmask[temp[:,0]] + pmask[temp[:,1]])*(self._net.throat_conditions['Pc_entry']<inv_val)
            self._net.throat_conditions['Pc_invaded'][(self._net.throat_conditions['Pc_invaded']==0)*(tmask)] = inv_val
        #Remove temporary arrays and adjacency matrices
        del self._net.throat_conditions['invaded']
        del self._net.adjacency_matrix['csr']['invaded']

    def _do_one_inner_iteration(self,inv_val):
        r"""
        Determines which throats are invaded at a given applied capillary pressure

        This function uses the scipy.csgraph module for the graph theory cluster
        labeling algorithm (connected_components)

        """
        #Generate a tlist containing boolean values for throat state
        self._net.throat_conditions['invaded'] = self._net.throat_conditions['Pc_entry']<inv_val
        sum(self._net.throat_conditions['invaded'])
        #Fill adjacency matrix with invasion state info
        I = {'invaded': self._net.throat_conditions['Pc_entry']<inv_val}
        self._net.create_adjacency_matrix(I,sprsfmt='csr',dropzeros=True)
        clusters = sprs.csgraph.connected_components(self._net.adjacency_matrix['csr']['invaded'])[1]
        #Find all pores with at least 1 invaded throat
        Pinvaded = sp.zeros_like(clusters,dtype=bool)
        temp = self._net.get_connected_pores(np.r_[0:self._net.get_num_throats()])
        temp = temp[self._net.throat_conditions['invaded']]
        temp = sp.hstack((temp[:,0],temp[:,1]))
        Pinvaded[temp] = True
        #Add injection sites to Pinvaded for ALOP
        if sp.shape(self._inv_sites)<sp.shape(clusters):
            Pinvaded[self._inv_sites] = True
        #Clean up clusters (not invaded = -1, invaded >=0)
        clusters = clusters*(Pinvaded) - (~Pinvaded)
        #Identify clusters connected to invasion sites
        inv_clusters = clusters[self._inv_sites]
        return np.in1d(clusters,inv_clusters)*(clusters>=0)


if __name__ == '__main__':
    print ''
    print ''    
    print '************Testing OrdinaryPercolation Algorithm**************'
    clock()
    print "="*50
    print "= Example: Create random network and run an ordinary\n= percolation algorithm"
    print "-"*50   
    params = {
        'domain_size'               : [1,1,1],  #physical network size [meters]
        'divisions'                 : [10,10,10], #Number of pores in each direction
        'lattice_spacing'           : [],  #spacing between pores [meters]
        'stats_pores'   : {  'name' : 'weibull_min', #Each statistical package takes different params, so send as dict
                            'shape' : 1.5,
                              'loc' : 6e-6,
                            'scale' : 2e-5},
        'stats_throats' : {  'name' : 'weibull_min',
                            'shape' : 1.5,
                              'loc' : 6e-6,
                            'scale' : 2e-5},
        'btype'                     : [0,0,0],  #boundary type to apply to opposing faces [x,y,z] (1=periodic)
        }
    
    print "- * Generate a simple cubic network" 
    pn = OpenPNM.Geometry.Cubic().generate(**params)
    print "- * Assign capillary pressures to throats"
    pn.throat_properties['Pc_entry'] = OpenPNM.Physics.CapillaryPressure.Washburn(pn,0.072,110)
    inlets = [0]
    print "- * Run Ordinary percolation algorithm"
    exp = OpenPNM.Algorithms.OrdinaryPercolation()
    exp.run(pn, npts=50, inv_sites=inlets)
    print "+"*50
    print "- * Completed OP algorithm in a",pn.get_num_pores(),'pore network \n-   with',exp._npts,'points in',np.round(clock(),decimals=2),'seconds.'
    print "+"*50
    print 