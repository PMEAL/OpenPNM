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

    def run(self,network,**params):
        r'''
        Parameters
        ----------

        npts: int
            number of simulation steps (pressures); default 25

        inv_sites: array_like
            invasion pores i.e. [1,2,3,4,10]; default = [0]
        '''
        super(OrdinaryPercolation,self).run(network,**params)
        return self

    def _setup(self, invading_fluid,defending_fluid, npts=25, inlets=[0],AL=True,**params):
        self._npts = npts
        self._AL = AL
        self._inv_sites = inlets
        self._fluid_inv = invading_fluid
        self._fluid_def = defending_fluid
        invading_fluid.set_pair(defending_fluid)
        self._fluid_inv.refresh()
        #Apply necessary pore scale physics models
        OpenPNM.Physics.CapillaryPressure.Washburn(self._net,self._fluid_inv)
        #Create a pore and throat conditions list to store inv_val at which each is invaded
        self._fluid_inv.pore_conditions['Pc_invaded'] = np.zeros(self._net.get_num_pores(),np.float)
        self._fluid_inv.throat_conditions['Pc_invaded'] = np.zeros(self._net.get_num_throats(),np.float)
        #Determine the invasion pressures to apply
        min_p = np.min(self._fluid_inv.throat_conditions['Pc_entry'])
        max_p = np.max(self._fluid_inv.throat_conditions['Pc_entry'])
        self._inv_points = np.logspace(np.log10(min_p),np.log10(max_p),self._npts)

    def _do_outer_iteration_stage(self):
        #Generate curve from points
        for inv_val in self._inv_points:
#            self._logger.info("Applying Pc = "+str(int(inv_val)))
            #Apply one applied pressure and determine invaded pores
            self._do_one_inner_iteration(inv_val)
        #Remove temporary arrays and adjacency matrices
        del self._net.adjacency_matrix['csr']['invaded']

    def _do_one_inner_iteration(self,inv_val):
        r"""
        Determines which throats are invaded at a given applied capillary pressure

        This function uses the scipy.csgraph module for the graph theory cluster
        labeling algorithm (connected_components)

        """
        #Generate a tlist containing boolean values for throat state
        Tinvaded = self._fluid_inv.throat_conditions['Pc_entry']<=inv_val
        #Fill adjacency matrix with invasion state info
        I = {'invaded': Tinvaded}
        self._net.create_adjacency_matrix(I,sprsfmt='csr',dropzeros=True)
        clusters = sprs.csgraph.connected_components(self._net.adjacency_matrix['csr']['invaded'])[1]
        #Find all pores with at least 1 invaded throat (invaded)
        Pinvaded = sp.zeros_like(clusters,dtype=bool)
        temp = self._net.get_connected_pores(self._net.throat_properties['numbering'])
        temp = temp[Tinvaded]
        temp = sp.hstack((temp[:,0],temp[:,1]))
#        temp = self._net.get_connected_pores(self._net.throat_properties['numbering'][Tinvaded],flatten=True)
        Pinvaded[temp] = True
        if self._AL:
            #Add injection sites to Pinvaded
            Pinvaded[self._inv_sites] = True
            #Clean up clusters (not invaded = -1, invaded >=0)
            clusters = clusters*(Pinvaded) - (~Pinvaded)
            #Identify clusters connected to invasion sites
            inv_clusters = sp.unique(clusters[self._inv_sites])
        else:
            #Clean up clusters (not invaded = -1, invaded >=0)
            clusters = clusters*(Pinvaded) - (~Pinvaded)
            #All clusters are invasion sites
            inv_clusters = sp.r_[0:self._net.get_num_pores()]
        #Store invasion pressure in pores and throats
        pmask = np.in1d(clusters,inv_clusters)
        #Store result of invasion step
        self._fluid_inv.pore_conditions['Pc_invaded'][(self._fluid_inv.pore_conditions['Pc_invaded']==0)*(pmask)] = inv_val
        #Determine Pc_invaded for throats as well
        temp = self._net.throat_properties['connections']
        tmask = (pmask[temp[:,0]] + pmask[temp[:,1]])*(self._fluid_inv.throat_conditions['Pc_entry']<=inv_val)
        self._fluid_inv.throat_conditions['Pc_invaded'][(self._fluid_inv.throat_conditions['Pc_invaded']==0)*(tmask)] = inv_val

    def evaluate_trapping(self,network,invading_fluid,outlets):
        Np = network.get_num_pores()
        Nt = network.get_num_throats()
        fluid_inv = invading_fluid
        fluid_inv.pore_conditions['Pc_trapped'] = sp.zeros((Np,),dtype=float)
        inv_points = sp.unique(fluid_inv.throat_conditions['Pc_invaded'])
        for inv_val in inv_points[0:-1]:
            #Find clusters of defender pores
            Pinvaded = fluid_inv.pore_conditions['Pc_invaded']<=inv_val
            temp = network.get_connected_pores(sp.r_[0:Nt])
            PTPstate = sp.sum(Pinvaded[temp],1)
            Tinvaded = (PTPstate>0)*(fluid_inv.throat_conditions['Pc_entry']<=inv_val)
            PTPstate = PTPstate + Tinvaded #0 = all open, 1=1 pore filled, 2=2 pores filled 3=2 pores + 1 throat filled
            I = {'defended': (PTPstate==0)}
            network.create_adjacency_matrix(I,sprsfmt='csr',dropzeros=True)
            clusters = sprs.csgraph.connected_components(network.adjacency_matrix['csr']['defended'])[1]
            ##Clean up clusters (invaded = -1, defended >=0)
            clusters = clusters*(~Pinvaded) - (Pinvaded)
            #Identify clusters connected to outlet sites
            out_clusters = sp.unique(clusters[outlets])
            trapped_clusters = (~sp.in1d(clusters,out_clusters))*(clusters>=0)
            pmask = trapped_clusters
            fluid_inv.pore_conditions['Pc_trapped'][(fluid_inv.pore_conditions['Pc_trapped']==0)*(pmask)] = inv_val
        fluid_inv.pore_conditions['Pc_invaded'][fluid_inv.pore_conditions['Pc_trapped']>0]=0

    def update_occupancy(fluid,Pc=0):
        r"""
        Updates the fluid occupancy status of invading and defending fluids as determined by the OP algorithm

        Parameters
        ----------
        fluid : OpenPNM Fluid Object
            This can be either the invading or defending fluid used.
        """
        #Apply occupancy to given fluid
        try:
            fluid.pore_conditions['occupancy'] = sp.array(fluid.pore_conditions['Pc_invaded']<=Pc,ndmin=1)
            fluid.throat_conditions['occupancy'] = sp.array(fluid.throat_conditions['Pc_invaded']<=Pc,ndmin=1)
        except:
            print ('OP has not been run with this fluid, checking partner fluid')
            try:
                #Apply occupancy to given fluid
                fluid.pore_conditions['occupancy'] = sp.array(~(fluid.partner.pore_conditions['Pc_invaded']<=Pc),ndmin=1)
                fluid.throat_conditions['occupancy'] = sp.array(~(fluid.partner.throat_conditions['Pc_invaded']<=Pc),ndmin=1)
            except:
                raise Exception('It seems that OP has not been run on either fluid')
        #Apply occupancy to partner fluid
        fluid.partner.pore_conditions['occupancy'] = sp.array(~fluid.pore_conditions['occupancy'],ndmin=1)
        fluid.partner.throat_conditions['occupancy'] = sp.array(~fluid.throat_conditions['occupancy'],ndmin=1)

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
