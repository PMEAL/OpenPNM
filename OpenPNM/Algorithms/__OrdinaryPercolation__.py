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

import scipy as sp
import numpy as np
import scipy.sparse as sprs
import matplotlib.pyplot as plt

from .__GenericAlgorithm__ import GenericAlgorithm

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

    def _setup(self, invading_fluid, defending_fluid, npts=25, inlets=[0],AL=True,**params):
        self._npts = npts
        self._AL = AL
        self._inv_sites = inlets
        self._fluid_inv = invading_fluid
        self._fluid_def = defending_fluid
        #Create a pore and throat conditions list to store inv_val at which each is invaded
        self._p_inv = sp.zeros((self._net.get_num_pores(),))
        self._p_seq = sp.zeros_like(self._p_inv)
        self._t_inv = sp.zeros((self._net.get_num_throats(),))
        self._t_seq = sp.zeros_like(self._t_inv)        
        #Determine the invasion pressures to apply
        self._t_cap = self._net.get_throat_conditions(self._fluid_inv,'capillary_pressure')
        min_p = sp.amin(self._t_cap)
        max_p = sp.amax(self._t_cap)
        self._inv_points = sp.logspace(sp.log10(min_p),sp.log10(max_p),self._npts)

    def _do_outer_iteration_stage(self):
        #Generate curve from points
        for inv_val in self._inv_points:
            #Apply one applied pressure and determine invaded pores
            self._logger.debug('Applying capillary pressure: '+str(inv_val))
            self._do_one_inner_iteration(inv_val)
        #Store results using networks' get/set method
        self._net.set_pore_conditions('water','inv_Pc',self._p_inv)
        self._net.set_throat_conditions('water','inv_Pc',self._t_inv)
        #Find invasion sequence values (to correspond with IP algorithm)
        self._p_seq = sp.searchsorted(sp.unique(self._p_inv),self._p_inv)
        self._t_seq = sp.searchsorted(sp.unique(self._t_inv),self._t_inv)
        self._net.set_pore_conditions('water','inv_seq',self._p_seq)
        self._net.set_throat_conditions('water','inv_seq',self._t_seq)
        #Remove temporary arrays and adjacency matrices
        del self._net.adjacency_matrix['csr']['invaded']
        
    def _do_one_inner_iteration(self,inv_val):
        r"""
        Determine which throats are invaded at a given applied capillary pressure

        This function uses the scipy.csgraph module for the cluster labeling algorithm (connected_components)

        """
        #Generate a tlist containing boolean values for throat state
        Tinvaded = self._t_cap<=inv_val
        #Fill adjacency matrix with invasion state info
        I = {'invaded': Tinvaded}
        self._net.create_adjacency_matrix(I,sprsfmt='csr',dropzeros=True)
        clusters = sprs.csgraph.connected_components(self._net.adjacency_matrix['csr']['invaded'])[1]
        #Find all pores with at least 1 invaded throat (invaded)
        Pinvaded = sp.zeros_like(clusters,dtype=bool)
        temp = self._net.get_connected_pores(self._net.throat_properties['numbering'])
        temp = temp[Tinvaded]
        temp = sp.hstack((temp[:,0],temp[:,1]))
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
        self._p_inv[(self._p_inv==0)*(pmask)] = inv_val
        #Determine Pc_invaded for throats as well
        temp = self._net.throat_properties['connections']
        tmask = (pmask[temp[:,0]] + pmask[temp[:,1]])*(self._t_inv<=inv_val)
        self._t_inv[(self._t_inv==0)*(tmask)] = inv_val

    def evaluate_trapping(self,outlets=[0]):
        Np = self._net.get_num_pores()
        Nt = self._net.get_num_throats()
        self._p_trap = sp.zeros((Np,),dtype=float)
        inv_points = sp.unique(self._p_inv)
        for inv_val in inv_points[0:-1]:
            #Find clusters of defender pores
            Pinvaded = self._p_inv<=inv_val
            temp = self._net.get_connected_pores(sp.r_[0:Nt])
            PTPstate = sp.sum(Pinvaded[temp],1)
            Tinvaded = (PTPstate>0)*(self._t_inv<=inv_val)
            PTPstate = PTPstate + Tinvaded #0 = all open, 1=1 pore filled, 2=2 pores filled 3=2 pores + 1 throat filled
            I = {'defended': (PTPstate==0)}
            self._net.create_adjacency_matrix(I,sprsfmt='csr',dropzeros=True)
            clusters = sprs.csgraph.connected_components(self._net.adjacency_matrix['csr']['defended'])[1]
            ##Clean up clusters (invaded = -1, defended >=0)
            clusters = clusters*(~Pinvaded) - (Pinvaded)
            #Identify clusters connected to outlet sites
            out_clusters = sp.unique(clusters[outlets])
            trapped_clusters = (~sp.in1d(clusters,out_clusters))*(clusters>=0)
            pmask = trapped_clusters
            self._p_trap[(self._p_trap==0)*(pmask)] = inv_val
        self._p_trap[self._p_trap>0] = 0
        self._net.set_pore_conditions(self._fluid_inv,'trap_Pc',self._p_trap)

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
            print('OP has not been run with this fluid, checking partner fluid')
            try:
                #Apply occupancy to given fluid
                fluid.pore_conditions['occupancy'] = sp.array(~(fluid.partner.pore_conditions['Pc_invaded']<=Pc),ndmin=1)
                fluid.throat_conditions['occupancy'] = sp.array(~(fluid.partner.throat_conditions['Pc_invaded']<=Pc),ndmin=1)
            except:
                raise Exception('It seems that OP has not been run on either fluid')
        #Apply occupancy to partner fluid
        fluid.partner.pore_conditions['occupancy'] = sp.array(~fluid.pore_conditions['occupancy'],ndmin=1)
        fluid.partner.throat_conditions['occupancy'] = sp.array(~fluid.throat_conditions['occupancy'],ndmin=1)
        
    def plot_drainage_curve(self):
          r"""
          Plot drainage capillary pressure curve 
          """
          try:
            PcPoints = sp.unique(self._net.get_pore_conditions(self._fluid_inv,'inv_Pc'))
          except:
            raise ('Cannot print curve: ordinary percolation simulation has not been run')
          Snwp_t = sp.zeros_like(PcPoints)
          Snwp_p = sp.zeros_like(PcPoints)
          Pvol = sum(self._net.pore_properties['volume'])
          Tvol = sum(self._net.throat_properties['volume'])
          for i in range(1,sp.size(PcPoints)):
              Pc = PcPoints[i]
              Snwp_p[i] = sum((self._net.pore_properties['volume'][self._p_inv<=Pc]))/Pvol
              Snwp_t[i] = sum((self._net.throat_properties['volume'][self._t_inv<=Pc]))/Tvol
          plt.plot(PcPoints,Snwp_p,'r.-')
          plt.plot(PcPoints,Snwp_t,'b.-')

if __name__ == '__main__':
    print('no tests yet')
