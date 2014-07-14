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
    Simulates a capillary drainage experiment by looping through a list of 
    capillary pressures.  

    Parameters
    ----------
    network : OpenPNM Network Object
        The network upon which the simulation will be run
    name : string, optional
        The name to assign to the Algorithm Object

    Notes
    -----
    To run this algorithm, use 'setup()' to provide the necessary simulation 
    parameters, and then use 'run()' to execute it.  Use 'update()' to send
    the results of the simulation out of the algorithm object.
    """

    def __init__(self, **kwargs):
        r"""

        """
        super(OrdinaryPercolation,self).__init__(**kwargs)
        self._logger.debug("Create Drainage Percolation Algorithm Object")
        
    def setup(self,invading_fluid = None,defending_fluid = None,inlets = [0],npts = 25,capillary_pressure = 'capillary_pressure',AL=True,**params):
        r'''
        '''
        # Parse params
        self._fluid_inv = invading_fluid
        self._fluid_def = defending_fluid
        try: self._inv_sites
        except: self._inv_sites = inlets
        self._npts = npts
        self._AL = AL
        self._p_cap = capillary_pressure

    def run(self):
        r'''
        '''
        #See if setup has been run
        try: capillary_pressure = self._p_cap
        except: 
            raise Exception('setup has not been run, cannot proceed')
        #Create pore and throat conditions lists to store inv_val at which each is invaded
        self._p_inv = sp.zeros((self._net.num_pores(),))
        self._p_seq = sp.zeros_like(self._p_inv)
        self._t_inv = sp.zeros((self._net.num_throats(),))
        self._t_seq = sp.zeros_like(self._t_inv)
        #Determine the invasion pressures to apply
        self._t_cap = self._fluid_inv['throat.capillary_pressure']
        min_p = sp.amin(self._t_cap)*0.98  # nudge min_p down slightly
        max_p = sp.amax(self._t_cap)*1.02  # bump max_p up slightly
        self._inv_points = sp.logspace(sp.log10(min_p),sp.log10(max_p),self._npts)
        self._do_outer_iteration_stage()

    def _do_outer_iteration_stage(self):
        #Generate curve from points
        for inv_val in self._inv_points:
            #Apply one applied pressure and determine invaded pores
            self._logger.info('Applying capillary pressure: '+str(inv_val))
            self._do_one_inner_iteration(inv_val)
        #Store results using networks' get/set method
        self.set_pore_data(prop='inv_Pc',data=self._p_inv)
        self.set_throat_data(prop='inv_Pc',data=self._t_inv)
        #Find invasion sequence values (to correspond with IP algorithm)
        self._p_seq = sp.searchsorted(sp.unique(self._p_inv),self._p_inv)
        self._t_seq = sp.searchsorted(sp.unique(self._t_inv),self._t_inv)
        self.set_pore_data(prop='inv_seq',data=self._p_seq)
        self.set_throat_data(prop='inv_seq',data=self._t_seq)

    def _do_one_inner_iteration(self,inv_val):
        r"""
        Determine which throats are invaded at a given applied capillary pressure

        This function uses the scipy.csgraph module for the cluster labeling algorithm (connected_components)

        """
        #Generate a tlist containing boolean values for throat state
        Tinvaded = self._t_cap<=inv_val
        #Fill adjacency matrix with invasion state info
        clusters = self._net.find_clusters(Tinvaded)
        #Find all pores with at least 1 invaded throat (invaded)
        Pinvaded = sp.zeros_like(clusters,dtype=bool)
        nums = self._net.throats()
        temp = self._net.find_connected_pores(nums)
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
            inv_clusters = sp.r_[0:self._net.num_pores()]
        #Store invasion pressure in pores and throats
        pmask = np.in1d(clusters,inv_clusters)
        #Store result of invasion step
        self._p_inv[(self._p_inv==0)*(pmask)] = inv_val
        #Determine Pc_invaded for throats as well
        temp = self._net.get_throat_data(prop='conns')
        tmask = (pmask[temp[:,0]] + pmask[temp[:,1]])*(Tinvaded)
        self._t_inv[(self._t_inv==0)*(tmask)] = inv_val

    def evaluate_trapping(self, outlets):
        r"""
        Finds trapped pores and throats after a full ordinary
        percolation drainage has been run

        Parameters
        ---------
        outlets : array_like
            A list of pores that define the wetting phase outlets.
            Disconnection from these outlets results in trapping.

        """
        self._p_trap = sp.zeros_like(self._p_inv, dtype=float)
        self._t_trap = sp.zeros_like(self._t_inv, dtype=float)
        try:
            inv_points = sp.unique(self._p_inv)  # Get points used in OP
        except:
            self._logger.error('Orindary percolation has not been run!')
            raise Exception('Aborting algorithm')
        tind = self._net.get_throat_indices()
        conns = self._net.find_connected_pores(tind)
        for inv_val in inv_points[0:-1]:
            #Find clusters of defender pores
            Pinvaded = self._p_inv <= inv_val
            Cstate = sp.sum(Pinvaded[conns], axis=1)
            Tinvaded = self._t_inv <= inv_val
            Cstate = Cstate + Tinvaded #0 = all open, 1=1 pore filled, 2=2 pores filled 3=2 pores + 1 throat filled
#            self._net.create_adjacency_matrix(data=(Cstate == 0), prop='defended', sprsfmt='csr', dropzeros=True)
#            clusters = sprs.csgraph.connected_components(self._net.adjacency_matrix['csr']['defended'])[1]
            clusters = self._net.find_clusters(Cstate==0)
            ##Clean up clusters (invaded = -1, defended >=0)
            clusters = clusters*(~Pinvaded) - (Pinvaded)
            #Identify clusters connected to outlet sites
            out_clusters = sp.unique(clusters[outlets])
            trapped_pores = ~sp.in1d(clusters, out_clusters)
            self._p_trap[(self._p_trap == 0)[trapped_pores]] = inv_val
            trapped_throats = self._net.find_neighbor_throats(trapped_pores)
            self._t_trap[(self._t_trap == 0)[trapped_throats]] = inv_val       
            trapped_throats = sp.where(Cstate==2)[0]
            self._t_trap[(self._t_trap == 0)[trapped_throats]] = inv_val
        self._p_inv[self._p_trap > 0] = sp.inf
        self._t_inv[self._t_trap > 0] = sp.inf        
        self.set_pore_data(prop='inv_Pc', data=self._p_inv)
        self.set_throat_data(prop='inv_Pc', data=self._t_inv)

    def update(self, Pc=0, occupancy='occupancy'):
        r"""
        Updates the occupancy status of invading and defending fluids
        as determined by the OP algorithm

        """
        #Apply invasion pressure to invading fluid
        p_inv = self.get_pore_data(prop='inv_Pc')
        self._fluid_inv.set_pore_data(prop='inv_Pc',data=p_inv)
        t_inv = self.get_throat_data(prop='inv_Pc')    
        self._fluid_inv.set_throat_data(prop='inv_Pc',data=t_inv)
        #Find invasion sequence values (to correspond with IP algorithm)
        p_seq = self.get_pore_data(prop='inv_seq')
        self._fluid_inv.set_pore_data(prop='inv_seq',data=p_seq)
        t_seq = self.get_throat_data(prop='inv_seq')
        self._fluid_inv.set_throat_data(prop='inv_seq',data=t_seq)
        
        
        p_inv = self.get_pore_data(prop='inv_Pc')<=Pc
        t_inv = self.get_throat_data(prop='inv_Pc')<=Pc
        #Apply occupancy to invading fluid
        temp = sp.array(p_inv,dtype=sp.float64,ndmin=1)
        self._fluid_inv.set_pore_data(prop=occupancy,data=temp)
        temp = sp.array(t_inv,dtype=sp.float64,ndmin=1)
        self._fluid_inv.set_throat_data(prop=occupancy,data=temp)
        #Apply occupancy to defending fluid
        if self._fluid_def != None:
            temp = sp.array(~p_inv,dtype=sp.float64,ndmin=1)
            self._fluid_def.set_pore_data(prop=occupancy,data=temp)
            temp = sp.array(~t_inv,dtype=sp.float64,ndmin=1)
            self._fluid_def.set_throat_data(prop=occupancy,data=temp)

    def plot_drainage_curve(self,
                            pore_volume='volume',
                            throat_volume='volume'):
          r"""
          Plot drainage capillary pressure curve
          """
          try:
            PcPoints = sp.unique(self.get_pore_data(prop='inv_Pc'))
          except:
            raise Exception('Cannot print drainage curve: ordinary percolation simulation has not been run')
          Snwp_t = sp.zeros_like(PcPoints)
          Snwp_p = sp.zeros_like(PcPoints)
          Pvol = self._net.get_pore_data(prop=pore_volume)
          Tvol = self._net.get_throat_data(prop=throat_volume)
          Pvol_tot = sum(Pvol)
          Tvol_tot = sum(Tvol)
          for i in range(0,sp.size(PcPoints)):
              Pc = PcPoints[i]
              Snwp_p[i] = sum(Pvol[self._p_inv<=Pc])/Pvol_tot
              Snwp_t[i] = sum(Tvol[self._t_inv<=Pc])/Tvol_tot
          plt.plot(PcPoints,Snwp_p,'r.-')
          plt.plot(PcPoints,Snwp_t,'b.-')
          plt.xlim(xmin=0)
          plt.show()

if __name__ == '__main__':
    print('no tests yet')
