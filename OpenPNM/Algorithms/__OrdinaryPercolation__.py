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

    loglevel : integer, optional
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)
    loggername : string, optional
        Set the name of the logger to be output on the console. Defaults to class name.

    Note
    ----
    Some info about OP?
    """

    def __init__(self, **kwargs):
        r"""

        """
        super(OrdinaryPercolation,self).__init__(**kwargs)
        self._logger.debug("Create Drainage Percolation Algorithm Object")

    def run(self, invading_fluid, defending_fluid, npts=25, inlets=[0],AL=True,**params):
        self._npts = npts
        self._AL = AL
        self._inv_sites = inlets
        self._fluid_inv = invading_fluid
        self._fluid_def = defending_fluid
        #Create a pore and throat conditions list to store inv_val at which each is invaded
        self._p_inv = sp.zeros((self._net.num_pores(),))
        self._p_seq = sp.zeros_like(self._p_inv)
        self._t_inv = sp.zeros((self._net.num_throats(),))
        self._t_seq = sp.zeros_like(self._t_inv)
        #Determine the invasion pressures to apply
        self._t_cap = self._net.get_throat_data(phase=self._fluid_inv,prop='capillary_pressure')
        min_p = sp.amin(self._t_cap)
        max_p = sp.amax(self._t_cap)
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
        self._net.create_adjacency_matrix(data=Tinvaded,prop='invaded',sprsfmt='csr',dropzeros=True)
        clusters = sprs.csgraph.connected_components(self._net.adjacency_matrix['csr']['invaded'])[1]
        #Find all pores with at least 1 invaded throat (invaded)
        Pinvaded = sp.zeros_like(clusters,dtype=bool)
        nums = self._net.get_throat_data(prop='numbering')
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
        temp = self._net.get_throat_data(prop='connections')
        tmask = (pmask[temp[:,0]] + pmask[temp[:,1]])*(self._t_inv<=inv_val)
        self._t_inv[(self._t_inv==0)*(tmask)] = inv_val

    def evaluate_trapping(self,outlets=[0]):
        r"""
        Finds trapped pores and throats after a full ordinary percolation/drainage has been run

        Parameters
        ---------
        outlets : array_like
            A list of pores that define the wetting phase outlets.  Disconnection from these outlets results in trapping.

            TODO: Ideally this should update the inv_Pc property and set trapped pores to a value of inf.
            This would allow update_occupancy and plotting to work.
        """
        Np = self._net.num_pores()
        Nt = self._net.num_throats()
        self._p_trap = sp.zeros((Np,),dtype=float)
        inv_points = sp.unique(self._p_inv)
        conns = self._net.find_connected_pores(sp.r_[0:Nt])
        for inv_val in inv_points[0:-1]:
            #Find clusters of defender pores
            Pinvaded = self._p_inv<=inv_val
            PTPstate = sp.sum(Pinvaded[conns],1)
            Tinvaded = (PTPstate>0)*(self._t_inv<=inv_val)
            PTPstate = PTPstate + Tinvaded #0 = all open, 1=1 pore filled, 2=2 pores filled 3=2 pores + 1 throat filled
            self._net.create_adjacency_matrix(data=(PTPstate==0),prop='defended',sprsfmt='csr',dropzeros=True)
            clusters = sprs.csgraph.connected_components(self._net.adjacency_matrix['csr']['defended'])[1]
            ##Clean up clusters (invaded = -1, defended >=0)
            clusters = clusters*(~Pinvaded) - (Pinvaded)
            #Identify clusters connected to outlet sites
            out_clusters = sp.unique(clusters[outlets])
            trapped_clusters = (~sp.in1d(clusters,out_clusters))*(clusters>=0)
            pmask = trapped_clusters
            self._p_trap[(self._p_trap==0)*(pmask)] = inv_val
        self._p_trap[self._p_trap>0] = 0
        self.set_pore_data(phase=self._fluid_inv,prop='inv_Pc',data=self._p_trap)

    def update(self,Pc=0):
        r"""
        Updates the occupancy status of invading and defending fluids as determined by the OP algorithm

        """
        #Apply occupancy to invading fluid
        p_inv = self.get_pore_data(prop='inv_Pc')
        self._net.set_pore_data(phase=self._fluid_inv,prop='inv_Pc',data=p_inv)
        t_inv = self.get_throat_data(prop='inv_Pc')
        self.set_throat_data(phase=self._fluid_inv,prop='inv_Pc',data=t_inv)
        #Find invasion sequence values (to correspond with IP algorithm)
        p_seq = self.get_pore_data(prop='inv_seq')
        self._net.set_pore_data(phase=self._fluid_inv,prop='inv_seq',data=p_seq)
        t_seq = self.get_throat_data(prop='inv_seq')
        self._net.set_throat_data(phase=self._fluid_inv,prop='inv_seq',data=t_seq)
        #Remove temporary arrays and adjacency matrices
        p_inv = self.get_pore_data(prop='inv_Pc')<=Pc
        t_inv = self.get_throat_data(prop='inv_Pc')<=Pc
        self._net.set_pore_data(phase=self._fluid_inv,prop='occupancy',data=p_inv)
        self._net.set_throat_data(phase=self._fluid_inv,prop='occupancy',data=t_inv)
        #Apply occupancy to defending fluid
        self._net.set_pore_data(phase=self._fluid_def,prop='occupancy',data=~p_inv)
        self._net.set_throat_data(phase=self._fluid_def,prop='occupancy',data=~t_inv)

    def plot_drainage_curve(self):
          r"""
          Plot drainage capillary pressure curve
          """
          try:
            PcPoints = sp.unique(self._net.get_pore_data(phase=self._fluid_inv,prop='inv_Pc'))
          except:
            raise Exception('Cannot print drainage curve: ordinary percolation simulation has not been run')
          Snwp_t = sp.zeros_like(PcPoints)
          Snwp_p = sp.zeros_like(PcPoints)
          Pvol = self._net.get_pore_data(prop='volume')
          Tvol = self._net.get_throat_data(prop='volume')
          Pvol_tot = sum(Pvol)
          Tvol_tot = sum(Tvol)
          for i in range(0,sp.size(PcPoints)):
              Pc = PcPoints[i]
              Snwp_p[i] = sum(Pvol[self._p_inv<=Pc])/Pvol_tot
              Snwp_t[i] = sum(Tvol[self._t_inv<=Pc])/Tvol_tot
          plt.plot(PcPoints,Snwp_p,'r.-')
          plt.plot(PcPoints,Snwp_t,'b.-')
          plt.xlim(xmin=0)

if __name__ == '__main__':
    print('no tests yet')
