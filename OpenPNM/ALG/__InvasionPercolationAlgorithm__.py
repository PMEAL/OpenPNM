#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: CEF PNM Team
# License: TBD
# Copyright (c) 2012

#from __future__ import print_function

"""
module __InvasionPercolationAlgorithm__: Invasion Percolation Algorithm
========================================================================

.. warning:: The classes of this module should be loaded through the 'ALG.__init__.py' file.

"""

import OpenPNM
import scipy as sp
import numpy as np
import scipy.sparse as sprs
from time import clock
import heapq
import itertools

from __GenericAlgorithm__ import GenericAlgorithm

class InvasionPercolationAlgorithm(GenericAlgorithm):
    r"""   
    
    Invasion_Percolation - Class to run IP algorithm on constructed networks
    
    Parameters
    ----------
    inlets : list of integers (default: [0])
        list of inlet nodes
    outlets : list of integers (default: [1])
        list of outlet nodes
    end_condition : string('breaktrhough')
        choice between 'breakthrough' and 'total'
            
    Examples
    --------
    
    
    TODO:
    """
    
    def __init__(self,net=OpenPNM.NET.GenericNetwork,inlets=[0],outlets=[1],end_condition='breakthrough',**kwords):
        r"""
        
        """
        super(InvasionPercolationAlgorithm,self).__init__(net = net,**kwords)
        self._logger.info("Create IP Algorithm Object")
        self._logger.info("\t end condition: "+end_condition)
        self.inlets = inlets
        self.outlets = outlets
        self.end_condition = end_condition
        self.counter = 0

    def _do_outer_iteration_stage(self):
        r"""
        Executes the outer iteration stage
        """
        self._logger.info("Outer Iteration Stage ")
        self.pseq = 1
        self.tseq = 1
        self._setup_for_IP()
        self._condition_update()
        #self.Tinv = np.zeros(self._net.get_num_throats())
        while self.condition:
            self._do_one_outer_iteration()
        self._net.set_pore_property(name="IP_Pinv",ndarray=self.Pinv,columns=None)
        self._net.set_throat_property(name="IP_Tinv",ndarray=self.Tinv,columns=None)
        self._net.set_pore_property(name="IP_Pseq",ndarray=self.Psequence,columns=None)
        self._net.set_throat_property(name="IP_Tseq",ndarray=self.Tsequence,columns=None)

    def _do_one_outer_iteration(self):
        r"""
        One iteration of an outer iteration loop for an algorithm 
        (e.g. time or parametric study)
        """
        if (sp.mod(self.counter,500)==False):
            self._logger.info("Outer Iteration (counter = "+str(self.counter)+")")
        self._update_network_state()
        self._do_inner_iteration_stage()
        self._condition_update()
        self.counter += 1
#        if np.mod(self.counter,100)==0:
#            print 'on a multiple of 100'
#            print self.counter
#            print len(np.nonzero(self.Tinv)[0])
#            print len(self.Tinv)
     
    def _do_inner_iteration_stage(self):
        r"""
        Executes the inner iteration stage
        """
        self._logger.debug("  Inner Iteration Stage: ")        
        
        self.plast = len(np.nonzero(self.Pinv)[0])
        for i in np.unique(self.Pinv[self.Pinv>0]):
            self.current_cluster = i
            self._do_one_inner_iteration()
        self.pnew = len(np.nonzero(self.Pinv)[0])
        self.tseq += 1
        if self.pnew>self.plast:
            self.pseq += 1
      

    def _do_one_inner_iteration(self):
        r"""
        Executes one inner iteration
        """
        self._logger.debug("    Inner Iteration")
        #if self.plists[self.current_cluster-1]!=0:                    
        neighbors = self._net.get_neighbor_throats(self.plists[self.current_cluster-1])
        for j in neighbors:
            if (j not in self.tlists[self.current_cluster-1]):
                heapq.heappush(self.tpoints[self.current_cluster-1],(self._net.throat_properties['Pc_entry'][j],j))
                self.tlists[self.current_cluster-1].append(j)
        tinvade = heapq.heappop(self.tpoints[self.current_cluster-1])[1]
        self.Tsequence[tinvade] = self.tseq
        Pores = self._net.get_connected_pores(tinvade)
        if np.in1d(Pores,np.nonzero(self.Pinv)[0]).all():
            clusterNumber = min(self.Pinv[Pores])
            self.Tinv[tinvade] = clusterNumber
            if self.Pinv[Pores[0]]!=self.Pinv[Pores[1]] :
                maxCluster = np.max(self.Pinv[Pores])
                self.Pinv[np.where(self.Pinv==maxCluster)[0]] = clusterNumber
                self.Tinv[np.where(self.Tinv==maxCluster)[0]] = clusterNumber
                if self.current_cluster==clusterNumber:
                    self.plists[clusterNumber-1] = self.plists[maxCluster-1]
                self.plists[maxCluster-1] = 0                               
                self.tlists[clusterNumber-1] = self.tlists[clusterNumber-1] + self.tlists[maxCluster-1]
                self.tlists[maxCluster-1] = []
                self.tpoints[clusterNumber-1] = heapq.merge(self.tpoints[clusterNumber-1],self.tpoints[maxCluster-1])
                self.tpoints[clusterNumber-1] = list(k for k,v in itertools.groupby(self.tpoints[clusterNumber-1]))
                self.tpoints[maxCluster-1] = []
        else:
            self.Tinv[tinvade] = self.current_cluster
            self.NewPore = Pores[self.Pinv[Pores][:,0]==0]
            self.plists[self.current_cluster-1] = self.NewPore
            self.Pinv[self.NewPore] = self.current_cluster
            self.Psequence[self.NewPore] = self.pseq 

        
    def _setup_for_IP(self):
        r"""
        Determines cluster labelling and condition for completion
        """
        # if empty, add Pc_entry to throat_properties
        tdia = self._net.throat_properties['diameter']
        # calculate Pc_entry from diameters
        surface_tension = 0.486
        contact_angle = 180
        Pc_entry = -4*surface_tension*np.cos(np.deg2rad(contact_angle))/(tdia/1000/1000)  
        self._net.set_throat_property(name="Pc_entry",ndarray=Pc_entry,columns=None)
        # Creating an array for invaded Pores(Np long, 0 for uninvaded, cluster number for inaveded)
        self.Pinv = np.zeros((self._net.get_num_pores(),1),dtype=np.int32)
        # Creating an array for invaded throats(Nt long, 0 for uninvaded, cluster number for inaveded)
        self.Tinv = np.zeros((self._net.get_num_throats(),1),dtype=np.int32)
        # Creating an array for tracking invaded Pores(Np long, 0 for uninvaded,sequence for inaveded)
        self.Psequence = np.zeros((self._net.get_num_pores(),1),dtype=np.int32)
        # Creating an array for tracking invaded throats(Nt long, 0 for uninvaded,sequence for inaveded)
        self.Tsequence = np.zeros((self._net.get_num_throats(),1),dtype=np.int32)
        # Creating an array for tracking the last invaded pore in each cluster.
        # its length is equal to the maximum number of possible clusters.
        self.plists = np.zeros((len(self.inlets),1),dtype=np.int32)
        # Iterator variables for sequences and cluster numbers
        clusterNumber = 1
        # Creating an empty list to store the list of potential throats for invasion in each cluster.
        # its length is equal to the maximum number of possible clusters. 
        self.tlists = []
        # Creating a list for each cluster to store both potential throat and corresponding throat value
        self.tpoints = []
        # Initializing invasion percolation for each possible cluster
        for i in self.inlets:
            self.Pinv[i] = clusterNumber
            interface_throat_numbers = self._net.get_neighbor_throats(np.where(self.Pinv==clusterNumber)[0])
            interface_throat_pressures = self._net.throat_properties["Pc_entry"][interface_throat_numbers]            
            Interface= zip(interface_throat_pressures,interface_throat_numbers)
            heapq.heapify(Interface)
            self.tlists.append(interface_throat_numbers.tolist())
            tinvade = heapq.heappop(Interface)[1]
            self.tpoints.append(Interface)
            self.Tinv[tinvade] = clusterNumber            
            self.Tsequence[tinvade] = self.tseq
            self.Pores = self._net.get_connected_pores(tinvade)
            self.NewPore = self.Pores[self.Pinv[self.Pores][:,0]==0]
            self.plists[clusterNumber-1] = self.NewPore
            self.Pinv[self.NewPore] = clusterNumber
            self.Psequence[self.NewPore] = self.pseq
            clusterNumber += 1
        self.tseq += 1
        self.pseq += 1
            
    def _condition_update(self):
        if self.end_condition == 'breakthrough':
            self.condition = self.NewPore not in self.outlets
        elif self.end_condition == 'total':
            self.condition = not self.Tinv.all()
        
    def _update_network_state(self):
        r"""
        This does nothing
        """

class InvasionPercolationAlgorithmTiming(InvasionPercolationAlgorithm):
    r"""   
    
    Invasion_Percolation with cluster growth timing - Class to run IP algorithm on constructed networks
    
    Parameters
    ----------
    inlets : list of integers (default: [0])
        list of inlet nodes
    outlets : list of integers (default: [1])
        list of outlet nodes
    end_condition : string('breaktrhough')
        choice between 'breakthrough' and 'total'
            
    Examples
    --------
    
    
    TODO:
    1)  Currently requires all clusters to start out with identical flow rates, currently a value of 1 unit of volume per unit of time
    2)  Currently requires cap volume function to be a linear function of pressure. Things will get a bit more complicated if we shouldn't assume this.

    Suggested Improvements:
    a) Allow input of cluster flow-rates (condensation rates)
    b) Allow updating of cluster flow-rates (this will require a delta-t calculation at each step, instead of a total t calculation).
    c) Allow for a non-linear relationship between pressure and throat-cap volume.
        
    """
    
    def __init__(self,net=OpenPNM.NET.GenericNetwork,inlets=[0],outlets=[1],end_condition='breakthrough',**kwords):
        r"""
        
        """
        super(InvasionPercolationAlgorithm,self).__init__(net = net,**kwords)
        self._logger.info("Create IP Algorithm Object")
        self._logger.info("\t end condition: "+end_condition)
        self.inlets = inlets
        if sp.size(inlets) == 1:
            self.inlets = [inlets]
        self.outlets = outlets
        self.end_condition = end_condition
        self.counter = 0
        self.condition = 1

    def _setup_for_IP(self):
        r"""
        Determines cluster labelling and condition for completion
        """
        self._logger.debug( '+='*25)
        self._logger.debug( 'INITIAL SETUP (STEP 1)')
        # if empty, add Pc_entry to throat_properties
        tdia = self._net.throat_properties['diameter']
        # calculate Pc_entry from diameters
        surface_tension = 0.486
        contact_angle = 180
        Pc_entry = -4*surface_tension*np.cos(np.deg2rad(contact_angle))/(tdia)#/1000/1000)  
        self._net.set_throat_property(name="Pc_entry",ndarray=Pc_entry,columns=None)
        # calculate Volume_coef for each throat
        self.Tvol_coef = tdia*tdia*tdia*np.pi/6/Pc_entry
        # Creating an array for invaded Pores(Np long, 0 for uninvaded, cluster number for inaveded)
        self.Pinv = np.zeros((self._net.get_num_pores(),1),dtype=np.int32)
        self.Pinv_original = np.zeros((self._net.get_num_pores(),1),dtype=np.int32)
        # Creating an array for invaded throats(Nt long, 0 for uninvaded, cluster number for inaveded)
        self.Tinv = np.zeros((self._net.get_num_throats(),1),dtype=np.int32)
        # Creating arrays for tracking invaded Pores(Np long, 0 for uninvaded, sequence for inaveded)
        self.Psequence = np.zeros((self._net.get_num_pores(),1),dtype=np.int32)
        # Creating arrays for tracking invaded Pores(Np long, -1 for uninvaded, simulation time for inaveded)
        self.Ptime = np.zeros((self._net.get_num_pores(),1),dtype=np.float64)-1
        # Creating arrays for tracking invaded throats(Nt long, 0 for uninvaded, sequence for inaveded)
        self.Tsequence = np.zeros((self._net.get_num_throats(),1),dtype=np.int32)
        # Creating arrays for tracking invaded Pores(Np long, -1 for uninvaded, simulation time for inaveded)
        self.Ttime = np.zeros((self._net.get_num_throats(),1),dtype=np.float64)-1
        # Creating an array for tracking the last invaded pore in each cluster.
        # its length is equal to the maximum number of possible clusters.
        #self.plists = np.zeros((len(self.inlets),1),dtype=np.int32)
        # Iterator variables for sequences and cluster numbers
        clusterNumber = 1
        # Determine how many clusters there are
        clusterCount = 0
        for i in self.inlets:
            clusterCount += 1
        # Storage for cluster information
        self.cluster_data = {}
        self.cluster_data['flow_rate'] = np.ones((clusterCount),dtype=np.float64)
        self.cluster_data['haines_pressure'] = np.zeros((clusterCount),dtype=np.float64)
        self.cluster_data['haines_time'] = np.zeros((clusterCount),dtype=np.float64)
        self.cluster_data['vol_coef'] = np.zeros((clusterCount),dtype=np.float64)
        self.cluster_data['cap_volume'] = np.zeros((clusterCount),dtype=np.float64)
        self.cluster_data['pore_volume'] = np.zeros((clusterCount),dtype=np.float64)
        self.cluster_data['haines_throat'] = np.zeros((clusterCount),dtype=np.int32)
        self.cluster_data['active'] = np.ones((clusterCount),dtype=np.int8)
        self.cluster_data['transform'] = np.zeros((clusterCount),dtype=np.int16)
        for i in range(clusterCount):
            self.cluster_data['transform'][i] = i+1
        # Creating an empty list to store the list of potential throats for invasion in each cluster.
        # its length is equal to the maximum number of possible clusters. 
        self.tlists = []
        # Creating a list for each cluster to store both potential throat and corresponding throat value
        self.tpoints = []
        # Initializing invasion percolation for each possible cluster
        for i in self.inlets:
            # Calculate total volume in all invaded pores
            self.cluster_data['pore_volume'][clusterNumber-1] = np.sum(self._net.pore_properties['volume'][i])
            # Label all invaded pores with their cluster
            self.Pinv[i] = clusterNumber
            self.Pinv_original[i] = clusterNumber
            # Label all inlet pores as invaded
            self.Psequence[i] = self.tseq
            self.Ptime[i] = self.sim_time
            # Find all throats that border invaded pores
            interface_throat_numbers = self._net.get_neighbor_throats(np.where(self.Pinv==clusterNumber)[0])
            # Sum all interfacial throats' volume coeffients for throat cap volume calculation
            self.cluster_data['vol_coef'][clusterNumber-1] = np.sum(self.Tvol_coef[interface_throat_numbers])
            # Make a list of all entry pressures of the interfacial throats            
            interface_throat_pressures = self._net.throat_properties["Pc_entry"][interface_throat_numbers]#[0]         
            # Zip pressures and numbers together so that HeapQ can work its magic
            self._logger.debug('interface throat(s) found:')
            self._logger.debug(interface_throat_numbers)                   
            self._logger.debug( 'interface throat pressure(s):')
            self._logger.debug(interface_throat_pressures)
            Interface= zip(interface_throat_pressures,interface_throat_numbers)
            # Turn the zipped throat interfaces object into a heap            
            heapq.heapify(Interface)
            # Add to the total list of interface throats in the system
            self.tlists.append(interface_throat_numbers.tolist())
            # Add to the total list of invaded interface throats in the system         
            self.tpoints.append(Interface)   
            # Pop off the first entry (lowest pressure) on the throat info list
            invaded_throat_info = Interface[0]
            # Determine pressure at Haines Jump
            self.cluster_data['haines_pressure'][clusterNumber-1] = invaded_throat_info[0]
            # Calculate cap_volume at Haines Jump
            self.cluster_data['cap_volume'][clusterNumber-1] = self.cluster_data['haines_pressure'][clusterNumber-1]*self.cluster_data['vol_coef'][clusterNumber-1]
            # Calculate time at Haines Jump
            self.cluster_data['haines_time'][clusterNumber-1] = (self.cluster_data['pore_volume'][clusterNumber-1]+
                                        self.cluster_data['cap_volume'][clusterNumber-1])/self.cluster_data['flow_rate'][clusterNumber-1]
            # Record invaded throat
            self.cluster_data['haines_throat'][clusterNumber-1] = invaded_throat_info[1]
            ##self.Tinv[tinvade] = clusterNumber            
            ##self.Tsequence[tinvade] = self.tseq
            ##self.Pores = self._net.get_connected_pores(tinvade)
            ##self.NewPore = self.Pores[self.Pinv[self.Pores][:,0]==0]
            ##self.plists[clusterNumber-1] = self.NewPore
            ##self.Pinv[self.NewPore] = clusterNumber
            ##self.Psequence[self.NewPore] = self.pseq
            clusterNumber += 1
        self._logger.debug(self.cluster_data['pore_volume'])
        self._logger.debug( 'cap volumes')
        self._logger.debug( self.cluster_data['cap_volume'])
        self._logger.debug( 'haines_throats')
        self._logger.debug( self.cluster_data['haines_throat'])
        self._logger.debug( 'max throat cap volumes')
        self._logger.debug( self.Tvol_coef*self._net.throat_properties["Pc_entry"])
        self._logger.debug( '+='*25)
        
        
        self.tseq += 1
        self.pseq += 1

    def _do_outer_iteration_stage(self):
        r"""
        Executes the outer iteration stage
        """
        self._logger.info("Outer Iteration Stage ")
        self.pseq = 1
        self.tseq = 1
        self.NewPore = -1
        # Time keeper
        self.sim_time = 0
        self._setup_for_IP()
        self._condition_update()
        #self.Tinv = np.zeros(self._net.get_num_throats())
        while self.condition:
            self._do_one_outer_iteration()
        self._net.set_pore_property(name="IP_Pinv",ndarray=self.Pinv,columns=None)
        self._net.set_pore_property(name="IP_Pinv_original",ndarray=self.Pinv_original,columns=None)
        self._net.set_throat_property(name="IP_Tinv",ndarray=self.Tinv,columns=None)
        self._net.set_pore_property(name="IP_Pseq",ndarray=self.Psequence,columns=None)
        self._net.set_throat_property(name="IP_Tseq",ndarray=self.Tsequence,columns=None)
        self._net.set_pore_property(name="IP_Ptime",ndarray=self.Ptime,columns=None)
        self._net.set_throat_property(name="IP_Ttime",ndarray=self.Ttime,columns=None)

    def _do_one_outer_iteration(self):
        r"""
        One iteration of an outer iteration loop for an algorithm 
        (e.g. time or parametric study)
        """
        if (sp.mod(self.counter,500)==False):
            self._logger.info("Outer Iteration (counter = "+str(self.counter)+")")
        self._update_network_state()
        self._do_inner_iteration_stage()
        self._condition_update()
        self.counter += 1
#        if np.mod(self.counter,100)==0:
#            self._logger.debug( 'on a multiple of 100'
#            self._logger.debug( self.counter
#            self._logger.debug( len(np.nonzero(self.Tinv)[0])
#            self._logger.debug( len(self.Tinv)
     
    def _do_inner_iteration_stage(self):
        r"""
        Executes the inner iteration stage
        """
        self._logger.debug("  Inner Iteration Stage: ")        
        
        self.plast = len(np.nonzero(self.Pinv)[0])
        # determine the cluster with the earliest Haines time
        self.current_cluster = 1 + self.cluster_data['haines_time'].tolist().index(min(self.cluster_data['haines_time']))
        # update simulation clock
        self._logger.debug( 'sim time = ')
        self._logger.debug(self.sim_time)
        self._logger.debug(' haines time:')
        self._logger.debug( self.cluster_data['haines_time'])
        # The code really messes up when the [0] isn't in the next line. sim_time seems to just point to a place on the haines time array
        self.sim_time = min(self.cluster_data['haines_time'])
        self._logger.debug( 'sim time after update= ')
        self._logger.debug(self.sim_time)
        # run through the Haines Jump steps        
        self._do_one_inner_iteration()
        self.pnew = len(np.nonzero(self.Pinv)[0])
        self.tseq += 1
        if self.pnew>self.plast:
            self.pseq += 1
      

    def _do_one_inner_iteration(self):
        r"""
        Executes one inner iteration
        """
        self._logger.debug("    Inner Iteration")             
        # Fill throat and connecting pore
        # Pop out the largest throat (lowest Pcap) in the list, read the throat number
        tinvade = heapq.heappop(self.tpoints[self.current_cluster-1])[1]
        self._logger.debug( ' ')
        self._logger.debug( '--------------------------------------------------')
        self._logger.debug( 'STEP')
        self._logger.debug(self.tseq)
        self._logger.debug( 'trying to access cluster: ')
        self._logger.debug(self.current_cluster)
        self._logger.debug( 'when these clusters are active active: ')
        self._logger.debug(sp.nonzero(self.cluster_data['active'])[0])
        self._logger.debug( 'Haines at throat,time: ')
        self._logger.debug(tinvade)
        self._logger.debug(self.sim_time)

        # Mark throat as invaded
        self.Tsequence[tinvade] = self.tseq
        self.Ttime[tinvade] = self.sim_time
        # Remove throat's contribution to the vol_coef
        self.cluster_data['vol_coef'][self.current_cluster-1] = self.cluster_data['vol_coef'][self.current_cluster-1]-self.Tvol_coef[tinvade]
        # Mark pore as invaded
        Pores = self._net.get_connected_pores(tinvade)
        # If both pores are already invaded:
        if np.in1d(Pores,np.nonzero(self.Pinv)[0]).all():
            self.NewPore = -1
            # Label invaded throat with smaller cluster number
            #find cluster 1
            clusters = self.cluster_data['transform'][self.Pinv[Pores]-1]
            self._logger.debug('clusters = ')
            self._logger.debug(clusters)
            self.current_cluster = min(clusters)[0]
            self.Tinv[tinvade] = self.current_cluster
            # if pores are from 2 different clusters:
            if self.Pinv[Pores[0]]!=self.Pinv[Pores[1]] :
                # find name of larger cluster number                
                maxCluster = max(clusters)[0]
                self._logger.info(' ')
                self._logger.info('CLUSTERS COMBINING:')
                self._logger.info(self.current_cluster)
                self._logger.info(maxCluster)
                self._logger.info('at time')
                self._logger.info(self.sim_time)
                # update the cluster transform
                self.cluster_data['transform'][self.cluster_data['transform']==maxCluster] = [self.current_cluster][0]  
                # relabel all pores and throats from larger number with smaller number
                self.Pinv[np.where(self.Pinv==maxCluster)[0]] = self.current_cluster
                self.Tinv[np.where(self.Tinv==maxCluster)[0]] = self.current_cluster
                # append the list of throats for the other cluster to the current cluster                              
                self.tlists[self.current_cluster-1] = self.tlists[self.current_cluster-1] + self.tlists[maxCluster-1]
                # delete the throat lists on the other cluster     
                self.tlists[maxCluster-1] = []
                # merge the heaps of throat information
                self.tpoints[self.current_cluster-1] = list(heapq.merge(self.tpoints[self.current_cluster-1],self.tpoints[maxCluster-1]))
                # update the clusters' vol_coefs
                self.cluster_data['vol_coef'][self.current_cluster-1] += self.cluster_data['vol_coef'][maxCluster-1]
                self.cluster_data['vol_coef'][maxCluster-1] = 0  
                # update the clusters' pore volume
                self.cluster_data['pore_volume'][self.current_cluster-1] += self.cluster_data['pore_volume'][maxCluster-1]
                self.cluster_data['pore_volume'][maxCluster-1] = 0
                # update the clusters' flowrates
                self.cluster_data['flow_rate'][self.current_cluster-1] += self.cluster_data['flow_rate'][maxCluster-1]
                self.cluster_data['flow_rate'][maxCluster-1] = 0
                self._logger.debug( 'new flowrate for cluster ')
                self._logger.debug(self.current_cluster)
                self._logger.debug('is')
                self._logger.debug(self.cluster_data['flow_rate'][self.current_cluster-1])
                # check if either was inactive (broke through already)
                if self.cluster_data['active'][maxCluster-1] + self.cluster_data['active'][self.current_cluster-1]<2:
                    self._logger.debug('making clusters ')
                    self._logger.debug(self.current_cluster)
                    self._logger.debug('and')
                    self._logger.debug(maxCluster)
                    self._logger.debug('inactive due to one being inactive already')
                    self._logger.debug(self.cluster_data['active'][self.current_cluster-1])
                    self._logger.debug(self.cluster_data['active'][maxCluster-1])
                    self.cluster_data['active'][maxCluster-1] = 0
                    self.cluster_data['active'][self.current_cluster-1] = 0
                    self.cluster_data['haines_time'][self.current_cluster-1] = 100000000000000000000000000000000
                    self._logger.info(' ')
                    self._logger.info('CLUSTER MERGED WITH A BREAKTHROUGH CLUSTER')
                self._logger.info('making cluster ')
                self._logger.info(maxCluster)
                self._logger.info('inactive due to merge')
                # update the old cluster's activity and time
                self.cluster_data['haines_time'][maxCluster-1] = 100000000000000000000000000000000
                self.cluster_data['active'][maxCluster-1] = 0 
                # NO IDEA WHAT THIS LINE DOES PLEASE HELP MAHMOUD
                #self.tpoints[self.current_cluster-1] = list(k for k,v in itertools.groupby(self.tpoints[self.current_cluster-1]))
                self.tpoints[maxCluster-1] = []

        else:
            # label invaded throat with current cluster
            self.Tinv[tinvade] = self.current_cluster
            # find univaded pore, NewPore
            self.NewPore = Pores[self.Pinv[Pores][:,0]==0][0]
            self._logger.debug( ' ')            
            self._logger.debug( 'INVADING PORE: ')
            self._logger.debug(self.NewPore)
            self._logger.debug('the other pore is one of: ')
            self._logger.debug(Pores)
            self._logger.debug( 'position: ')
            self._logger.debug(self._net.pore_properties['coords'][self.NewPore])
            # label that pore as invaded
            self.Pinv[self.NewPore] = self.current_cluster
            self.Pinv_original[self.NewPore] = self.current_cluster
            self.Ptime[self.NewPore] = self.sim_time
            self.Psequence[self.NewPore] = self.tseq 
            # update self.cluster_data.['pore_volume']
            self.cluster_data['pore_volume'][self.current_cluster-1] += self._net.pore_properties['volume'][self.NewPore]
            # Make a list of all throats neighboring pores in the cluster
            # Update interface list        
            neighbors = self._net.get_neighbor_throats(self.NewPore)
            for j in neighbors:
                # If a throat is not labelled as invaded by the cluster, it must be an interfacial throat
                if (j not in self.tlists[self.current_cluster-1]):
                    self._logger.debug( 'new throat:')
                    self._logger.debug(j)
                    self._logger.debug('connecting pores:')
                    self._logger.debug(self._net.get_connected_pores(j))
                    # Add this throat data (pressure, number) to this cluster's "heap" of throat data.
                    heapq.heappush(self.tpoints[self.current_cluster-1],(self._net.throat_properties['Pc_entry'][j],j))
                    # Add new throat number to throat list for this cluster
                    self.tlists[self.current_cluster-1].append(j)
                    # Update the cluster's vol_coef
                    self.cluster_data['vol_coef'][self.current_cluster-1] = self.cluster_data['vol_coef'][self.current_cluster-1]+self.Tvol_coef[j]
        # Find next Haines Jump info
        # Make sure you are not re-invading a throat
        while self.Tinv[self.tpoints[self.current_cluster-1][0][1]] > 0:
            if self.tpoints[self.current_cluster-1] == []:
                self._logger.debug( 'making cluster ')
                self._logger.debug(self.current_cluster)
                self._logger.debug('inactive due to tpoints = [] ')
                self.cluster_data['active'][self.current_cluster-1] = 0
                self.cluster_data['haines_time'][self.current_cluster-1] = 100000000000000000000000000000000
                break
            tremove = heapq.heappop(self.tpoints[self.current_cluster-1])[1]
            self.cluster_data['vol_coef'][self.current_cluster-1] = self.cluster_data['vol_coef'][self.current_cluster-1]-self.Tvol_coef[tremove]
        next_throat = self.tpoints[self.current_cluster-1][0][1]
        self.cluster_data['haines_throat'][self.current_cluster-1] = next_throat
        self.cluster_data['haines_pressure'][self.current_cluster-1] = self.tpoints[self.current_cluster-1][0][0]
        self.cluster_data['cap_volume'][self.current_cluster-1] = self.cluster_data['haines_pressure'][self.current_cluster-1]*self.cluster_data['vol_coef'][self.current_cluster-1]
            
        # Calculate the new Haines jump time
        self._logger.debug( 'haines time before last stage:')
        self._logger.debug( self.cluster_data['haines_time'])
        if self.tpoints[self.current_cluster-1] == []:
            self._logger.debug('making cluster ')
            self._logger.debug(self.current_cluster)
            self._logger.debug('inactive due to self.tpoints being empty for that cluster')
            self.cluster_data['active'][self.current_cluster-1] = 0
            self.cluster_data['haines_time'][self.current_cluster-1] = 100000000000000000000000000000000
        if self.cluster_data['active'][self.current_cluster-1] == 1:
            self.cluster_data['haines_time'][self.current_cluster-1] = (self.cluster_data['pore_volume'][self.current_cluster-1]+self.cluster_data['cap_volume'][self.current_cluster-1])/self.cluster_data['flow_rate'][self.current_cluster-1]
        if self.cluster_data['haines_time'][self.current_cluster-1] < self.sim_time:
            self.cluster_data['haines_time'][self.current_cluster-1] = self.sim_time + 0.01
        self._logger.debug('haines time at the end of the throat stuff')
        self._logger.debug(self.cluster_data['haines_time'])
            
    def _condition_update(self):
        if self.end_condition == 'breakthrough':
            if self.NewPore in self.outlets:
                self._logger.info( ' ')
                self._logger.info( 'BREAKTHROUGH AT PORE: ')
                self._logger.info(self.NewPore)
                self._logger.info('in cluster ')
                self._logger.info(self.current_cluster)
                self._logger.info('at time')
                self._logger.info(self.sim_time)
                self.cluster_data['active'][self.current_cluster-1] = 0
                self.cluster_data['haines_time'][self.current_cluster-1] = 100000000000000000000000000000000
            if np.sum(self.cluster_data['active']) == 0:
                self._logger.info( ' ')
                self._logger.info( 'SIMULATION FINISHED; no more active clusters')
                self._logger.info('at time')
                self._logger.info(self.sim_time)
                self.condition = 0
        elif self.end_condition == 'total':
            self.condition = not self.Tinv.all()            
            
            
if __name__ =="__main__":
    
    from time import clock
    start=clock()
    print "="*50
    print "= Example: Create random network and run an invasion\n= percolation algorithm"
    print "-"*50
    print "- * generate a simple cubic network"    
    #sp.random.seed(1)
    pn = OpenPNM.GEN.Cubic(domain_size=[200,1,20],lattice_spacing=1,btype = [0,1,0]).generate()
    #pn = OpenPNM.GEN.Delaunay(domain_size=[30,30,10],num_pores = 5000 ,btype = [1,1,0]).generate()
    print "+"*50
    print "Sample generated at t =",clock()-start,"seconds."
    print "+"*50
    #print "- * Assign pore volumes"
    #pore_volumes=sp.random.rand(pn.get_num_pores())
    #pore_volumes[range(pn.get_num_pores([0]),pn.get_num_pores())]=0
    #pn.set_pore_property(name='volume',ndarray=pore_volumes,columns = None)  
    print '- * Assign boundary pore volumes = 0'
    pn.pore_properties['diameter'][pn.pore_properties['type']>0] = 0
        
    
    print "- * Define inlet and outlet faces"
    inlets = sp.nonzero(pn.pore_properties['type']==1)[0]
    outlets = sp.nonzero(pn.pore_properties['type']==6)[0]
    inlets2 = sp.unique((inlets[sp.random.randint(sp.size(inlets,0))],inlets[sp.random.randint(sp.size(inlets,0))],
                       inlets[sp.random.randint(sp.size(inlets,0))],inlets[sp.random.randint(sp.size(inlets,0))],
                       inlets[sp.random.randint(sp.size(inlets,0))],inlets[sp.random.randint(sp.size(inlets,0))],
                       inlets[sp.random.randint(sp.size(inlets,0))],inlets[sp.random.randint(sp.size(inlets,0))],
                       inlets[sp.random.randint(sp.size(inlets,0))],inlets[sp.random.randint(sp.size(inlets,0))]))
    print inlets2
    #print "- * assign random pore and throat diameters"
    #pn.pore_properties['diameter'] = sp.random.rand(pn.get_num_pores(),1)
    #pn.throat_properties['diameter'] = sp.random.rand(pn.get_num_throats(),1)
    
    print "- * Run Invasion percolation algorithm"
    IP = InvasionPercolationAlgorithmTiming(net=pn,inlets=inlets2,outlets=outlets,loglevel=20,loggername="TestInvPercAlg")
    IP.run()
    print "+"*50
    print "IP completed at t =",clock()-start,"seconds."
    print "+"*50
    print "- * Save output to IP.vtp"
    OpenPNM.IO.NetToVtp(net = pn,filename="IP.vtp")
    
    print "="*50
    print "Program Finished at t = ",clock()-start,"seconds."
    print "="*50