#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: CEF PNM Team
# License: TBD
# Copyright (c) 2012

#from __future__ import print_function

"""
module __InvasionPercolation__: Invasion Percolation Algorithm
========================================================================

.. warning:: The classes of this module should be loaded through the 'Algorithms.__init__.py' file.

"""
import OpenPNM
import scipy as sp
import numpy as np
import heapq
from OpenPNM.Utilities import misc

from OpenPNM.Algorithms.__GenericAlgorithm__ import GenericAlgorithm


class InvasionPercolation(GenericAlgorithm):
    r"""
    Invasion percolation with cluster growth timing - Class to run IP algorithm on constructed networks

    Parameters
    ----------
    network : Descendent of OpenPNM.Network.GenericNetwork
        A valid network for this algorithm
    name : string
        The name this algorithm will go by
  
    Examples
    --------
    >>> pn = OpenPNM.Network.TestNet()
    >>> geo = OpenPNM.Geometry.TestGeometry(network=pn,pores=pn.pores(),throats=pn.throats())
    >>> phase1 = OpenPNM.Phases.TestPhase(network=pn)
    >>> phase2 = OpenPNM.Phases.TestPhase(network=pn)
    >>> phys1 = OpenPNM.Physics.TestPhysics(network=pn, phase=phase1,pores=pn.pores(),throats=pn.throats())
    >>> phys2 = OpenPNM.Physics.TestPhysics(network=pn, phase=phase2,pores=pn.pores(),throats=pn.throats())
    >>> IP = OpenPNM.Algorithms.InvasionPercolation(network=pn, name='IP')
    >>> IP.run(invading_phase=phase1, defending_phase=phase2, inlets=pn.pores('top'), outlets=pn.pores('bottom'),report=0)
         IP algorithm at 0 % completion at 0.0 seconds
         IP algorithm at 100% completion at  0.0  seconds
    >>> IP.update_results()
    >>> max(phase1['pore.IP_inv_seq']) #unless something changed with our test objects, this should print "60"
    60

    Suggested Improvements ::

        a) Allow updating of cluster flow-rates (this will require a delta-t calculation at each step, instead of a total t calculation).
        b) Allow for a non-linear relationship between pressure and throat-cap volume.


    """
    def __init__(self,**kwords):
        r'''
        '''
        super(InvasionPercolation,self).__init__(**kwords)
        self._logger.info("Create IP Algorithm Object")

    def run(self,invading_phase,
               defending_phase,
               inlets=[0],
                outlets=[-1],
                end_condition='breakthrough',
                capillary_pressure='capillary_pressure',
                pore_volume_name='volume',
                throat_volume_name='volume',
                throat_diameter_name='diameter',
                timing='ON',
                inlet_flow=1e-12, #default flowrate is 1 nanoliter/sec/cluster
                report=20):
        r"""
        Runs the IP algorithm

        Parameters
        ----------
        invading_phase : OpenPNM Phase Object
            phase which will displace defending phase
        defending_phase : OpenPNM Phase Object
            phase which will be displaced by invading phase
        inlets : list of integers (default: [0])
            list of inlet nodes
        outlets : list of integers (default: [-1])
            list of outlet nodes
        end_condition : string('breakthrough')
            choice between 'breakthrough' and 'total'
        capillary_pressure : string('capillary_pressure')
            name given to throat capillary pressure property
        pore_volume_name : string('volume')
            name given to pore volume property
        throat_diameter_name : string('diameter')
            name given to throat diameter property
        timing : string ('ON')
            turns volume and flowrate calculations 'ON' or 'OFF'
        inlet_flow : float (1)
            m3/s for each cluster (affects timestamp of pore filling)
        report : int (20)
            percentage multiple at which a progress report is printed

    
        Returns
        -------
        The algorithm will aquire the following pore data ::
    
            invaded          : True for invaded, False for uninvaded
            defended         : True for uninvaded, False for invaded
            cluster_final    : 0 for uninvaded, merged cluster number for invaded
            cluster_original : 0 for uninvaded, original cluster number for invaded
            inv_seq          : 0 for uninvaded, simulation step for invaded
            inv_time         : 0 for uninvaded, simulation time for invaded
            inv_sat          : 0 for uninvaded, simulation saturation for invaded
    
        and throat data ::
        
            invaded          : True for invaded, False for uninvaded
            defended         : True for uninvaded, False for invaded
            cluster_final    : 0 for uninvaded, merged cluster number for invaded
            inv_seq          : 0 for uninvaded, simulation step for invaded
            inv_time         : 0 for uninvaded, simulation time for invaded
            inv_sat          : 0 for uninvaded, simulation saturation for invaded
            inv_Pc           : throat capillary pressures 
      
        """

        self._logger.info("\t end condition: "+end_condition)
        self._inlets = inlets
        self._outlets = outlets
        if end_condition=='total':
            self._brkevent = []
        self._inlet_flow = inlet_flow
        try:    self._phase = self._net._phases[invading_phase]
        except: self._phase = invading_phase
        try:    self._phase_def = self._net._phases[defending_phase]
        except: self._phase_def = defending_phase

        if sp.size(inlets) == 1:
            self._inlets = [inlets]
        if sp.size(outlets) == 1:
            self._outlets = [outlets]
        self._end_condition = end_condition
        self._counter = 0
        self._condition = 1
        self._rough_increment = report
        if report == 0:
            self._rough_increment = 100
        self._timing = timing=='ON'
        self._capillary_pressure_name = capillary_pressure
        self._pore_volume_name = pore_volume_name
        self._throat_volume_name = throat_volume_name
        self._throat_diameter_name = throat_diameter_name

        super(InvasionPercolation,self).run()

    def _setup_for_IP(self):
        r"""
        Determines cluster labelling and condition for completion
        """
        self._clock_start = misc.tic()
        self._logger.debug( '+='*25)
        self._logger.debug( 'INITIAL SETUP (STEP 1)')
        # if empty, add Pc_entry to throat_properties
        tdia = self._net['throat.'+self._throat_diameter_name]
        # calculate Pc_entry from diameters
        try:
            self['throat.inv_Pc'] = self._phase['throat.'+self._capillary_pressure_name]
        except:
            self._logger.error('Capillary pressure not assigned to invading phase '+self._phase.name
                +', check for capillary pressure in defending phase '+self._phase_def.name +' instead')
            try:
                self['throat.inv_Pc'] = self._phase_def['throat.'+self._capillary_pressure_name]
                self._phase['throat.'+self._capillary_pressure_name] = self._phase_def['throat.'+self._capillary_pressure_name]
            except:
                self._logger.error('Capillary pressure neither assigned to defending phase '+self._phase_def.name
                    +' nor to invading phase '+self._phase.name)
        if self._timing:
            # calculate Volume_coef for each throat
            self._Tvol_coef = tdia*tdia*tdia*np.pi/12/self['throat.inv_Pc']
        # Creating an array for invaded Pores(Np long, 0 for uninvaded, cluster number for inaveded)
        self['pore.cluster_final'] = 0
        self['pore.cluster_original'] = 0
        # Creating an array for invaded throats(Nt long, 0 for uninvaded, cluster number for inaveded)
        self['throat.cluster_final'] = 0
        # Creating arrays for tracking invaded Pores(Np long, 0 for uninvaded, sequence for inaveded)
        self['pore.inv_seq'] =0
        if self._timing:
            # Creating arrays for tracking invaded Pores(Np long, -1 for uninvaded, simulation time for inaveded)
            self['pore.inv_time'] = -1.
        # Creating arrays for tracking invaded throats(Nt long, 0 for uninvaded, sequence for inaveded)
        self['throat.inv_seq'] = 0
        if self._timing:
            # Creating arrays for tracking invaded Pores(Np long, -1 for uninvaded, simulation time for inaveded)
            self['throat.inv_time'] = -1.
        # Iterator variables for sequences and cluster numbers
        clusterNumber = 1
        # Determine how many clusters there are
        self._clusterCount = 0
        for i in self._inlets:
            self._clusterCount += 1
        # Storage for cluster information
        self._cluster_data = {}
        if self._timing:
            self._cluster_data['flow_rate'] = np.ones((self._clusterCount),dtype=float)*self._inlet_flow
            self._cluster_data['haines_pressure'] = np.zeros((self._clusterCount),dtype=float)
            self._cluster_data['haines_time'] = np.zeros((self._clusterCount),dtype=float)
            self._cluster_data['vol_coef'] = np.zeros((self._clusterCount),dtype=float)
            self._cluster_data['cap_volume'] = np.zeros((self._clusterCount),dtype=float)
            self._cluster_data['pore_volume'] = np.zeros((self._clusterCount),dtype=float)
            self._cluster_data['throat_volume'] = np.zeros((self._clusterCount),dtype=float)
        self._cluster_data['haines_throat'] = np.zeros((self._clusterCount),dtype=int)
        self._cluster_data['active'] = np.ones((self._clusterCount),dtype=int)
        self._cluster_data['transform'] = np.zeros((self._clusterCount),dtype=int)
        for i in range(self._clusterCount):
            self._cluster_data['transform'][i] = i+1
        # Creating an empty list to store the list of potential throats for invasion in each cluster.
        # its length is equal to the maximum number of possible clusters.
        self._tlists = [[] for i in self._inlets]
        # Creating a list for each cluster to store both potential throat and corresponding throat value
        self._tpoints = [[] for i in self._inlets]
        # Initializing invasion percolation for each possible cluster
        self._pore_volumes = self._net['pore.'+self._pore_volume_name]
        self._throat_volumes = self._net['throat.'+self._throat_volume_name]
        for pores in self._inlets:
            if sp.shape(pores) == ():
                pores = [pores]
            # Label all invaded pores with their cluster
            self['pore.cluster_original'][pores] = clusterNumber
            # Label all inlet pores as invaded
            self['pore.inv_seq'][pores] = self._tseq
            if self._timing:
                self['pore.inv_time'][pores] = self._sim_time
            # Find all throats that border invaded pores
            interface_throat_numbers = self._net.find_neighbor_throats(pores)
            self.cluster_update(clusterNumber,pores,[],interface_throat_numbers)
            clusterNumber += 1
        if self._timing:
            self._logger.debug( 'pore volumes')
            self._logger.debug(self._cluster_data['pore_volume'])
            self._logger.debug( 'cap volumes')
            self._logger.debug( self._cluster_data['cap_volume'])
        self._logger.debug( 'haines_throats')
        self._logger.debug( self._cluster_data['haines_throat'])
        self._tseq += 1
        self._pseq += 1
        self._current_cluster = 0
        # Calculate the distance between the inlet and outlet pores
        self._outlet_position = np.average(self._net.get_data(prop='coords',pores='all')[self._outlets],0)
        if any([sp.shape(i) > () for i in self._inlets]): 
            inlets = []
            for i in self._inlets:
                inlets = sp.union1d(inlets,i)
            inlets = sp.array(inlets,int)
        else:
            inlets = self._inlets
        inlet_position = np.average(self._net.get_data(prop='coords',pores='all')[inlets],0)
        dist_sqrd = (self._outlet_position-inlet_position)*(self._outlet_position-inlet_position)
        self._initial_distance = np.sqrt(dist_sqrd[0]+dist_sqrd[1]+dist_sqrd[2])
        self._logger.debug( 'initial distance')
        self._logger.debug( self._initial_distance)
        self._current_distance = self._initial_distance
        self._percent_complete = np.round((self._initial_distance-self._current_distance)/self._initial_distance*100, decimals = 1)
        self._logger.info( 'percent complete')
        self._logger.info( self._percent_complete)
        self._rough_complete = 0
        print('     IP algorithm at',np.int(self._rough_complete),'% completion at',np.round(misc.toc(quiet=True)),'seconds')
        self._logger.debug( '+='*25)

    def _do_outer_iteration_stage(self):
        r"""
        Executes the outer iteration stage
        """
        self._logger.info("Outer Iteration Stage ")
        self._pseq = 1
        self._tseq = 1
        self._NewPore = -1
        # Time keeper
        self._sim_time = 0
        self._setup_for_IP()
        self._condition_update()
        #self['throat.cluster_final'] = np.zeros(self._net.num_throats())
        while self._condition:
            self._do_one_outer_iteration()
        
        #Calculate Saturations
        v_total = sp.sum(self._net['pore.volume'])+sp.sum(self._net['throat.volume'])
        sat = 0.
        self['pore.inv_sat'] = 1.
        self['throat.inv_sat'] = 1.    
        for i in range(1,self._tseq+1):
            inv_pores = sp.where(self['pore.inv_seq']==i)[0]
            inv_throats = sp.where(self['throat.inv_seq']==i)[0]
            new_sat = (sum(self._pore_volumes[inv_pores])+sum(self._throat_volumes[inv_throats]))/v_total
            sat += new_sat
            self['pore.inv_sat'][inv_pores] = sat
            self['throat.inv_sat'][inv_throats] = sat
        self.sat = sat

    def _do_one_outer_iteration(self):
        r"""
        One iteration of an outer iteration loop for an algorithm
        (e.g. time or parametric study)
        """
        if (sp.mod(self._counter,500)==False):
            self._logger.info("Outer Iteration (counter = "+str(self._counter)+")")
        self._do_inner_iteration_stage()
        self._condition_update()
        self._counter += 1

    def _do_inner_iteration_stage(self):
        r"""
        Executes the inner iteration stage
        """
        self._logger.debug("  Inner Iteration Stage: ")

        self._plast = len(np.nonzero(self['pore.cluster_final'])[0])
        if self._timing:
            # determine the cluster with the earliest Haines time
            self._current_cluster = 1 + self._cluster_data['haines_time'].tolist().index(min(self._cluster_data['haines_time']))
            # update simulation clock
            self._logger.debug( 'sim time = ')
            self._logger.debug(self._sim_time)
            self._logger.debug(' haines time:')
            self._logger.debug( self._cluster_data['haines_time'])
            # The code really messes up when the [0] isn't in the next line. sim_time seems to just point to a place on the haines time array
            self._sim_time = min(self._cluster_data['haines_time'])
            self._logger.debug( 'sim time after update= ')
            self._logger.debug(self._sim_time)
        else:
            # Cycle to the next active cluster
            condition = 0
            loop_count = 0
            original_cluster = self._current_cluster
            cnum = original_cluster+1
            while condition == 0:
                if cnum > self._clusterCount:
                    cnum = 1
                if self._cluster_data['active'][cnum-1] == 1:
                    condition = 1
                    self._current_cluster = cnum
                if cnum == original_cluster:
                    loop_count = loop_count+1
                if loop_count > 1:
                    self._logger.error('No clusters active. Stuck in infinite loop.')
                cnum = cnum + 1

        # run through the Haines Jump steps
        self._do_one_inner_iteration()
        self._pnew = len(np.nonzero(self['pore.cluster_final'])[0])
        self._tseq += 1
        if self._pnew>self._plast:
            self._pseq += 1


    def _do_one_inner_iteration(self):
        r"""
        Executes one inner iteration
        """
        self._logger.debug("    Inner Iteration")
        # Fill throat and connecting pore
        # Pop out the largest throat (lowest inv_Pc) in the list, read the throat number
        tinvade = heapq.heappop(self._tpoints[self._current_cluster-1])[1]
        emptyCluster = -1
        fullCluster =  self._current_cluster
        if self._tpoints[self._current_cluster-1] == []:
            emptyCluster = self._current_cluster
        self._logger.debug( ' ')
        self._logger.debug( '--------------------------------------------------')
        self._logger.debug( 'STEP')
        self._logger.debug(self._tseq)
        self._logger.debug( 'trying to access cluster: ')
        self._logger.debug(self._current_cluster)
        self._logger.debug( 'when these clusters are active active: ')
        self._logger.debug(sp.nonzero(self._cluster_data['active'])[0])
        self._logger.debug( 'Haines at throat,time: ')
        self._logger.debug(tinvade)
        if self._timing:
            self._logger.debug(self._sim_time)

        # Mark throat as invaded
        self['throat.inv_seq'][tinvade] = self._tseq
        if self._timing:
            self['throat.inv_time'][tinvade] = self._sim_time
            # update self._cluster_data.['pore_volume']
            self._cluster_data['throat_volume'][self._current_cluster-1] += self._throat_volumes[tinvade]
            # Remove throat's contribution to the vol_coef
            self._cluster_data['vol_coef'][self._current_cluster-1] = self._cluster_data['vol_coef'][self._current_cluster-1] - self._Tvol_coef[tinvade]
        # Mark pore as invaded
        Pores = self._net.find_connected_pores(tinvade)
        # If both pores are already invaded:
        if np.in1d(Pores,np.nonzero(self['pore.cluster_final'])[0]).all():
            self._NewPore = -1
            # Label invaded throat with smaller cluster number
            #find cluster 1
            clusters = self._cluster_data['transform'][self['pore.cluster_final'][Pores]-1]
            self._logger.debug('clusters = ')
            self._logger.debug(clusters)
            self._current_cluster = min(clusters)
            self['throat.cluster_final'][tinvade] = self._current_cluster
            # if pores are from 2 different clusters:
            if self['pore.cluster_final'][Pores[0]]!=self['pore.cluster_final'][Pores[1]] :
                # find name of larger cluster number
                maxCluster = max(clusters)
                curCluster = self._current_cluster
                if emptyCluster == maxCluster:
                    fullCluster = curCluster
                if emptyCluster == curCluster:
                    fullCluster = maxCluster
                self._logger.info(' ')
                self._logger.info('CLUSTERS COMBINING:')
                self._logger.info(curCluster)
                self._logger.info(maxCluster)
                if self._timing:
                    self._logger.info('at time')
                    self._logger.info(self._sim_time)
                # update the cluster transform
                self._cluster_data['transform'][self._cluster_data['transform']==maxCluster] = [curCluster][0]
                # check if either was inactive (broke through already)
                if self._cluster_data['active'][maxCluster-1] + self._cluster_data['active'][self._current_cluster-1]<2:
                    self._logger.debug('making clusters ')
                    self._logger.debug(self._current_cluster)
                    self._logger.debug('and')
                    self._logger.debug(maxCluster)
                    self._logger.debug('inactive due to one being inactive already')
                    self._logger.debug(self._cluster_data['active'][curCluster-1])
                    self._logger.debug(self._cluster_data['active'][maxCluster-1])
                    self.cluster_remove(curCluster)
                    self._logger.info(' ')
                    self._logger.info('CLUSTER MERGED WITH A BREAKTHROUGH CLUSTER')
                else:
                    # relabel all pores and throats from larger number with smaller number
                    cluster_pores = self.toindices((self['pore.cluster_final']==maxCluster) + (self['pore.cluster_final']==curCluster))
                    cluster_throats = self.toindices((self['throat.cluster_final']==maxCluster) + (self['throat.cluster_final']==curCluster))
                    if emptyCluster == -1: 
                        cluster_int_throats = list(zip(*self._tpoints[curCluster-1]))[1] + list(zip(*self._tpoints[maxCluster-1]))[1]
                    else: 
                        cluster_int_throats = list(zip(*self._tpoints[fullCluster-1]))[1]
                    self._cluster_data['flow_rate'][curCluster-1] += self._cluster_data['flow_rate'][maxCluster-1]
                    self.cluster_update(curCluster,cluster_pores,cluster_throats,cluster_int_throats,tinvade)
                self._logger.info('making cluster ')
                self._logger.info(maxCluster)
                self._logger.info('inactive due to merge')
                # update the old cluster's activity and time
                self.cluster_remove(maxCluster)                


        else:
            # label invaded throat with current cluster
            self['throat.cluster_final'][tinvade] = self._current_cluster
            # find univaded pore, NewPore
            self._NewPore = Pores[self['pore.cluster_final'][Pores]==0][0]
            self._logger.debug( ' ')
            self._logger.debug( 'INVADING PORE: ')
            self._logger.debug(self._NewPore)
            self._logger.debug('the other pore is one of: ')
            self._logger.debug(Pores)
            self._logger.debug( 'position: ')
            self._logger.debug(self._net.get_data(prop='coords',pores='all')[self._NewPore])
            # label that pore as invaded
            self['pore.cluster_final'][self._NewPore] = self._current_cluster
            self['pore.cluster_original'][self._NewPore] = self._current_cluster
            if self._timing:
                self['pore.inv_time'][self._NewPore] = self._sim_time
            self['pore.inv_seq'][self._NewPore] = self._tseq
            if self._timing:
                # update self._cluster_data.['pore_volume']
                self._cluster_data['pore_volume'][self._current_cluster-1] += self._pore_volumes[self._NewPore]
            # Make a list of all throats neighboring pores in the cluster
            # Update interface list
            neighbors = self._net.find_neighbor_throats(self._NewPore)
            for j in neighbors:
                # If a throat is not labelled as invaded by the cluster, it must be an interfacial throat
                if (j not in self._tlists[self._current_cluster-1]):
                    self._logger.debug( 'new throat:')
                    self._logger.debug(j)
                    self._logger.debug('connecting pores:')
                    self._logger.debug(self._net.find_connected_pores(j))
                    # Add this throat data (pressure, number) to this cluster's "heap" of throat data.
                    heapq.heappush(self._tpoints[self._current_cluster-1],(self._phase.get_data(prop=self._capillary_pressure_name,throats='all')[j],j))
                    # Add new throat number to throat list for this cluster
                    self._tlists[self._current_cluster-1].append(j)
                    if self._timing:
                        # Update the cluster's vol_coef
                        self._cluster_data['vol_coef'][self._current_cluster-1] = self._cluster_data['vol_coef'][self._current_cluster-1]+self._Tvol_coef[j]
        if self._tpoints[self._current_cluster-1] != []:
            # Make sure you are not re-invading a throat in the next step (might never happen with new cluster routines)
            while self['throat.cluster_final'][self._tpoints[self._current_cluster-1][0][1]] > 0:
                tremove = heapq.heappop(self._tpoints[self._current_cluster-1])[1]
                if self._tpoints[self._current_cluster-1] == []:
                    self._logger.debug( 'making cluster ')
                    self._logger.debug(self._current_cluster)
                    self._logger.debug('inactive due to tpoints = [] ')
                    self.cluster_remove(self._current_cluster)
                    print('still happening!')
                    break
            # Find next Haines Jump info
            if self._tpoints[self._current_cluster-1] != []:
                next_throat = self._tpoints[self._current_cluster-1][0][1]
                self._cluster_data['haines_throat'][self._current_cluster-1] = next_throat
                if self._timing:
                    self._cluster_data['haines_pressure'][self._current_cluster-1] = self._tpoints[self._current_cluster-1][0][0]
                    self._cluster_data['cap_volume'][self._current_cluster-1] = self._cluster_data['haines_pressure'][self._current_cluster-1]*self._cluster_data['vol_coef'][self._current_cluster-1]
                    # Calculate the new Haines jump time
                    self._logger.debug( 'haines time before last stage:')
                    self._logger.debug( self._cluster_data['haines_time'])
        if self._tpoints[self._current_cluster-1] == []:
            self._logger.debug('making cluster ')
            self._logger.debug(self._current_cluster)
            self._logger.debug('inactive due to self._tpoints being empty for that cluster')
            self.cluster_remove(self._current_cluster)
        if self._timing:
            if self._cluster_data['active'][self._current_cluster-1] == 1:
                self._cluster_data['haines_time'][self._current_cluster-1] = (self._cluster_data['pore_volume'][self._current_cluster-1]+self._cluster_data['throat_volume'][self._current_cluster-1]+self._cluster_data['cap_volume'][self._current_cluster-1])/self._cluster_data['flow_rate'][self._current_cluster-1]
            if self._cluster_data['haines_time'][self._current_cluster-1] < self._sim_time:
                self._cluster_data['haines_time'][self._current_cluster-1] = self._sim_time
            self._logger.debug('haines time at the end of the throat stuff')
            self._logger.debug(self._cluster_data['haines_time'])

    def _condition_update(self):
         # Calculate the distance between the new pore and outlet pores
        if self._end_condition == 'breakthrough':
            newpore_position = self._net.get_data(prop='coords',pores='all')[self._NewPore]
            dist_sqrd = (self._outlet_position-newpore_position)*(self._outlet_position-newpore_position)
            if dist_sqrd[0].shape==(3,):     # need to do this for MatFile networks because newpore_position is a nested array, not a vector (?)
                dist_sqrd = dist_sqrd[0]
            newpore_distance = np.sqrt(dist_sqrd[0]+dist_sqrd[1]+dist_sqrd[2])
            self._logger.debug( 'newpore distance')
            self._logger.debug( newpore_distance)
            if newpore_distance < self._current_distance:
                self._percent_complete = np.round((self._initial_distance-newpore_distance)/self._initial_distance*100, decimals = 1)
                self._logger.info( 'percent complete')
                self._logger.info( self._percent_complete)
                self._current_distance = newpore_distance
        elif self._end_condition == 'total':
            self._percent_complete = np.round((np.sum(self['pore.cluster_final']>0)/self._net.num_pores())*100, decimals = 1)
        if self._percent_complete > self._rough_complete + self._rough_increment:
            self._rough_complete = np.floor(self._percent_complete/self._rough_increment)*self._rough_increment
            print('     IP algorithm at',np.int(self._rough_complete),'% completion at',np.round(misc.toc(quiet=True)),'seconds')

        # Determine if a new breakthrough position has occured
        if self._NewPore in self._outlets:
            self._logger.info( ' ')
            self._logger.info( 'BREAKTHROUGH AT PORE: ')
            self._logger.info(self._NewPore)
            self._logger.info('in cluster ')
            self._logger.info(self._current_cluster)
            if self._timing:
                self._logger.info('at time')
                self._logger.info(self._sim_time)
            if self._end_condition == 'breakthrough':
                self.cluster_remove(self._current_cluster)
            elif self._end_condition == 'total':
                self._brkevent.append(self._NewPore)
        if np.sum(self._cluster_data['active']) == 0:
            self._logger.info( ' ')
            self._logger.info( 'SIMULATION FINISHED; no more active clusters')
            if self._timing:
                self._logger.info('at time')
                self._logger.info(self._sim_time)
            self._condition = 0
            print('     IP algorithm at 100% completion at ',np.round(misc.toc(quiet=True)),' seconds')
            
    def cluster_update(self,cl_num,pores,throats,int_throats,bad_throat=-1):
        r'''
        '''
        int_throats = sp.unique(int_throats)
        int_throats = int_throats[int_throats!=bad_throat]
        pores = sp.unique(pores)
        throats = sp.unique(throats)
        #label all pores as invaded
        self['pore.cluster_final'][pores] = cl_num
        if sp.shape(throats) != (0,):
            self['throat.cluster_final'][throats] = cl_num
        if self._timing:
            # Calculate total volume in all invaded pores
            self._cluster_data['pore_volume'][cl_num-1] = np.sum(self._pore_volumes[pores])
            # Calculate total volume in all invaded throats
            if sp.shape(throats) != (0,):
                self._cluster_data['throat_volume'][cl_num-1] = np.sum(self._throat_volumes[throats])
            # Sum all interfacial throats' volume coeffients for throat cap volume calculation
            self._cluster_data['vol_coef'][cl_num-1] = np.sum(self._Tvol_coef[int_throats])
        # Make a list of all entry pressures of the interfacial throats
        interface_throat_pressures = self['throat.inv_Pc'][int_throats]#[0]
        # Zip pressures and numbers together so that HeapQ can work its magic
        Interface= list(zip(interface_throat_pressures,int_throats))
        # Turn the zipped throat interfaces object into a heap
        heapq.heapify(Interface)
        # Add to the total list of interface throats in the system
        self._tlists[cl_num-1] = int_throats.tolist()
        # Add to the total list of invaded interface throats in the system
        self._tpoints[cl_num-1] = Interface
        # Pop off the first entry (lowest pressure) on the throat info list
        invaded_throat_info = Interface[0]
        if self._timing:
            # Determine pressure at Haines Jump
            self._cluster_data['haines_pressure'][cl_num-1] = invaded_throat_info[0]
            # Calculate cap_volume at Haines Jump
            self._cluster_data['cap_volume'][cl_num-1] = self._cluster_data['haines_pressure'][cl_num-1]*self._cluster_data['vol_coef'][cl_num-1]
            # Calculate throat_volume at Haines Jump
            self._cluster_data['throat_volume'][cl_num-1] = self._cluster_data['throat_volume'][cl_num-1]+self._throat_volumes[invaded_throat_info[1]]
            # Calculate time at Haines Jump
            self._cluster_data['haines_time'][cl_num-1] = (self._cluster_data['pore_volume'][cl_num-1]+self._cluster_data['throat_volume'][cl_num-1]+
                                        self._cluster_data['cap_volume'][cl_num-1])/self._cluster_data['flow_rate'][cl_num-1]
        # Record invaded throat
        self._cluster_data['haines_throat'][cl_num-1] = invaded_throat_info[1]
    
    def cluster_remove(self,cl_num):
        if self._timing:
            self._cluster_data['haines_time'][cl_num-1] = 1e32
        self._cluster_data['active'][cl_num-1] = 0
        self._tpoints[cl_num-1] = []
        
        
    def update_results(self,occupancy='occupancy',IPseq=None,IPsat=None):
        r"""
        
        Returns
        -------
        The invading phase will aquire the following pore data ::
        
            occupancy           : 0. for univaded, 1. for invaded
            IP_cluster_final    : 0 for uninvaded, merged cluster number for invaded
            IP_cluster_original : 0 for uninvaded, original cluster number for invaded
            IP_inv_seq          : 0 for uninvaded, simulation step for invaded
            IP_inv_time         : 0 for uninvaded, simulation time for invaded
        
        and throat data ::
        
            occupancy           : 0 for univaded, 1 for invaded
            IP_cluster_final    : 0 for uninvaded, merged cluster number for invaded
            IP_inv_seq          : 0 for uninvaded, simulation step for invaded
            IP_inv_time         : 0 for uninvaded, simulation time for invaded
        
        """
        self._phase['pore.IP_cluster_final']=self['pore.cluster_final']
        self._phase['pore.IP_cluster_original']=self['pore.cluster_original']
        self._phase['throat.IP_cluster_final']=self['throat.cluster_final']
        self._phase['pore.IP_inv_seq']=self['pore.inv_seq']
        self._phase['throat.IP_inv_seq']=self['throat.inv_seq']
        if self._timing:
            self._phase['pore.IP_inv_time']=self['pore.inv_time']
            self._phase['throat.IP_inv_time']=self['throat.inv_time']
            
        if IPseq==None:
            if IPsat != None:
                sat_pores = self['pore.inv_sat']<=IPsat
                if sum(sat_pores) == 0:
                    IPseq = 0
                else:
                    IPseq = max(self['pore.inv_seq'][sat_pores])
            else:
                IPseq = self._tseq

        try:
            self._phase['pore.'+occupancy] = 0.
            inv_pores = (self['pore.inv_seq']>0)&(self['pore.inv_seq']<=IPseq)
            self._phase['pore.'+occupancy][inv_pores] = 1.
            self['pore.invaded'] = inv_pores
            self._phase['throat.'+occupancy] = 0.
            inv_throats = (self['throat.inv_seq']>0)&(self['throat.inv_seq']<=IPseq)
            self._phase['throat.'+occupancy][inv_throats] = 1.
            self['throat.invaded'] = inv_throats
            self.sat = max(self['throat.inv_sat'][inv_throats])
            
        except:
            print('Something bad happened while trying to update phase',self._phase.name)
        try:
            self._phase_def['pore.'+occupancy]=~inv_pores
            self['pore.defended'] = ~inv_pores
            self._phase_def['throat.'+occupancy]=~inv_throats
            self['throat.defended'] = ~inv_throats
        except:
            print('A partner phase has not been set so inverse occupancy cannot be set')

        
if __name__ == '__main__':
    import doctest
    doctest.testmod(verbose=True)
    