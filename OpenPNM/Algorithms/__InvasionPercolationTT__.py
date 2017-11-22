# -*- coding: utf-8 -*-
"""
===============================================================================
InvasionPercolationBasic: Simple IP
===============================================================================

"""
import heapq as hq
import scipy as sp
import numpy as np
from OpenPNM.Algorithms import GenericAlgorithm
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import OpenPNM.Utilities.vertexops as vo
import time


class InvasionPercolationTT(GenericAlgorithm):
    r"""
    A classic/basic invasion percolation algorithm optimized for speed.

    Parameters
    ----------
    network : OpenPNM Network object
        The Network upon which the invasion should occur.

    Notes
    ----
    n/a

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def setup(self, phase, def_phase, **kwargs):
        r"""
        Set up the required parameters for the algorithm

        Parameters
        ----------
        phase : OpenPNM Phase object
            The phase to be injected into the Network.  The Phase must have the
            capillary entry pressure values for the system.

        prop : string
            The name of the property containing the capillary entry
            pressure.  This should be applied to both throats and pores.

        """
        self._phase = phase
        self._def_phase = def_phase
        # Setup arrays and info
        if np.shape(np.shape(phase['throat.capillary_pressure']))[0] > 1:
            self._bi_directional = True
            #If bi-directional throats, get the first one for now, switched
            #later when algorithm runs and works out which one to apply
            self['throat.entry_pressure'] = phase['throat.capillary_pressure'][:,0]
        else:
            self['throat.entry_pressure'] = phase['throat.capillary_pressure']
            self._bi_directional = False
        self['pore.entry_pressure'] = phase['pore.capillary_pressure']
        self.reset_invasion_info()
        self._key_words = kwargs
   
    def reset_invasion_info(self):
        for elem in ['pore', 'throat']:
            for prop in ['inv_Pc',
                         'inv_sat',
                         'occupancy',
                         'action']:
                try:
                    del self[elem+'.'+prop]
                except:
                    pass
            for prop in ['inv_seq',
                         'inv_Pc',
                         'cluster',
                         'action']:
                try:
                    self[elem+'.'+prop] = -1
                except:
                    pass
        # Masks for tracking pores and throats at the interface
        # Basically a quick way of getting to all the elements in the queues
        self._interface_Ts = np.zeros(self.Nt, dtype=bool)
        self._interface_Ps = np.zeros(self.Np, dtype=bool)

    def set_inlets(self, pores=None, clusters=None):
        r"""

        Parameters
        ----------
        pores : array_like
            The list of inlet pores from which the Phase can enter the Network
        """
        if 'inlets' in self._key_words.keys():
            pores = self._key_words['inlets']
        
        if clusters is None:
            clusters = []
            clusters.append(pores)
        self.queue={}
        try:
            inlet_inv_seq = self._key_words['inlet_inv_seq']
        except:
            inlet_inv_seq = -1
        for i, cluster in enumerate(clusters):
            self.queue[i]=[]
            # Perform initial analysis on input pores
            self['pore.inv_seq'][cluster] = inlet_inv_seq
            self['pore.cluster'][cluster] = i
            if np.size(cluster) > 1:
                for elem_id in cluster:
                    self._add_ts2q(elem_id, self.queue[i], action=0)
            elif np.size(cluster) == 1:
                self._add_ts2q(cluster, self.queue[i], action=0)
            else:
                logger.warning("Some inlet clusters have no pores")
        try:
            if self._key_words['snap_off']:
                self.apply_snap_off()
        except:
            pass
        try:
            if self._key_words['partial']:
                self.apply_partial_sat()
        except:
            self._phase['pore.occuancy']=False
            self._phase['throat.occupancy']=False
            self._def_phase['pore.occupancy']=True
            self._def_phase['throat.occupancy']=True

    
    def _add_ts2q(self, pore, queue, action=-1):
        """
        Helper method to add throats to the queue
        """
        elem_type = 'throat'
        # Find throats connected to newly invaded pore
        Ts = self._net.find_neighbor_throats(pores=pore)
        # Remove already invaded throats from Ts
        Ts = Ts[self['throat.inv_seq'][Ts] <= 0]
        if len(Ts) > 0:
            self._interface_Ts[Ts]=True
            for T in Ts:
                if self._bi_directional:
                    #new code for bi-directional Pc
                    #get index of pore being invaded next and apply correct entry pressure
                    pind = list(self._net['throat.conns'][T] != pore).index(True)
                    self['throat.entry_pressure'][T] = self._phase['throat.capillary_pressure'][T][pind]
                data = []
                # Pc
                data.append(self['throat.entry_pressure'][T])
                # Element Index
                data.append(T)
                # Element Type (Pore of Throat)
                data.append(elem_type)
                # Invasion Action - 0=Burst, 1=Coop, 2=Snap, 3=Touch??
                data.append(action)
                hq.heappush(queue, data)
    
    def _add_ps2q(self, throat, queue, action=-1):
        """
        Helper method to add pores to the queue
        """
        elem_type = 'pore'
        # Find pores connected to newly invaded throat
        Ps = self._net['throat.conns'][throat]
        # Remove already invaded pores from Ps
        Ps = Ps[self['pore.inv_seq'][Ps] <= 0]
        if len(Ps) > 0:
            self._interface_Ps[Ps]=True
            for P in Ps:
                data = []
                # Pc
                data.append(self["pore.entry_pressure"][P])
                # Element Index
                data.append(P)
                # Element Type (Pore of Throat)
                data.append(elem_type)
                # Invasion Action - 0=Burst, 1=Coop, 2=Snap, 3=Touch??
                data.append(action)
                hq.heappush(queue, data)
        
    def run(self, max_pressure=None, outlets=None, **kwargs):
        r"""
        Perform the algorithm

        Parameters
        ----------
        n_steps : int
            The number of throats to invaded during this step

        """
        if 'throat.entry_pressure' not in self.keys():
            logger.error("Setup method must be run first")
        #if sp.all(self['pore.inv_seq'] == -1):
        if 'inlets' in self._key_words.keys():
            logger.info("Setting inlet pores at shared pressure")
            self.set_inlets(pores=self._key_words['inlets'])
        elif 'clusters' in self._key_words.keys():
            logger.info("Setting inlet clusters at individual pressures")
            self.set_inlets(clusters=self._key_words['clusters'])
        else:
            logger.error("Either 'inlets' or 'clusters' must be passed to"+
                         " setup method")
        self._coop_fill = False
        if 'coop_fill' in self._key_words.keys():
            if self._key_words['coop_fill']:
                self._coop_fill = True
                try:
                    inv_points = kwargs['inv_points']
                except:
                    inv_points = None
                self._setup_coop_filling(filling_model='throat.alpha', inv_points=inv_points)
        if max_pressure is None:
            max_pressure = sp.inf
        #queue = self.queue
        if len(self.queue.items()) == 0:
            logger.warn('queue is empty, this network is fully invaded')
            return

        max_p_reached = [False]*len(self.queue.items())
        count = -1
        invasion_running = [True]*len(self.queue.items())
        high_Pc = np.zeros(len(self.queue.items()))
        while np.any(invasion_running) and not np.all(max_p_reached):

            for cluster_num in self.queue.keys():
                if invasion_running[cluster_num]:
                    queue = self.queue[cluster_num]
                    pressure, elem_id, elem_type, action = hq.heappop(queue)
                    if elem_type == 'pore':
                        self._interface_Ps[elem_id]=False
                    else:
                        self._interface_Ts[elem_id]=False
                    if pressure > max_pressure:
                        max_p_reached[cluster_num]=True
                    else:
                        if self[elem_type+'.inv_seq'][elem_id] == -1:
                            count += 1
                            #record highest Pc cluster has reached
                            if high_Pc[cluster_num]< pressure:
                                high_Pc[cluster_num] = pressure
                            "The newly invaded element is available for invasion"
                            self[elem_type+'.inv_seq'][elem_id] = count
                            self[elem_type+'.cluster'][elem_id] = cluster_num
                            self[elem_type+'.inv_Pc'][elem_id] = high_Pc[cluster_num]
                            self[elem_type+'.action'][elem_id] = action
                            if elem_type == 'throat':
                                self._add_ps2q(elem_id, queue, action=0)
                            elif elem_type == 'pore':
                                self._add_ts2q(elem_id, queue, action=0)
                                if self._coop_fill:
                                    self._check_coop(elem_id, queue)
                            
                        elif (self[elem_type+'.cluster'][elem_id] != cluster_num and 
                              invasion_running[self[elem_type+'.cluster'][elem_id].astype(int)]):
                            "The newly invaded element is part of an invading cluster"
                            "Merge the clusters using the existing cluster number"
                            existing_cluster = self[elem_type+'.cluster'][elem_id]
                            logger.info("Merging cluster "+str(cluster_num)+
                                        " into cluster "+str(existing_cluster)+
                                        " at sequence "+str(count))

                        elif (self[elem_type+'.cluster'][elem_id] != cluster_num and 
                              ~invasion_running[self[elem_type+'.cluster'][elem_id].astype(int)]):
                            "The newly invaded element is part of cluster that has stopped invading"
                            logger.info("Cluster " +str(cluster_num)+ " terminated")

                        elif self[elem_type+'.cluster'][elem_id]==cluster_num:
                            "Self intersecting or repeating elements"
                            pass
                        else:
                            logger.warning("Clusters " +str(cluster_num)+" and "+
                                        str(self[elem_type+'.cluster'][elem_id])+
                                        " did something funny")

                    if len(queue) == 0 or max_p_reached[cluster_num]:
                        # If the cluster contains no more entries invasion has finished
                        invasion_running[cluster_num]=False
            self._invade_isolated_Ts()
            if outlets is not None:
                tcs = np.unique(self['pore.cluster'][outlets]).astype(int) # terminated clusters
                tcs = tcs[tcs >= 0]
                if len(tcs) > 0:
                    for tc in tcs:
                        if invasion_running[tc] == True:
                            invasion_running[tc]=False
                            logger.info("Cluster " +str(tc)+ " reached outlet at sequence "+str(count))
        

    def return_results(self, pores=[], throats=[], Pc=None):
        r"""
        Places the results of the IP simulation into the Phase object.

        Parameters
        ----------
        pores and throats : array_like
            The list of pores and throats whose values should be returned to
            the Phase object.  Default is all of them.

        """
        pores = sp.array(pores, ndmin=1)
        throats = sp.array(throats, ndmin=1)
        if len(pores) == 0:
            pores = self.Ps
        if len(throats) == 0:
            throats = self.Ts
        self._phase['throat.inv_seq'] = sp.nan
        self._phase['pore.inv_seq'] = sp.nan
        self._phase['throat.inv_seq'][throats] = \
            self['throat.inv_seq'][throats]
        self._phase['pore.inv_seq'][pores] = \
            self['pore.inv_seq'][pores]
        self._phase['throat.cluster'] = self['throat.cluster']
        self._phase['pore.cluster'] = self['pore.cluster']
        if Pc is not None:
            if 'pore.inv_Pc' in self.keys():
                self['throat.occupancy'] = self['throat.inv_Pc'] <= Pc
                self['pore.occupancy'] = self['pore.inv_Pc'] <= Pc
                self._phase['throat.occupancy'] = self['throat.occupancy'].astype(float)
                self._phase['pore.occupancy'] = self['pore.occupancy'].astype(float)
                self._def_phase['throat.occupancy'] = 1.0 - self._phase['throat.occupancy']
                self._def_phase['pore.occupancy'] = 1.0 - self._phase['pore.occupancy']
            else:
                logger.warning("Occupancy not updated, please run " +
                               "extract_drainage() method to populate" +
                               " the inv_Pc data")
        if "pore.inv_Pc" in self.props():
            self._phase['pore.inv_Pc'] = \
                self['pore.inv_Pc']
            self._phase['throat.inv_Pc'] = \
                self['throat.inv_Pc']
        if "pore.inv_sat" in self.props():
            self._phase['pore.inv_sat'] = \
                self['pore.inv_sat']
            self._phase['throat.inv_sat'] = \
                self['throat.inv_sat']
        if "pore.trapped" in self.labels():
            self._phase['pore.trapped'] = self['pore.trapped']
        if "pore.trapping_sequence" in self.props():
            self._phase['pore.trapping_sequence'] = \
                self['pore.trapping_sequence']
        if "throat.trapped" in self.labels():
            self._phase['throat.trapped'] = self['throat.trapped']
        if "throat.trapping_sequence" in self.props():
            self._phase['throat.trapping_sequence'] = \
                self['throat.trapping_sequence']

    def apply_flow(self, flowrate):
        r"""
        Convert the invaded sequence into an invaded time for a given flow rate
        considering the volume of invaded pores and throats.

        Parameters
        ----------
        flowrate : float
            The flow rate of the injected fluid

        Returns
        -------
        Creates a throat array called 'invasion_time' in the Algorithm
        dictionary

        """
        P12 = self._net['throat.conns']
        a = self['throat.inv_seq']
        b = sp.argsort(self['throat.inv_seq'])
        P12_inv = self['pore.inv_seq'][P12]
        # Find if the connected pores were invaded with or before each throat
        P1_inv = P12_inv[:, 0] == a
        P2_inv = P12_inv[:, 1] == a
        c = sp.column_stack((P1_inv, P2_inv))
        d = sp.sum(c, axis=1, dtype=bool)  # List of Pores invaded with each throat
        # Find volume of these pores
        P12_vol = sp.zeros((self.Nt,))
        P12_vol[d] = self._net['pore.volume'][P12[c]]
        # Add invaded throat volume to pore volume (if invaded)
        T_vol = P12_vol + self._net['throat.volume']
        # Cumulative sum on the sorted throats gives cumulated inject volume
        e = sp.cumsum(T_vol[b] / flowrate)
        t = sp.zeros((self.Nt,))
        t[b] = e  # Convert back to original order
        self._phase['throat.invasion_time'] = t
    
    def apply_pc(self, Pc):
        r"""
        Invasion sequence runs concurrently for all clusters so must assess
        each one individually to determine the pressure threshold condition.
        N.B This may not work as clusters could have merged before pressure
        condition is met...
        """
        # If occupancy data not found on object initialise it ready for
        # filling with cluster data that will need to index into it
        if 'throat.occupancy' not in self.keys():
            self['throat.occupancy'] = False
        if 'pore.occupancy' not in self.keys():
            self['pore.occupancy'] = False
        for cluster in np.unique(self['pore.cluster']):
            self._apply_pc_to_cluster(Pc, cluster)
            
    def _apply_pc_to_cluster(self, Pc, cluster):
        r"""
        Inspect the invasion sequence and apply a capillary pressure.
        Update the fluid occupancy based on invading up to this pressure.
        Effectively discount any pores or throats with invasion sequence
        after the first invaded element with pressure exceeding this value
        including that element.
        """
        test_pressure = Pc
        Ps = self['pore.cluster'] == cluster
        Ts = self['throat.cluster'] == cluster
        #ignore totally uninvaded
        Ps *= self['pore.inv_seq']>=0
        Ts *= self['throat.inv_seq']>=0
        if (np.sum(Ts) + np.sum(Ps)) > 0:
            # Filter univaded pores and throats for this cluster only
            uninvaded_ps = self['pore.inv_seq'][Ps][self['pore.entry_pressure'][Ps] > test_pressure]
            uninvaded_ts = self['throat.inv_seq'][Ts][self['throat.entry_pressure'][Ts] > test_pressure]
            # Assess whether the cluster is partially invaded at this pressure
            # Break points in the pore and throat invasion sequence will exist
            # Otherwise
            break_p = []
            break_t = []
            #remove inlets and those that have still not been invaded by running the algorithm
            #i.e. those with sequence -1 may not have been reached if max steps or max pressure was reached
            if np.size(uninvaded_ps) > 0:
                break_p = uninvaded_ps[uninvaded_ps>0]
            if np.size(uninvaded_ts) > 0:   
                break_t = uninvaded_ts[uninvaded_ts>0]
    
            if (np.size(break_p) > 0) or (np.size(break_t) > 0):
            #if np.size(break_t) > 0:
                if (np.size(break_p) == 0):
                    break_seq = break_t.min()
                elif (np.size(break_t) == 0):
                    break_seq = break_p.min()
                else:
                    break_seq = np.min((break_p.min(),break_t.min()))
                #print(Pc,break_seq)
                self['throat.occupancy'][Ts] = (self['throat.inv_seq'][Ts] < break_seq)
                self['pore.occupancy'][Ps] = (self['pore.inv_seq'][Ps] < break_seq)
            else:
                # Applied pressure exceeds all entry pressures in network
                if np.sum(Ts) > 0:
                    self['throat.occupancy'][Ts] = True
                    self['throat.occupancy'][Ts][self['throat.inv_seq'][Ts] < 0] = False

                if np.sum(Ps) > 0:
                    self['pore.occupancy'][Ps] = True
                    self['pore.occupancy'][Ps][self['pore.inv_seq'][Ps] < 0] = False
        else:
            logger.warning("Cluster "+str(cluster)+" has no invaded elements at Pressure: "+str(Pc))

    def extract_drainage(self, inv_points=None, npts=100, trapping_outlets=None):
        r"""
        Extract quasi-static drainage data by applying a set of capillary pressures
        """
        if "pore.inv_seq" not in self.props():
            logger.error("Cannot plot drainage curve. Please run algorithm first")
        if inv_points is None:
            t_entry_pc = self['throat.entry_pressure'][self['throat.inv_seq']>=0]
            p_entry_pc = self['pore.entry_pressure'][self['pore.inv_seq']>=0]
            pmin = np.min((t_entry_pc.min(),p_entry_pc.min()))
            pmax = np.max((t_entry_pc.max(),p_entry_pc.max()))
            inv_points = np.linspace(pmin,pmax,npts)
        inv_Pc_p = np.ones(self._net.Np)*sp.inf
        inv_Pc_t = np.ones(self._net.Nt)*sp.inf
        inv_Sat_p = np.ones(self._net.Np)*sp.inf
        inv_Sat_t = np.ones(self._net.Nt)*sp.inf
        pvol = np.sum(self._net['pore.volume'])
        tvol = np.sum(self._net['throat.volume'])
        tot_vol = pvol + tvol
        saturation = np.zeros(len(inv_points))
        for i, Pc in enumerate(inv_points):
            
            print("Extracting Drainage at Pc: "+str(Pc))
            ### Method 1
            #self.reset_invasion_info()
            #self.run(max_pressure=Pc)
            #self.check_snap_off(Pc)
            #self['pore.occupancy'] = self["pore.inv_seq"] > -1
            #self['throat.occupancy'] = self["throat.inv_seq"] > -1
            ### Method 2
            #self.apply_pc(Pc)
            #if trapping_outlets is not None:
            #    self.apply_trapping(outlets=trapping_outlets, partial=True)
            ### Method 3 taken from monitoring the highest pressure cluster
            ### has passed through at point of invasion
            self['pore.occupancy'] = self['pore.inv_Pc'] <= Pc
            self['throat.occupancy'] = self['throat.inv_Pc'] <= Pc
            #self['pore.occupancy'][uninv_Ps] = False
            #self['throat.occupancy'][uninv_Ts] = False
            #record pressure at which element is first invaded
            inv_Pc_p[sp.isinf(inv_Pc_p) * self['pore.occupancy']]=Pc
            inv_Pc_t[sp.isinf(inv_Pc_t) * self['throat.occupancy']]=Pc
            sat_p_vol = np.sum(self._net['pore.volume'][self['pore.occupancy']])
            sat_t_vol = np.sum(self._net['throat.volume'][self['throat.occupancy']])
            Sat = (sat_p_vol + sat_t_vol)/tot_vol
            saturation[i] = Sat
            inv_Sat_p[sp.isinf(inv_Sat_p) * self['pore.occupancy']]=Sat
            inv_Sat_t[sp.isinf(inv_Sat_t) * self['throat.occupancy']]=Sat

        self['pore.inv_Pc']=inv_Pc_p
        self['throat.inv_Pc']=inv_Pc_t
        self['pore.inv_sat']=inv_Sat_p
        self['throat.inv_sat']=inv_Sat_t
        return saturation

    
    def plot_drainage_curve(self, fig=None, inv_points=None, npts=100, lpf=False, trapping_outlets=None):
        r"""
        Plot a simple drainage curve
        """
        if "pore.inv_Pc" not in self.props():
            self.extract_drainage(inv_points, npts, trapping_outlets)
        if inv_points is None:
            inv_points = np.unique(self['throat.inv_Pc'][~sp.isnan(self['throat.inv_Pc'])])
        sat_p = np.zeros(len(inv_points))
        sat_t = np.zeros(len(inv_points))
        inv_p = self['pore.inv_Pc']
        inv_t = self['throat.inv_Pc']
        # Handle trapped pores and throats
        if np.sum(self['pore.inv_seq'] == -1) > 0:
            inv_p[self['pore.inv_seq'] == -1] = 2*inv_points.max()
        if np.sum(self['throat.inv_seq'] == -1) > 0:
            inv_t[self['throat.inv_seq'] == -1] = 2*inv_points.max()
        num_p = np.zeros(len(inv_points),dtype=int)
        num_t = np.zeros(len(inv_points),dtype=int)
        for i, Pc in enumerate(inv_points):
            if lpf:
                frac = self.evaluate_late_pore_filling(Pc, Swp_init=0.25, eta=2.5)
                p_vol = self._net['pore.volume']*frac
            else:
                p_vol = self._net['pore.volume']
            sat_p[i] = np.sum(p_vol[inv_p<=Pc])
            sat_t[i] = np.sum(self._net['throat.volume'][inv_t<=Pc])
            num_p[i] = np.sum(inv_p<=Pc)
            num_t[i] = np.sum(inv_t<=Pc)

        pvol = np.sum(self._net['pore.volume'])
        tvol = np.sum(self._net['throat.volume'])
        tot_vol = pvol + tvol
        tot_sat = sat_p + sat_t
        #Normalize
        sat_p /= tot_vol
        sat_t /= tot_vol
        tot_sat /= tot_vol
        if fig is None:
            fig = plt.figure()
        a = fig.add_subplot(111)
        a.plot(inv_points, sat_p, 'r*-',label='pore')
        a.plot(inv_points, sat_t, 'b*-',label='throat')
        a.plot(inv_points, tot_sat, 'g*-',label='total')
        a.legend(bbox_to_anchor=(0,1.02,1,0.102),loc=3, ncol=3, borderaxespad=0)
        a.set_xlabel("Capillary Pressure [Pa]")
        a.set_ylabel("Saturation")
        return fig, tot_sat, num_p, num_t


    def evaluate_late_pore_filling(self, Pc, Swp_init=0.75, eta=3.0,
                                   wetting_phase=False):
        r"""
        Compute the volume fraction of the phase in each pore given an initial
        wetting phase fraction (Swp_init) and a growth exponent (eta)
        returns the fraction of the pore volume occupied by wetting or
        non-wetting phase.
        Assumes Non-wetting phase displaces wetting phase
        """
        try:
            Swp = np.ones(len(self['pore.inv_Pc']))
            mask = self['pore.inv_Pc'] < Pc
            Swp[mask] = Swp_init*(self['pore.inv_Pc'][mask]/Pc)**eta
            Swp[np.isnan(Swp)] = 1.0
        except:
            Swp = np.ones(self._net.Np)
        Snwp = 1-Swp
        if wetting_phase:
            return Swp
        else:
            return Snwp

    def apply_trapping(self, outlets, partial=False):
        """
        Apply trapping based on algorithm described by Y. Masson [1].
        It is applied as a post-process and runs the percolation algorithm in
        reverse assessing the occupancy of pore neighbors.
        3 situations can happen on invasion without trapping:
            The number of defending clusters stays the same and clusters can
            shrink
            A cluster of size one is suppressed
            A cluster is split into multiple clusters
        In reverse the following situations can happen:
            The number of defending clusters stays the same and clusters can
            grow
            A cluster of size one is created
            Mutliple clusters merge into one cluster
        With trapping the reversed rules are adjusted so that:
            Only clusters that do not connect to a sink can grow and merge.
            At the point that a neighbor connected to a sink is touched the
            trapped cluster stops growing as this is the point of trapping in
            forward invasion time.

        Logger info displays the invasion Sequence and pore index and a message
        with condition number based on the modified trapping rules and the
        assignment of the pore to a given cluster.

        Initially all invaded pores are given cluster label -1
        Outlets / Sinks are given -2
        New clusters that grow into fully trapped clusters are either
        identified at the point of breakthrough or grow from nothing if the
        full invasion sequence is run, they are assigned numbers from 0 up.

        Ref:
        [1] Masson, Y., 2016. A fast two-step algorithm for invasion
        percolation with trapping. Computers & Geosciences, 90, pp.41-48

        Parameters
        ----------
        outlets : list or array of pore indices for defending fluid to escape
        through
        
        partial : Boolean indicating whether partially filled network

        Returns
        -------
        Creates a throat array called 'pore.clusters' in the Algorithm
        dictionary. Any positive number is a trapped cluster

        Also creates 2 boolean arrays Np and Nt long called '<element>.trapped'
        """

        if partial:
            # Set occupancy
            invaded_ps = self['pore.inv_seq'] > -1
            # Put defending phase into clusters
            clusters = self._net.find_clusters2(~invaded_ps)
            # Identify clusters that are connected to an outlet and set to -2
            # -1 is the invaded fluid
            # -2 is the defender fluid able to escape
            # All others now trapped clusters which grow as invasion is reversed
            out_clusters = sp.unique(clusters[outlets])
            for c in out_clusters:
                if c >= 0:
                    clusters[clusters == c] = -2
        else:
            # Go from end
            clusters = np.ones(self._net.Np, dtype=int)*-1
            clusters[outlets] = -2

        # Turn into a list for indexing
        inv_seq = np.vstack((self['pore.inv_seq'].astype(int),
                             np.arange(0, self._net.Np, dtype=int))).T
        # Reverse sort list
        inv_seq = inv_seq[inv_seq[:, 0].argsort()][::-1]
        next_cluster_num = np.max(clusters)+1
        # For all the steps after the inlets are set up to break-through
        # Reverse the sequence and assess the neighbors cluster state
        stopped_clusters = np.zeros(self._net.Np, dtype=bool)
        all_neighbors = self._net.find_neighbor_pores(self._net.pores(),
                                                      flatten=False)
        for un_seq, pore in inv_seq:
            if pore not in outlets and un_seq > -1:  # Don't include outlets
                nc = clusters[all_neighbors[pore]]  # Neighboring clusters
                unique_ns = np.unique(nc[nc != -1])  # Unique Neighbors
                seq_pore = "S:"+str(un_seq)+" P:"+str(pore)
                if np.all(nc == -1):
                    # This is the start of a new trapped cluster
                    clusters[pore] = next_cluster_num
                    next_cluster_num += 1
                    msg = (seq_pore+" C:1 new cluster number: " +
                           str(clusters[pore]))
                    logger.info(msg)
                elif len(unique_ns) == 1:
                    # Grow the only connected neighboring cluster
                    if not stopped_clusters[unique_ns[0]]:
                        clusters[pore] = unique_ns[0]
                        msg = (seq_pore+" C:2 joins cluster number: " +
                               str(clusters[pore]))
                        logger.info(msg)
                    else:
                        clusters[pore] = -2
                elif -2 in unique_ns:
                    # We have reached a sink neighbor, stop growing cluster
                    msg = (seq_pore+" C:3 joins sink cluster")
                    logger.info(msg)
                    clusters[pore] = -2
                    # Stop growth and merging
                    stopped_clusters[unique_ns[unique_ns > -1]] = True
                else:
                    # We might be able to do some merging
                    # Check if any stopped clusters are neighbors
                    if np.any(stopped_clusters[unique_ns]):
                        msg = (seq_pore+" C:4 joins sink cluster")
                        logger.info(msg)
                        clusters[pore] = -2
                        # Stop growing all neighboring clusters
                        stopped_clusters[unique_ns] = True
                    else:
                        # Merge multiple un-stopped trapped clusters
                        new_num = unique_ns[0]
                        clusters[pore] = new_num
                        for c in unique_ns:
                            clusters[clusters == c] = new_num
                            msg = (seq_pore + " C:5 merge clusters: " +
                                   str(c) + " into "+str(new_num))
                            logger.info(msg)

        # And now return clusters
        #self['pore.trap_cluster'] = clusters
        clusters[self._net.pores('boundary')] = -2
        num_trap = np.sum(np.unique(clusters) >= 0)
        if num_trap > 0:
            logger.info("Number of trapped clusters" + str(num_trap))
            self['pore.trapped'] = clusters > -1
            print("Number of trapped pores: ",np.sum(self['pore.trapped']))
            self['pore.inv_seq'][self['pore.trapped']] = -1
            self['throat.trapped'] = np.zeros([self._net.Nt], dtype=bool)
            for c in np.unique(clusters[clusters >= 0]):
                c_ts = self._net.find_neighbor_throats(clusters == c,
                                                       mode='intersection')
                self['throat.trapped'][c_ts] = True
            print("Number of trapped throats: ",np.sum(self['throat.trapped']))
            self['throat.inv_seq'][self['throat.trapped']] = -1
            #Assumes invasion has run to the end
            self._phase['pore.occupancy']= ~self['pore.trapped']
            self._def_phase['pore.occupancy']= self['pore.trapped']
            self._phase['throat.occupancy']= ~self['throat.trapped']
            self._def_phase['throat.occupancy']= self['throat.trapped']
        else:
            logger.info("No trapped clusters found")


    def apply_snap_off(self, snap_off='throat.snap_off', queue=None):
        r"""
        Add all the throats to the queue with snap off pressure
        """
        if queue is None:
            queue = self.queue[0]
        try:
            Pc_snap_off = self._phase[snap_off]
            logger.info("Adding snap off pressures to queue")
            for T in self._net.throats():
                if not np.isnan(Pc_snap_off[T]):
                    hq.heappush(queue, [Pc_snap_off[T], T, 'throat', 2])
        except:
            logger.warning("Phase "+self._phase.name+" doesn't have property "+
                           snap_off)
    
    def apply_partial_sat(self, queue=None):
        r"""
        Method to start invasion from a partially saturated state
        """
        if queue is None:
            invading_cluster = 0
            queue = self.queue[invading_cluster]
        occ_type = self._phase['pore.occupancy'].dtype
        occupied = np.array([1], dtype=occ_type)
        occ_Ps = self._phase['pore.occupancy'] == occupied
        occ_Ts = self._phase['throat.occupancy'] == occupied
        low_val = -999999
        if np.sum(occ_Ps) > 0:
            logger.warn("Applying partial saturation to "+str(np.sum(occ_Ps))+" pores")
            self['pore.inv_seq'][occ_Ps] = 0
            for P in self._net.pores()[occ_Ps]:
                self._add_ts2q(P, queue, action=0)
                self['pore.cluster'][P] = invading_cluster
                self['pore.inv_Pc'][P] = low_val
        if np.sum(occ_Ts) > 0:
            logger.warn("Applying partial saturation to "+str(np.sum(occ_Ts))+" throats")
        self['throat.inv_seq'][occ_Ts] = 0
        for T in self._net.throats()[occ_Ts]:
            #self._add_ps2q(T, queue)
            self['throat.cluster'][T] = invading_cluster
            self['throat.inv_Pc'][T] = low_val
    

    def _perpendicular_vector(self, vec):
        return np.cross(vec, [1,1,1])

    
    def _circ_points(self, a, b, c):
        t = np.arange(0, 2*np.pi, 0.1)
        pts = []
        for i in range(3):
            pts.append(c[i] + a[i]*np.sin(t) + b[i]*np.cos(t))
        return np.asarray(pts).T
    
    
    def _check_intersection(self, c1, c2, r1, r2, pore_center, pore_rad):
        r"""
        Helper function to determine whether the intersetion between two
        Spheres formed by growing menisci in the throats lies within the pore
        body that they connect to
        """
        intersection = np.zeros(len(r1), dtype=bool)
        
        vec_n = c2 - c1
        dist = np.linalg.norm(vec_n, axis=1)
        vec_a = self._perpendicular_vector(vec_n)
        dist_a = np.linalg.norm(vec_a, axis=1)
        vec_a *= 1/(np.vstack((dist_a, dist_a, dist_a)).T)
        #vec_a = self._norm_vector(vec_a)
        vec_b = np.cross(vec_a, vec_n)
        dist_b = np.linalg.norm(vec_b, axis=1)
        vec_b *= 1/(np.vstack((dist_b, dist_b, dist_b)).T)
        #vec_b = self._norm_vector(vec_b)
        x = (dist**2 - r2**2 + r1**2)/(2*dist)
        p = c1 + (vec_n.T*(x/dist)).T # intersection centre

        sq = 4 * dist**2 * r1**2 - (dist**2 - r2**2 + r1**2)**2
        
            #Could try to vectorize this but it would be pretty complicated!!!
        for i in range(len(r1)):
            if sq[i] > 0:
                # An intersection is found
                rad = (1/(2*dist[i]))*np.sqrt(sq[i])
                # Points around circle
                circle = self._circ_points(vec_a[i]*rad, vec_b[i]*rad, p[i])
                # If any points on the circle defining the intersection
                # Are within the pore radius distance from the pore center
                # The intersection is inside the pore - Not exact as pores
                # Are irregular shapes
                c2p = np.linalg.norm(circle-pore_center, axis=1)
                if np.any(c2p<pore_rad):
                    intersection[i]=True

        return intersection

    def _draw_sphere(self, centre, rad, c='r', fig=None, surface=False):
        if fig is None:
            fig = plt.figure()
        ax = fig.gca(projection='3d')
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        x = rad*np.cos(u)*np.sin(v)
        y = rad*np.sin(u)*np.sin(v)
        z = rad*np.cos(v)
    
        x += centre[0]
        y += centre[1]
        z += centre[2]
        if surface:
            ax.plot_surface(x, y, z, color=c)
        else:
            ax.plot_wireframe(x, y, z, color=c)
        return [x,y,z]
    
    def _setup_coop_filling(self, inv_points=None, filling_model=None, pores=None,
                            plot=False):
        r"""
        Evaluate the cooperative pore filling condition that the combined
        filling angle in next neighbor throats cannot exceed the geometric
        angle between their throat planes.
        This is used when the invading fluid has access to multiple throats 
        connected to a pore
        """
        from OpenPNM.Physics import models as pm
        if filling_model is None:
            logger.error('Please supply a model to calculate the filling angle')
        else:
            try:
                phys = self._phase.physics(self._phase.physics()[0])[0]
            except:
                logger.error('Phase has no Physics object associated')
        try:
            geom = self._net.geometries(self._net.geometries()[0])[0]
        except:
            logger.error('Model only works with Voronoi geometry... so far')
        if inv_points is None:
            inv_points = np.arange(0,30100, 500)
        if pores is None:
            pores = range(self.Np)

        self.tt_Pc = np.ndarray([self.Nt, self.Nt], dtype=float)
        self.tt_Pc.fill(sp.nan)
        # Dictionary with keys of capillary pressure 
        self.coop_data = {}
        start = time.time()
        for Pc in inv_points:
            # Dictionary with keys of pore id
            pore_data = {}
            phys.add_model(propname=filling_model,
                           model=pm.capillary_pressure.filling_angle_new,
                           r_toroid=geom._fibre_rad,
                           Pc=Pc)
            phys.add_model(propname='throat.meniscus_radius',
                           model=pm.capillary_pressure.meniscus_radius,
                           r_toroid=geom._fibre_rad,
                           filling_angle=filling_model)
            phys.add_model(propname='throat.meniscus_center',
                           model=pm.capillary_pressure.meniscus_center,
                           r_toroid=geom._fibre_rad,
                           filling_angle=filling_model)

            for pore in pores:
                # Dictionary with keys of throat id
                throat_data = {}
                p_cen = self._net['pore.centroid'][pore]
                p_rad = self._net['pore.indiameter'][pore]/2
                throats = self._net.find_neighbor_throats(pores=pore, flatten=True)
                throat_centres = self._net['throat.centroid'][throats]
                throat_normals = self._net['throat.normal'][throats]
                unit = np.linalg.norm(throat_normals, axis=1)
                throat_normals /= np.vstack((unit,unit,unit)).T
                v = p_cen - throat_centres
                sign = np.sign(np.sum(v*throat_normals, axis=1))
                cen = phys['throat.meniscus_center'][throats]
                c3 = np.vstack((cen*sign,cen*sign,cen*sign)).T
                men_cen = throat_centres + c3*throat_normals
                pairs = []
                for i,T in enumerate(throats):
                    men_data = {}
                    men_data['cen']=men_cen[i]
                    men_data['rad']=phys['throat.meniscus_radius'][T]
                    men_data['offset']=phys['throat.meniscus_center'][T]
                    men_data['alpha']=phys[filling_model][T]
                    throat_data[T]=men_data
                for ni in range(len(throats)):
                    for nj in range(len(throats))[ni+1:]:
                        pairs.append([ni,nj])
                pairs = np.asarray(pairs)
                if len(pairs) > 0:
                    c1 = men_cen[pairs[:,0]]
                    c2 = men_cen[pairs[:,1]]
                    c2c = (c1-c2)
                    dist = np.linalg.norm(c2c, axis=1)
                    t1 = throats[pairs[:,0]]
                    t2 = throats[pairs[:,1]]
                    r1 = phys['throat.meniscus_radius'][t1]
                    r2 = phys['throat.meniscus_radius'][t2]
                    check_pos = np.logical_and(r1 > 0, r2 > 0)
                    check_alpha_t1 = ~sp.isnan(phys[filling_model][t1])
                    check_alpha_t2 = ~sp.isnan(phys[filling_model][t2])
                    check_alpha = check_alpha_t1*check_alpha_t2
                    check_nans = sp.isnan(self.tt_Pc[t1, t2])
                    check_rads = (r1+r2) >= dist
                    mask = check_pos*check_alpha*check_nans*check_rads
                    if np.any(mask):
                        " Check if intersecting circle lies within pore "
                        inter = self._check_intersection(c1=c1[mask],
                                                         c2=c2[mask],
                                                         r1=r1[mask],
                                                         r2=r2[mask],
                                                         pore_center=p_cen,
                                                         pore_rad=p_rad)
                        if np.any(inter):
                            self.tt_Pc[t1[mask][inter], t2[mask][inter]] = Pc
                            self.tt_Pc[t2[mask][inter], t1[mask][inter]] = Pc
                            if plot:
                                fig = plt.figure()
                                ax = fig.gca(projection='3d')
                                inter_Ts = throats[np.unique(pairs[mask][inter])]
                                print(str(Pc) +' '+ str(inter_Ts))
                                print('alpha', phys[filling_model][inter_Ts])
                                print('rad', phys['throat.meniscus_radius'][inter_Ts])
                                print('cen', phys['throat.meniscus_center'][inter_Ts])
                                vo.plot_pore(geometry=geom,
                                             pores=[pore],
                                             fig=fig,
                                             include_points=False)
                                menisci = men_cen[np.unique(pairs[mask][inter])]
                                rads = phys['throat.meniscus_radius'][throats[np.unique(pairs[mask][inter])]]
                                for i, meniscus in enumerate(menisci):
                                    self._draw_sphere(meniscus, rads[i], c='b', fig=fig, surface=True)
                                ax.scatter(menisci[:, 0], menisci[:, 1], menisci[:, 2], c='b')
                        #inter
                    #pairs
                #pore
                pore_data[pore]=throat_data
            #Pc
            self.coop_data[Pc]=pore_data
        #Finished
        print("COOP FILL TIME", time.time()-start)
        ###
        #PLOT THE COOP PRESSURES
        ###
        
        fig = plt.figure()
        fig.subplots_adjust(hspace=0.4)
        fig.subplots_adjust(wspace=0.4)    
        ax1 = fig.add_subplot(111)
        coop_Pc = self.tt_Pc[~np.isnan(self.tt_Pc)]
        ax1.hist(coop_Pc, 60, facecolor='red', normed=True)
        ax1.set_xlabel('Coop Pore Filling Capillary Pressure')
        ax1.set_ylabel('Frequency')
        ax1.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

        #from scipy.stats import itemfreq
        #items = itemfreq(coop_Pc)
        #for item in items:
        #    print(item)

    def _check_coop(self, pore, queue):
        r"""
        Method run in loop after every pore invasion. All connecting throats
        are now given access to the invading phase. Two throats with access to
        the invading phase can cooperatively fill any pores that they are both
        connected to, common pores.
        The invasion of the throats connected to the common pore is handled
        elsewhere.
        """
        for throat in self._net.find_neighbor_throats(pores=pore):
            # A pore has just been invaded, all it's throats now have
            # An interface residing inside them
            if self['throat.inv_seq'][throat] == -1:
                # If the throat is not the invading throat that gave access
                # To this pore
                # Get the pores that this throat connects with
                a = set(self._net['throat.conns'][throat])
                # Get a list of pre-calculated coop filling pressures for all
                # Throats this throat can coop fill with
                ts_Pc = self.tt_Pc[throat]
                ts = np.argwhere(~np.isnan(ts_Pc))
                # If there are any potential coop filling throats
                if len(ts) > 0:
                    # For each throat find the common pore and the uncommon pores
                    for t in ts.flatten():
                        # Find common pore and uncommon pores
                        b = set(self._net['throat.conns'][t])
                        common_pore = list(a.intersection(b))
                        uncommon_pores = list(a.symmetric_difference(b))
                        # If the common pore is not invaded but the others are
                        # The potential coop filling event can now happen
                        # Add the coop pressure to the queue
                        if ((np.all(self['pore.inv_seq'][uncommon_pores] > -1)) and
                            (self['pore.inv_seq'][common_pore] == -1)):
                            # Coop pore filling fills the common pore
                            # The throats that gave access are not invaded now
                            # However, isolated throats between two invaded pores
                            # Are taken care of elsewhere...
                            # This could result in cyclycally calling the function
                            # If all elements are added to the queue and this 
                            # Function does not have access the the invasion sequence
                            hq.heappush(queue, [ts_Pc[t], list(common_pore), 'pore', 1])


    def _invade_isolated_Ts(self):
        r"""
        Throats that are uninvaded connected to pores that are both invaded
        should be invaded too.
        """
        Ts = self._net['throat.conns'].copy()
        invaded_Ps = self['pore.inv_seq']>-1
        uninvaded_Ts = self['throat.inv_seq']==-1
        isolated_Ts = np.logical_and(invaded_Ps[Ts[:,0]],
                                     invaded_Ps[Ts[:,1]])
        isolated_Ts = np.logical_and(isolated_Ts, uninvaded_Ts)
        if np.any(isolated_Ts):
            Pc = np.max(self['pore.inv_Pc'][Ts], axis=1)
            seq = np.max(self['pore.inv_seq'][Ts], axis=1)
            max_array = Ts[:,0]
            second_higher = self['pore.inv_Pc'][Ts][:,1] > self['pore.inv_Pc'][Ts][:,0]
            max_array[second_higher] = Ts[:,1][second_higher]
            action = self['pore.action'][max_array]
            self['throat.inv_Pc'][isolated_Ts] = Pc[isolated_Ts]
            self['throat.inv_seq'][isolated_Ts] = seq[isolated_Ts]
            self['throat.action'][isolated_Ts] = action[isolated_Ts]
        