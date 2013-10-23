#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: CEF PNM Team
# License: TBD
# Copyright (c) 2012

#from __future__ import print_function
"""

module __Permeability__: Hagen Poiseuille permeability
========================================================================

.. warning:: The classes of this module should be loaded through the 'Algorithms.__init__.py' file.

"""

import OpenPNM
import scipy as sp
from __LinearSolver__ import LinearSolver

class Permeability(LinearSolver):
    r"""   
    
    Hagen Poiseuille permeability - Class to run an algorithm to find permeability on constructed networks
    
                        It returns pressure gradient inside the network.
                                  
                            
    Parameter
    ----------
    -loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)    
                
        
    """
    
    def __init__(self,loglevel=10,**kwargs):
        r"""
        Initializing the class
        """
        super(Permeability,self).__init__(**kwargs)
        self._logger.info("Create Hagen Poiseuille permeability Object")

            
    def _setup(self,**params):
        r"""
        This function executes the essential mathods before building matrices in Linear solution 
        """
        network = self._net
        self.fluid_name = params['fluid_name']
        network.refresh_fluid(self.fluid_name)
        # Building hydraulic conductance
        OpenPNM.Physics.FluidFlow.HydraulicConductance(network,self.fluid_name)
#        method = params['conduit_filling_method']
#        OpenPNM.Physics.MultiPhase.full_pore_filling(network)
#        OpenPNM.Physics.MultiPhase.calc_conduit_filling(network,method)
        g = network.throat_conditions['hydraulic_conductance'+'_'+self.fluid_name]
#        c = pn.throat_conditions['']                
        self._conductance = g
   
    def _do_inner_iteration_stage(self):
        r"""
                       
        """
        p = self._do_one_inner_iteration()       
        self._net.pore_conditions['partial_pressure'+'_'+self.fluid_name] = p
        print p
        
    def calc_total_permeability(self):
        
        network = self._net        
        pores1 = 0
        pores2 = 0
        area = 0
        length = 0 
                
        total_rate = 0 
        dp = network.pore_condtions['partial_pressure'][pores1] - network.pore_condtions['partial_pressure'][pores2]
        K = sp.absolute(total_rate*length/area*(sp.average(network.pore_conditions['viscosity'+'_'+self.fluid_name]/dp)))
        self._network_permeability = K