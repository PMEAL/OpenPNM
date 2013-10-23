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

            
    def _setup(self,fluid_name):
        r"""
        This function executes the essential mathods before building matrices in Linear solution 
        """
        network = self._net
        # Building diffusive conductance       
        OpenPNM.Physics.FluidFlow.HagenPoiseuille(network,fluid_name)        
        OpenPNM.Physics.MultiPhase.full_pore_filling(network)
        OpenPNM.Physics.MultiPhase.conduit_filled_state_calculator(network)
        self._net.throat_conditions['eff_conductance'] = OpenPNM.Physics.MultiPhase.apply_phase_state_to_conduit_conductance(network,fluid_name)
         
   
    def _do_inner_iteration_stage(self):
        r"""
                       
        """
        p = self._do_one_inner_iteration()       
        self._net.pore_conditions['partial_pressure'] = p
        print p