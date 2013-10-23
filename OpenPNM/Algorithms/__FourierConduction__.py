#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: CEF PNM Team
# License: TBD
# Copyright (c) 2012

#from __future__ import print_function
"""

module __FourierConduction__: 
========================================================================

.. warning:: The classes of this module should be loaded through the 'Algorithms.__init__.py' file.

"""

import OpenPNM
import scipy as sp
from __LinearSolver__ import LinearSolver

class FourierConduction(LinearSolver):
    r"""   
    
    Fourier Conduction - Class to run an algorithm for heat conduction through flow on constructed networks
    
                        It returns temperature gradient inside the network.
                                  
                            
    Parameter
    ----------
    -loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)    
                
        
    """
    
    def __init__(self,loglevel=10,**kwargs):
        r"""
        Initializing the class
        """
        super(FourierConduction,self).__init__(**kwargs)
        self._logger.info("Create Fourier Conduction Algorithm Object")
            
    def _setup(self,**params):
        r"""

        This function executes the essential mathods for building matrices in Linear solution 
        """
        network = self._net
        self.fluid_name = params['fluid_name']
        network.refresh_fluid(self.fluid_name)
        # Building hydraulic conductance
        OpenPNM.Physics.HeatConduction.ThermalConductance(network,self.fluid_name)
#        method = params['conduit_filling_method']
#        OpenPNM.Physics.MultiPhase.full_pore_filling(network)
#        OpenPNM.Physics.MultiPhase.calc_conduit_filling(network,method)
        g = network.throat_conditions['thermal_conductance'+'_'+self.fluid_name]
#        c = pn.throat_conditions['']
        self._conductance = g

    
    def _do_inner_iteration_stage(self):
        r"""
                       
        """
        T = self._do_one_inner_iteration()       
        self._net.pore_conditions['temperature'+'_'+self.fluid_name] = T
        print T