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

import scipy as sp
from .__LinearSolver__ import LinearSolver

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
        self._fluid = params['active_fluid']
        # Building thermal conductance
        g = self._fluid.throat_conditions['thermal_conductance']
        s = self._fluid.throat_conditions['occupancy']
        self._conductance = g*s+g*(-s)/1e3

    
    def _do_inner_iteration_stage(self):
        r"""
         main section of the algorithm              
        """
        T = self._do_one_inner_iteration()       
        self._fluid.pore_conditions['temperature'] = T
