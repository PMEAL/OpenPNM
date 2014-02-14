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

import scipy as sp
from .__LinearSolver__ import LinearSolver

class StokesFlow(LinearSolver):
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
        super(StokesFlow,self).__init__(**kwargs)
        self._logger.info("Create Hagen Poiseuille permeability Object")

            
    def _setup(self,**params):
        r"""
        This function executes the essential mathods before building matrices in Linear solution 
        """
      
        self._fluid = params['active_fluid']
        # Building hydraulic conductance
        g = self._fluid.get_throat_data(prop='hydraulic_conductance')
        s = self._fluid.get_throat_data(prop='occupancy')
        self._conductance = g*s+g*(-s)/1e3

    def _do_inner_iteration_stage(self):

        p = self._do_one_inner_iteration()
        self._fluid.set_pore_data(prop='pressure',data = p)

