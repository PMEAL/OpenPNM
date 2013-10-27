#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: CEF PNM Team
# License: TBD
# Copyright (c) 2012

#from __future__ import print_function
"""

module __ElectronConduction__: 
========================================================================

"""

import OpenPNM
import scipy as sp
from __LinearSolver__ import LinearSolver

class ElectronConduction(LinearSolver):
    r"""

    ElectronConduction - Class to run an algorithm for electron conduction on constructed networks

                        It returns ------- gradient inside the network.


    Parameter
    ----------
    -loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)


    """
    
    def __init__(self,loglevel=10,**kwargs):
        r"""
        Initializing the class
        """
        super(ElectronConduction,self).__init__(**kwargs)
        self._logger.info("Create Electrons Conduction Algorithm Object")
            
    def _setup(self,loglevel=10,**params):
        r"""

        This function executes the essential mathods for building matrices for Linear solution 
        """
        self._fluid = params['active_fluid']
        # Building electron conductance
        OpenPNM.Physics.ElectronConduction.ElectronicConductance(self._net,self._fluid)
        g = self._fluid.throat_conditions['electron_conductance']
        s = self._fluid.throat_conditions['occupancy']
        self._conductance = g*s

    
    def _do_inner_iteration_stage(self):
        v = self._do_one_inner_iteration()       
        self._fluid.pore_conditions['voltage'] = v
