#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: CEF PNM Team
# License: TBD
# Copyright (c) 2012

#from __future__ import print_function
"""

module __OhmicConduction__:
========================================================================

"""

import scipy as sp
from .__LinearSolver__ import LinearSolver

class OhmicConduction(LinearSolver):
    r"""

    OhmicConduction - Class to run an algorithm for electron conduction on constructed networks

                        It returns voltage gradient inside the network.

    """

    def __init__(self,**kwargs):
        r"""
        Initializing the class
        """
        super(OhmicConduction,self).__init__(**kwargs)
        self._logger.info("Create Ohmic Conduction Algorithm Object")

    def _setup(self,
               loglevel=10,
               electronic_conductance='electronic_conductance',
               occupancy='occupancy',
               voltage='voltage',
               **params):
        r"""

        This function executes the essential mathods for building matrices for Linear solution
        """
        self._fluid = params['active_fluid']
        try: self._fluid = self.find_object_by_name(self._fluid) 
        except: pass #Accept object
        self._X_name = voltage
        self._boundary_conditions_setup()
        g = self._fluid.get_throat_data(prop=electronic_conductance)
        s = self._fluid.get_throat_data(prop=occupancy)
        self._conductance = g*s+g*(-s)/1e3


    def _do_inner_iteration_stage(self):
        v = self._do_one_inner_iteration()
        self.set_pore_data(prop=self._X_name,data= v)
        self._logger.info('Solving process finished successfully!')
    
    def update(self):
        
        v = self.get_pore_data(prop=self._X_name)
        self._net.set_pore_data(phase=self._fluid,prop=self._X_name,data=v)
        self._logger.info('Results of ('+self.name+') algorithm have been updated successfully.')
        