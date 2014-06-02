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
        
    """
    
    def __init__(self,**kwargs):
        r"""
        Initializing the class
        """
        super(FourierConduction,self).__init__(**kwargs)
        self._logger.info("Create Fourier Conduction Algorithm Object")
            
    def _setup(self,
               thermal_conductance='thermal_conductance',
               occupancy='occupancy',
               **params):
        r"""

        This function executes the essential mathods for building matrices in Linear solution 
        """   

        self._logger.info("Setup for Fourier Algorithm")        
        self._fluid = params['active_fluid']
        try: self._fluid = self._net._fluids[self._fluid] 
        except: pass #Accept object
        success_1 = self._fluid.check_throat_health(props=occupancy)
        success_2 = self._fluid.check_throat_health(props=thermal_conductance)
        if not success_1:  
            self._fluid.set_data(prop=occupancy,throats='all',data=1)
            self._fluid.set_data(prop=occupancy,pores='all',data=1)
            self._logger.info('By default, it will be assumed that occupancy for '+self._fluid.name+' is equal to 1 in the entire network!')
        if success_2: 
            # Building thermal conductance based on occupancy
            g = self._fluid.get_throat_data(prop=thermal_conductance)
            s = self._fluid.get_throat_data(prop=occupancy)
            self._conductance = g*s+g*(-s)/1e3
            setup_conductance = True
        try:    setup_conductance
        except: raise Exception('There is an error for the throat property: '+thermal_conductance+'!')
        try:    self.existing_bc
        except: raise Exception('There is an error in applying boundary conditions!')

    
    def _do_inner_iteration_stage(self):
        r"""
         main section of the algorithm              
        """
        T = self._do_one_inner_iteration()       
        self.set_pore_data(prop='temperature',data= T)
        self._logger.info('Solving process finished successfully!')


    def update(self):
        
        T = self.get_pore_data(prop='temperature')
        self._fluid.set_pore_data(prop='temperature',data=T)
        self._logger.info('Results of ('+self.name+') algorithm have been updated successfully.')