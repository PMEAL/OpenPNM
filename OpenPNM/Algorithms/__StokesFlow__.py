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
        
    """
    
    def __init__(self,**kwargs):
        r"""
        Initializing the class
        """
        super(StokesFlow,self).__init__(**kwargs)
        self._logger.info("Create Hagen Poiseuille permeability Object")

            
    def _setup(self,
               hydraulic_conductance='hydraulic_conductance',
               x_term = 'pressure',
               occupancy='occupancy',
               **params):
        r"""
        This function executes the essential mathods before building matrices in Linear solution 
        """

        self._logger.info("Setup for Stokes Flow Algorithm")
        self._fluid = params['active_fluid']
        try: self._fluid = self._net._fluids[self._fluid] 
        except: pass #Accept object
        self._X_name = x_term
        success_1 = self._fluid.check_throat_health(props=occupancy)
        success_2 = self._fluid.check_throat_health(props=hydraulic_conductance)
        if not success_1:  
            self._fluid.set_data(prop=occupancy,throats='all',data=1)
            self._fluid.set_data(prop=occupancy,pores='all',data=1)
            self._logger.info('By default, it will be assumed that occupancy for '+self._fluid.name+' is equal to 1 in the entire network!')
        if success_2: 
            # Building hydraulic conductance based on occupancy
            g = self._fluid.get_throat_data(prop=hydraulic_conductance)
            s = self._fluid.get_throat_data(prop=occupancy)
            self._conductance = g*s+g*(-s)/1e3
            setup_conductance = True
        try:    setup_conductance
        except: raise Exception('There is an error for the throat property: '+ hydraulic_conductance+'!')
        try:    self.existing_bc
        except: raise Exception('There is an error in applying boundary conditions!')


    def _do_inner_iteration_stage(self):
        r"""
         main section of the algorithm              
        """
        p = self._do_one_inner_iteration()
        self.set_pore_data(prop=self._X_name,data = p)
        self._logger.info('Solving process finished successfully!')

    def update(self):
        
        p = self.get_pore_data(prop=self._X_name)
        self._fluid.set_pore_data(prop=self._X_name,data=p)
        self._logger.info('Results of ('+self.name+') algorithm have been updated successfully.')


    def effective_permeability(self,
                               fluid,
                               direction='',
                               d_term='viscosity',
                               x_term='pressure',
                               conductance='hydraulic_conductance',
                               occupancy='occupancy',
                               **params):
        r"""
        This function calculates effective diffusivity of a cubic network between face1 and face2.  
        face1 and face2 represent types of these two faces.

        """ 
        return self._calc_eff_prop(alg='Stokes',
                                  fluid=fluid,
                                  direction=direction,
                                  d_term=d_term,
                                  x_term=x_term,
                                  conductance=conductance,
                                  occupancy=occupancy)