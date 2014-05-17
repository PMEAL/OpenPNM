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
               occupancy='occupancy',
               **params):
        r"""
        This function executes the essential mathods before building matrices in Linear solution 
        """
        self._setup = 0 
        try :
            self.bc_setup
            if self.bc_setup==1:
                self._logger.info("Setup for Stokes Flow Algorithm")
                self._fluid = params['active_fluid']
                try: self._fluid = self._net._fluids[self._fluid] 
                except: pass #Accept object
                self._X_name = 'pressure'
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
                    self._setup = 1
                else: 
                    self._logger.error('In '+self._fluid.name+', there is an error for the property: '+hydraulic_conductance)
                    self._setup = 0
            else: 
                self._logger.error('There is an error in applying boundary conditions!')
                self._setup = 0
        except:
            self._logger.error('No boundary condition has been implemented for algorithm: '+self.name+'!') 


    def _do_inner_iteration_stage(self):
        
        if self._setup==1:
            p = self._do_one_inner_iteration()
            self.set_pore_data(prop='pressure',data = p)
            self._logger.info('Solving process finished successfully!')
        else: raise Exception('Error in setup section of the algorithm!'+self.name+' cannot be executed!')

    def update(self):
        
        p = self.get_pore_data(prop='pressure')
        self._fluid.set_pore_data(prop='pressure',data=p)
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