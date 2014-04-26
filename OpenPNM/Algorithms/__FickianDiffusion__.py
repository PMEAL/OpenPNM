#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: CEF PNM Team
# License: TBD
# Copyright (c) 2012

#from __future__ import print_function
"""

module __FickianDiffusion__: Fick's Law Diffusion
========================================================================

.. warning:: The classes of this module should be loaded through the 'Algorithms.__init__.py' file.

"""

import scipy as sp
from .__LinearSolver__ import LinearSolver

class FickianDiffusion(LinearSolver):
    r"""

    FickianDiffusion - Class to run Fick's law mass transfer diffusion on constructed networks

                        It returns conecentration gradient inside the network.

    """

    def __init__(self,**kwargs):
        r"""
        Initializing the class
        """
        super(FickianDiffusion,self).__init__(**kwargs)
        self._logger.info('Create Fickian Diffusion Algorithm Object')


    def _setup(self,
               conductance='diffusive_conductance',
               occupancy='occupancy',
               x_term='mole_fraction',               
               **params):
        r"""
        This function executes the essential methods specific to Fickian diffusion simulations
        """
        self._logger.info("Setup for Fickian Algorithm")        
        self._fluid = params['active_fluid']
        try: self._fluid = self._net._fluids[self._fluid]
        except: pass #Accept object
        self._X_name = x_term
        self._boundary_conditions_setup()
        # Variable transformation for Fickian Algorithm from xA to ln(xB)
        Dir_pores = self._net.get_pore_indices('all')[self._BCtypes==1]
        self._BCvalues[Dir_pores] = sp.log(1-self._BCvalues[Dir_pores])
        g = self._fluid.get_throat_data(prop=conductance)
        s = self._fluid.get_throat_data(prop=occupancy)
        self._conductance = g*s+g*(-s)/1e3
        

    def _do_inner_iteration_stage(self):

        X = self._do_one_inner_iteration()
        xA = 1-sp.exp(X)        
        self.set_pore_data(prop=self._X_name,data = xA)
        self._logger.info('Solving process finished successfully!')
              
    def update(self):
        
        x = self.get_pore_data(prop=self._X_name)        
        self._fluid.set_pore_data(prop=self._X_name,data=x)
        self._logger.info('Results of ('+self.name+') algorithm have been updated successfully.')
        

    def effective_diffusivity(self,
                              fluid,
                              direction='',                                                                               
                              conductance='diffusive_conductance',
                              occupancy='occupancy',
                              x_term='mole_fraction',
                              d_term='molar_density',
                              **params):
        r"""
        This function calculates effective diffusivity of a cubic network between face1 and face2.  
        face1 and face2 represent types of these two faces.

        """ 
        return self._calc_eff_prop(alg='Fickian',
                                  fluid=fluid,
                                  direction=direction,
                                  d_term=d_term,
                                  x_term=x_term,
                                  conductance=conductance,
                                  occupancy=occupancy)
