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

    Parameter
    ----------
    -loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)


    """

    def __init__(self,loglevel=10,**kwargs):
        r"""
        Initializing the class
        """
        super(FickianDiffusion,self).__init__(**kwargs)
        self._logger.info('Create Fickian Diffusion Algorithm Object')


    def _setup(self,
               diffusive_conductance='diffusive_conductance',
               occupancy='occupancy',
               mole_fraction='mole_fraction',
               **params):
        r"""
        This function executes the essential methods specific to Fickian diffusion simulations
        """
        self._logger.info("Setup for Fickian Algorithm")        
        self._fluid = params['active_fluid']
        self._X_name = mole_fraction
        self._boundary_conditions_setup()
        # Variable transformation for Fickian Algorithm from xA to ln(xB)
        Dir_pores = self._net.get_pore_indices('all')[self._BCtypes==1]
        self._BCvalues[Dir_pores] = sp.log(1-self._BCvalues[Dir_pores])
        g = self._fluid.get_throat_data(prop=diffusive_conductance)
        s = self._fluid.get_throat_data(prop=occupancy)
        self._conductance = g*s+g*(-s)/1e3

    def _do_inner_iteration_stage(self):

        X = self._do_one_inner_iteration()
        xA = 1-sp.exp(X)        
        self.set_pore_data(prop=self._X_name,data = xA)
        print('Solving process finished successfully!')
              
    def update(self):
        
        x = self.get_pore_data(prop=self._X_name)        
        self._net.set_pore_data(phase=self._fluid,prop=self._X_name,data=x)
        print('Results of this algorithm have been updated successfully.')
        

    def effective_diffusivity_cubic(self,
                                   fluid,
                                   face1='',
                                   face2='',
                                   d_term='molar_density',
                                   x_term='mole_fraction',
                                   conductance='diffusive_conductance',
                                   occupancy='',
                                   **params):
        r"""
        This function calculates effective diffusivity of a cubic network between face1 and face2.  
        face1 and face2 represent types of these two faces.

        """ 
        return self._calc_eff_prop_cubic(alg='Fickian',
                                  fluid=fluid,
                                  face1=face1,
                                  face2=face2,
                                  d_term=d_term,
                                  x_term=x_term,
                                  conductance=conductance,
                                  occupancy=occupancy)
