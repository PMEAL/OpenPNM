# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 10:50:28 2013

@author: Mahmoudreza Aghighi

module __FickianDiffusion__: Fick's Law Diffusion
========================================================================

.. warning:: The classes of this module should be loaded through the 'Algorithms.__init__.py' file.

"""

import OpenPNM
import scipy as sp
from __LinearSolver__ import LinearSolver

class FickianDiffusion(LinearSolver):
    r"""   
    
    FickianDiffusion - Class to run Fick's law mass transfer diffusion on constructed networks
    
                        It returns conecentration gradient inside the network.
                        An invasion algorithm should be used before running diffusion calculations.
                                   
                            
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
        self._logger.info("Create Fick's Diffusion Algorithm Object")

            
    def _setup(self,fluid_name):
        r"""
        This function executes the essential mathods before building matrices in Linear solution 
        """
        network = self._net        
        Dir_pores = network.pore_properties['numbering'][self.BCtypes==1]
        self.BCvalues[Dir_pores] = sp.log(1-self.BCvalues[Dir_pores])        
        OpenPNM.Physics.MassTransport.DiffusiveConductance(network,fluid_name)        
        OpenPNM.Physics.MultiPhase.full_pore_filling(network,Pc=0.0,Seq=0)
        OpenPNM.Physics.MultiPhase.conduit_filled_state_calculator(network)
        self._net.throat_conditions['eff_conductance'] = OpenPNM.Physics.MultiPhase.apply_phase_state_to_conduit_conductance(network,fluid_name)
         
   
    def _do_inner_iteration_stage(self):
        r"""
                       
        """
        print '_do_outer_iteration_stage'
        X = self._do_one_inner_iteration()
        xA = 1-sp.exp(X)        
        self._net.pore_conditions['mole_fractions'] = xA

