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
from __TransportsLinearSolver__ import TransportsLinearSolver

class Permeability(LinearSolver):
    r"""   
       
    """
    
    def __init__(self,net,loglevel=10,**kwargs):
        r"""
        Initializing the class
        """
        super(Permeability,self).__init__(**kwargs)
        self._logger.info("Create Permeability Algorithm Object")

            
    def _setup_for_Solver(self,                                    
                                    loglevel=10,
                                    DryNetwork=0,
                                    T=353.15,
                                    P=1.01325e5,                                    
                                    **params):
        r"""
        Main features: Applying Boundary Conditions & Creating Transport Conductivities 

        This function executes the essential mathods for building matrices in Linear solution 
        """
        print '_setup_for hydraulic'
        self._net.pore_properties['temperature'] = T
        self._net.pore_properties['GasPressure'] = P
        self.DryNetwork = DryNetwork
        if 'Chydra' not in self._net.throat_properties:
            self._logger.info("Creating hydraulic conductivities for dry conduits")
            self._net.throat_properties['Chydra_dry'] = OpenPNM.Physics.FluidFlow.Conduits_Hydraulic_Conductivity(self._net,**params) 
        else:
            self.logger.info("Hydraulic conductivities for dry conduits have already been created")     

        if not self.DryNetwork:
            self.logger.info("Applying non-wetting phase saturation states to conduit conductivities")
            OpenPNM.Physics.MultiPhase.Conduit_Filled_State_Calculator(self._net,**params)
            self._net.throat_properties['Conductivities_Exp'] = OpenPNM.Physics.MultiPhase.Apply_Phase_State_to_Conduit_Conductivity(self._net,**params)
        else:
            self._net.throat_properties['Conductivities_Exp'] = self._net.throat_properties['Chyrda_dry']

    
    def _do_inner_iteration_stage(self,Experiment='Exp1'):
        r"""
                       
        """
        print '_do_outer_iteration_stage'
        self._net.pore_properties[Experiment+'_pressure'] = self._do_one_inner_iteration()