# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 10:50:28 2013

@author: Mahmoudreza Aghighi

module __ElectronConduction__: 
========================================================================

.. warning:: The classes of this module should be loaded through the 'Algorithms.__init__.py' file.

"""

import OpenPNM
import scipy as sp
from __TransportsLinearSolver__ import TransportsLinearSolver

class ElectronConduction(TransportsLinearSolver):
    r"""   
            
    """
    
    def __init__(self,net,loglevel=10,**kwargs):
        r"""
        Initializing the class
        """
        super(ElectronConduction,self).__init__(net = net,**kwargs)
        self._logger.info("Create Electrons Conduction Algorithm Object")
        print 'init'
            
    def _setup_for_TransportSolver(self,                                    
                                    loglevel=10,
                                    T=353.15,                                    
                                    **params):
        r"""

        This function executes the essential mathods for building matrices in Linear solution 
        """
        print '_setup_'
        self._net.pore_properties['temperature'] = T
        if 'Cdiff' not in self._net.throat_properties:
            self._logger.info("Creating electric conductivities for for solid network")
            self._net.throat_properties['ECFiber'] = OpenPNM.Physics.ElectricTransport.Conduits_Conductivity(self._net,**params) 
        else:
            self.logger.info("Electric conductivities for solid network have already been created")     

        self._net.throat_properties['Conductivities_Exp'] = self._net.throat_properties['ECFiber']

    
    def _do_inner_iteration_stage(self,Experiment='Exp1'):
        r"""
                       
        """
        print '_do_outer_iteration_stage'
        self._net.pore_properties[Experiment+'_solid_voltage'] = self._do_one_inner_iteration()