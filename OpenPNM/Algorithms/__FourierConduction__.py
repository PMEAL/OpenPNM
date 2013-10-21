# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 10:50:28 2013

@author: Mahmoudreza Aghighi

module __FourierConduction__: 
========================================================================

.. warning:: The classes of this module should be loaded through the 'Algorithms.__init__.py' file.

"""

import OpenPNM
import scipy as sp
from __LinearSolver__ import LinearSolver

class FourierConduction(LinearSolver):
    r"""   
            
    """
    
    def __init__(self,loglevel=10,**kwargs):
        r"""
        Initializing the class
        """
        super(FourierConduction,self).__init__(**kwargs)
        self._logger.info("Create Fourier Conduction Algorithm Object")
        print 'init'
            
    def _setup_for_Solver(self,loglevel=10,
                               T=353.15,                                    
                               **params):
        r"""

        This function executes the essential mathods for building matrices in Linear solution 
        """
        print '_setup_'
        self._net.pore_properties['temperature'] = T
        if 'Cdiff' not in self._net.throat_properties:
            self._logger.info("Creating Fourier thermal conduction conductivities for for solid network")
            self._net.throat_properties['TCFiber'] = OpenPNM.Physics.ElectricTransport.Conduits_Conductivity(self._net,**params) 
        else:
            self.logger.info("Fourier thermal conductivities for solid network have already been created")     

        self._net.throat_properties['Conductivities_Exp'] = self._net.throat_properties['TCFiber']

    
    def _do_inner_iteration_stage(self,Experiment='Exp1'):
        r"""
                       
        """
        print '_do_outer_iteration_stage'
        self._net.pore_properties[Experiment+'_solid_temperature'] = self._do_one_inner_iteration()