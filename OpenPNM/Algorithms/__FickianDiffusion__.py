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
        self._logger.info("Create Fick's Diffusion Algorithm Object")


    def _setup(self,**params):
        r"""
        This function executes the essential methods specific to Fickian diffusion simulations
        """
        self._fluid = params['active_fluid']
        # Variable transformation for Fickian Algorithm from xA to ln(xB)
        Dir_pores = self._net.pore_properties['numbering'][self.BCtypes==1]
        self.BCvalues[Dir_pores] = sp.log(1-self.BCvalues[Dir_pores])
        g = self._fluid.throat_conditions['diffusive_conductance']
        s = self._fluid.throat_conditions['occupancy']
        self._conductance = g*s+g*(-s)/1e3

    def _do_inner_iteration_stage(self):

        X = self._do_one_inner_iteration()
        xA = 1-sp.exp(X)
        self._fluid.pore_conditions['mole_fraction'] = xA

    def calc_eff_diffusivity_cubic(self,face1=1,face2=6):
        r"""
        This function calculates effective diffusivity of a cubic network between face1 and face2.  
        face1 and face2 represent types of these two faces.

        """        
        face1_pores = self._net.pore_properties['numbering'][self._net.pore_properties['type']==face1]
        face2_pores = self._net.pore_properties['numbering'][self._net.pore_properties['type']==face2]
        xA = self._fluid.pore_conditions['mole_fraction']
        X1 = sp.log(1-xA[face1_pores])
        X2 = sp.log(1-xA[face2_pores])
        delta_X = sp.absolute(sp.average(X2)-sp.average(X1)) 
        C =sp.average(self._fluid.pore_conditions['molar_density'])
        
        coordx = self._net.pore_properties['coords'][:,0]
        coordy = self._net.pore_properties['coords'][:,1]
        coordz = self._net.pore_properties['coords'][:,2]
        if face1==1 or face1==6:
            L = sp.absolute(sp.unique(coordz[face1_pores])[0]-sp.unique(coordz[face2_pores])[0])*1e-6
            lx = (max(coordx[face1_pores]) - min(coordx[face1_pores]))*1e-6
            ly = (max(coordy[face1_pores]) - min(coordy[face1_pores]))*1e-6
            A = lx*ly            
        elif face1==2 or face1==5:
            L = sp.absolute(sp.unique(coordx[face1])[0]-sp.unique(coordx[face2])[0])*1e-6
            lz = (max(coordz[face1_pores]) - min(coordz[face1_pores]))*1e-6
            ly = (max(coordy[face1_pores]) - min(coordy[face1_pores]))*1e-6
            A = lz*ly  
        elif face1==3 or face1==4:
            L = sp.absolute(sp.unique(coordy[face1_pores])[0]-sp.unique(coordy[face2_pores])[0])*1e-6
            lx = (max(coordx[face1_pores]) - min(coordx[face1_pores]))*1e-6
            lz = (max(coordz[face1_pores]) - min(coordz[face1_pores]))*1e-6
            A = lx*lz
            
        fn = self._net.get_neighbor_pores(face1_pores)
        fn = fn[self._net.pore_properties['type'][fn]<1]
        ft = self._net.get_connecting_throat(face1_pores,fn)
        N = sp.sum(self._fluid.throat_conditions['diffusive_conductance'][ft]*sp.absolute(X1-sp.log(1-xA[fn])))
        Deff = N*L/(C*A*delta_X)
        return Deff