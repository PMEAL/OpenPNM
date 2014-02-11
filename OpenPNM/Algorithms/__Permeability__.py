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

class Permeability(LinearSolver):
    r"""   
    
    Hagen Poiseuille permeability - Class to run an algorithm to find permeability on constructed networks
    
                        It returns pressure gradient inside the network.                              
                            
    Parameter
    ----------
    -loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)    
                
        
    """
    
    def __init__(self,loglevel=10,**kwargs):
        r"""
        Initializing the class
        """
        super(Permeability,self).__init__(**kwargs)
        self._logger.info("Create Hagen Poiseuille permeability Object")

            
    def _setup(self,**params):
        r"""
        This function executes the essential mathods before building matrices in Linear solution 
        """
        self._fluid = params['active_fluid']
        # Building hydraulic conductance
        g = self._fluid.throat_conditions['hydraulic_conductance']
        s = self._fluid.throat_conditions['occupancy']
        self._conductance = g*s+g*(-s)/1e3
   
    def _do_inner_iteration_stage(self):
        r"""
        The main section of the algorithm               
        """
        p = self._do_one_inner_iteration()       
        self._fluid.pore_conditions['pressure'] = p

    def calc_eff_permeability_cubic(self,face1=1,face2=6):
        r"""
        This function calculates effective permeability of a cubic network between face1 and face2.  
        face1 and face2 represent types of these two faces.

        """        
        face1_pores = self._net.pore_properties['numbering'][self._net.pore_properties['type']==face1]
        face2_pores = self._net.pore_properties['numbering'][self._net.pore_properties['type']==face2]
        p = self._fluid.pore_conditions['pressure']
        P1 = p[face1_pores]
        P2 = p[face2_pores]
        delta_P = sp.absolute(sp.average(P2)-sp.average(P1)) 
        mu =sp.average(self._fluid.pore_conditions['viscosity'])
        
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
            
        fn = self._net.find_neighbor_pores(face1_pores)
        fn = fn[self._net.pore_properties['type'][fn]<1]
        ft = self._net.find_connecting_throat(face1_pores,fn)
        Q = sp.sum(self._fluid.throat_conditions['hydraulic_conductance'][ft]*sp.absolute(P1-p[fn]))
        Keff = Q*mu*L/(A*delta_P)
        return Keff