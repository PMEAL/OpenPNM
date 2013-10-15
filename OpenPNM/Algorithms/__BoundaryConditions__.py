# -*- coding: utf-8 -*-
"""
Created on Tue Jul 09 10:54:37 2013

@author: Mahmoudreza Aghighi

"""
import OpenPNM
import numpy as np
from __GenericAlgorithm__ import GenericAlgorithm

class BoundaryConditions(GenericAlgorithm):
        def __init__(self,net=OpenPNM.Network.GenericNetwork(),**kwords):
            r"""
            Initialize
            """
            super(BoundaryConditions,self).__init__(**kwords)
            self.indent = ""
            self._logger.debug("Construct class")
            self._net = net
        def _Boundary_Pores_Conditions(self,FaceTypes=[4,1,4,4,1,4],FaceValues=[0,0.2,0,0,0.5,0],BList=[]):
            r"""        
            - Assigning Type and Value for Boundary Pores.
            - Types of Boundary Values in FaceTypes:
               internal=0, Dirichlet=1, flux=2, periodic=3, insulated=4 mass_rate=5
            - In FaceValues, each item represents the value of that face according
              to the boundary condition type. Except for mass_rate condition which 
              the value represents the total mass rate applied to the entire face,
              not just a single pore. 
              For instance: if we have Facetypes[2]=1 and FaceValues[2]=0.75,it means
              that for face 3, concentration=0.75
              but if we have FaceTypes[4]=2 and FaceValues[4]=0.8, it means that
              for face 5, flux=0.8
            - Instead of applying boundary conditions to entire face, by using 
              BLists, it is possible to apply BC to arbitrary pores. For instance:
              BLists = [(459,2,0.097),([2043,36,443,2],1,0.5)]
              It means that for Pore 459: flux=0.097 & for Pores [2043,36,443,2]: Concentration=0.5
            """
            self._net.pore_properties['BCtype'] = np.zeros(self._net.get_num_pores())
            self._net.pore_properties['BCvalue'] = np.zeros(self._net.get_num_pores())
            for i in range(1,7):
                self._net.pore_properties['BCtype'][self._net.pore_properties['type']==i] = FaceTypes[i-1]
                self._net.pore_properties['BCvalue'][self._net.pore_properties['type']==i] = FaceValues[i-1]
            if BList:
                for item in BList:
                    pores = item[0]
                    self._net.pore_properties['BCtype'][pores] = item[1]
                    self._net.pore_properties['BCvalue'][pores] = item[2]