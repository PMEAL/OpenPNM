#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: CEF PNM Team
# License: TBD
# Copyright (c) 2012

#from __future__ import print_function
"""
module __LinearSolver__: Algorithm for solving transport processes
========================================================================

.. warning:: The classes of this module should be loaded through the 'Algorithms.__init__.py' file.

"""

import OpenPNM
import scipy as sp
import scipy.sparse as sprs
import scipy.sparse.linalg as splin
from __GenericAlgorithm__ import GenericAlgorithm
import matplotlib.pylab as plt 


class LinearSolver(GenericAlgorithm):
    r"""
    This class provides essential methods for building and solving matrices in a transport process.
    It will be inherited during Fickian diffusion, Fourier heat conduction, Hagen Poiseuille permeability and electron conduction.
        
    """
    
    def __init__(self,loglevel=10,**kwargs):
        r"""
        Initializing the class
        """
        super(LinearSolver,self).__init__(**kwargs)
        self._logger.info("Create Linear Solver Algorithm Object")
            

    def _do_one_inner_iteration(self):
        
        A = self._build_coefficient_matrix() 
        B = self._build_RHS_matrix()
        X = splin.spsolve(A,B) 
        return(X)
    
    def set_boundary_conditions(self,types=[],values=[]):
        r"""        
        Assigning Type and Value for Boundary Pores.
        
        - Types of Boundary Values:
           internal = 0, Dirichlet = 1, Nuemann_flux = 2, Nuemann_insulated = 3, Nuemann_rate = 4

        - For any network, it is possible to apply BC to arbitrary pores, by defining 
          two arrays: types and values. For example:
          
          types = array([2,4,1,0,...,1])
          values = array([0.1,0.0087,0.5,0,...,0.30])
          
          It means that:
              for pore 0: Nuemann_boundary flux = 0.1 
              for pore 1: Nuemann_boundary rate = 0.0087
              for pore 2: Dirichlet_boundary value = 0.5
              for pore 3: Internal pore without imposed boundary condition
              .
              .
              .
              for pore Np : Dirichlet_boundary value = 0.30
              
        Notes
        -----
        Nuemann_isolated is equivalent to Nuemann_flux boundary condition with flux = 0. 
        Negative value for Nuemann_rate or Nuemann_flux means that the quanitity of interest leaves the network. 

        """
        setattr(self,"BCtypes",types)
        setattr(self,"BCvalues",values)      
           
    def _build_coefficient_matrix(self):
        
        # Filling coefficient matrix

        tpore1 = self._net.throat_properties['connections'][:,0]
        tpore2 = self._net.throat_properties['connections'][:,1]
        row = tpore1
        col = tpore2
        data= self._net.throat_conditions['eff_conductance']
        
        pnum = self._net.pore_properties['numbering']      
        loc1 = sp.in1d(tpore2,tpore2[sp.in1d(tpore2,pnum[self.BCtypes!=1])])
        modified_tpore2 = tpore2[loc1]
        modified_tpore1 = tpore1[loc1]        
        row = sp.append(row,modified_tpore2)                
        col = sp.append(col,modified_tpore1)
        data = sp.append(data,self._net.throat_conditions['eff_conductance'][loc1])
        
        A_dim = self._net.get_num_pores()

        if (self.BCtypes==4).any():
            self.extera_Nuemann_equations = sp.unique(self.BCvalues[self.BCtypes==4])
            A_dim = A_dim + len(self.extera_Nuemann_equations)
            extera_neu = self.extera_Nuemann_equations
            g_super = 1e2
            for item in range(len(extera_neu)):
                loc_neu = sp.in1d(tpore2,pnum[self.BCvalues==extera_neu[item]])
                neu_tpore2 = tpore2[loc_neu]
                
                row = sp.append(row,neu_tpore2) 
                col = sp.append(col,len(neu_tpore2)*[A_dim-item-1])
                data = sp.append(data,len(neu_tpore2)*[-g_super])

                row = sp.append(row,len(neu_tpore2)*[A_dim-item-1]) 
                col = sp.append(col,neu_tpore2)
                data = sp.append(data,len(neu_tpore2)*[-g_super])
            
        else:
            self.extera_Nuemann_equations = 0
        
        self.Coeff_dimension = A_dim
        
        # Adding positions for diagonal
        
        row = sp.append(row,range(0,A_dim))
        col = sp.append(col,range(0,A_dim))
        data = sp.append(data,sp.zeros((A_dim,1)))

        a = sprs.coo.coo_matrix((data,(row,col)),(A_dim,A_dim))
        A = a.tocsr()        
       
        for i in range(0,A_dim):                    
                           
            if (i<self._net.get_num_pores() and self.BCtypes[i]==1):
                A[i,i] = 1 
            else:
                A[i,i] = -sp.sum(A[i,:][sp.nonzero(A[i,:])])

        return(A)
        
        
    def _build_RHS_matrix(self):
        
        extera_neu = self.extera_Nuemann_equations  
        A_dim = self.Coeff_dimension
        B = sp.zeros([A_dim,1])
        Dir_pores = self._net.pore_properties['numbering'][self.BCtypes==1]
        B[Dir_pores] = sp.reshape(self.BCvalues[Dir_pores],[len(Dir_pores),1])
        if (self.BCtypes==4).any():
            for item in range(len(extera_neu)):                
                B[A_dim-item-1,0] = extera_neu[item]
                
        return(B)