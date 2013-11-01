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
        return(X[range(self._net.get_num_pores())])

    def set_boundary_conditions(self,types=[],values=[]):
        r"""
        Assigning Type and Value for Boundary Pores.

        - Types of Boundary Conditions:
           internal = 0, Dirichlet = 1, Neumann_flux = 2, Neumann_insulated = 3, Neumann_rate = 4

        - For any network, it is possible to apply BC to arbitrary pores, by defining
          two arrays: types and values. For example:

            BCtypes = array([2, 1, 4 ,0, 4, 3, 4, 1]) 
            BCvalues = array([0.1, 0.5, 0.87, 0, -0.35, 0, 0.87, 0.30])
      
            It means that:
                  
                  for pore 0: Neumann, flux = 0.1
                  
                  for pore 1: Dirichlet, value = 0.5

                  for pore 2: Neumann, rate = 0.0087 
                  (hint: Since there are two pores (2,6) with Neumann_rate type which 
                  have the exact same amount of rate, algorithm assumes that 0.0087 is 
                  the rate of quantity of interest which leaves both pore 2 and 6)
                          
                  for pore 3: Internal pore without imposed boundary condition
                 (hint: If pore 3 is a boundary pore (a pore in boundary faces), algorithm 
                 by default assumes that, this is Neumann_insulated pore.)
            
                  for pore 4: Neumann, rate= -0.35
                  (hint: There is only one pore with Neumann_rate type and value of -0.35. So 
                  algorithm assumes that 0.35 is the rate of quantity of interest which is only entering pore 4)
                        
                  for pore 5: Neumann_insulated, value=0

                  for pore 6 : Neumann, rate=0.0087
                        (hint: Refer to pore 2)              .
                          
                  for pore 7 : Dirichlet, value = 0.30


        Notes
        -----
        - Neumann_insulated is equivalent to Neumann_flux boundary condition when flux
        is zero. Therefore, there is no need to define BCvalues for this kind of boundary condition.
        - In Fickian algorithm, positive value for Neumann_rate or Neumann_flux for a boundary pore means 
        that the quantity of interest leaves the pore, but for any other algorithms, positive Neumann value
        means that the quantity of interest enters this pore. 

        """
        setattr(self,"BCtypes",types)
        setattr(self,"BCvalues",values)
        

    def _build_coefficient_matrix(self):

        # Filling coefficient matrix

        pnum = self._net.pore_properties['numbering']        
        tpore1 = self._net.throat_properties['connections'][:,0]
        tpore2 = self._net.throat_properties['connections'][:,1]
        
        loc1 = sp.in1d(tpore1,tpore1[sp.in1d(tpore1,pnum[self.BCtypes!=1])])
        modified_tpore1 = tpore1[loc1]
        modified_tpore2 = tpore2[loc1]
        row = modified_tpore1
        col = modified_tpore2
        data_main = self._conductance
        data_main = data_main + 1e-30 #to prevent matrix singularities        
        data = data_main[loc1]
        
        loc2 = sp.in1d(tpore2,tpore2[sp.in1d(tpore2,pnum[self.BCtypes!=1])])
        modified_tpore2 = tpore2[loc2]
        modified_tpore1 = tpore1[loc2]
        row = sp.append(row,modified_tpore2)
        col = sp.append(col,modified_tpore1)
        data = sp.append(data,data_main[loc2])

        A_dim = self._net.get_num_pores()
        
        if (self.BCtypes==2).any():
            flux_pores = self._net.pore_properties['numbering'][self.BCtypes==2]
            flux_values = sp.unique(self.BCvalues[self.BCtypes==2])
            for i in range(len(flux_values)):
                f = flux_pores[sp.in1d(flux_pores,self._net.pore_properties['numbering'][self.BCvalues==flux_values[i]])]
                fn = self._net.get_neighbor_pores(f)
                fn = fn[self._net.pore_properties['type'][fn]<1]
                ft = self._net.get_connecting_throat(f,fn)
                self.BCtypes[f] = 4
                self.BCvalues[f] = sp.sum(self.BCvalues[f]*(self._net.throat_properties['diameter'][ft])**2)
        
        boundaries = self._net.pore_properties['numbering'][self._net.pore_properties['type']>0]
        if (boundaries[self.BCtypes[boundaries]==0]).any():
            self.BCtypes[boundaries[self.BCtypes[boundaries]==0]] = 3
        
        if (self.BCtypes==4).any():
            self.extera_Neumann_equations = sp.unique(self.BCvalues[self.BCtypes==4])
            A_dim = A_dim + len(self.extera_Neumann_equations)
            extera_neu = self.extera_Neumann_equations
            g_super = 1e2
            for item in range(len(extera_neu)):
                loc_neu = sp.in1d(tpore2,pnum[self.BCvalues==extera_neu[item]])
                neu_tpore2 = tpore2[loc_neu]

                row = sp.append(row,neu_tpore2)
                col = sp.append(col,len(neu_tpore2)*[A_dim-item-1])
                data = sp.append(data,len(neu_tpore2)*[g_super])

                row = sp.append(row,len(neu_tpore2)*[A_dim-item-1])
                col = sp.append(col,neu_tpore2)
                data = sp.append(data,len(neu_tpore2)*[g_super])

        else:
            self.extera_Neumann_equations = 0

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

        extera_neu = self.extera_Neumann_equations
        A_dim = self.Coeff_dimension
        B = sp.zeros([A_dim,1])
        Dir_pores = self._net.pore_properties['numbering'][self.BCtypes==1]
        B[Dir_pores] = sp.reshape(self.BCvalues[Dir_pores],[len(Dir_pores),1])
        if (self.BCtypes==4).any():
            for item in range(len(extera_neu)):
                B[A_dim-item-1,0] = extera_neu[item]
        
        return(B)