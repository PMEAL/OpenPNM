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
import scipy.sparse as sprs
import scipy.sparse.linalg as splin
from __GenericAlgorithm__ import GenericAlgorithm

class LinearSolver(GenericAlgorithm):
    r"""   

        
    """
    
    def __init__(self,loglevel=10,**kwargs):
        r"""
        Initializing the class
        """
        super(LinearSolver,self).__init__(**kwargs)
        self._logger.info("Create Linear Solver Algorithm Object")
            
    def _setup_for_Solver(self):
        r"""
        Main features: Applying Boundary Conditions & Creating Transport Conductivities 

        This function executes the essential mathods for building matrices in Linear solution 
        """
        print '_setup_for_LinearSolver'

    
    def _do_outer_iteration_stage(self,**params):
        r"""
                       
        """
        print '_do_outer_iteration_stage'
        self._logger.info("Outer Iteration Stage ")
        self._setup_for_Solver(**params)       
        self._do_inner_iteration(**params)

    def _do_inner_iteration_stage(self,**params):
        r"""
                       
        """
        print '_do_outer_iteration_stage'

    def _do_one_inner_iteration(self):
        
        print '_do_inner_iteration_stage'
        A = self._coefficient_matrix()   
        B = self._RHS_matrix()
        x = splin.spsolve(A,B) 
        print x
        return(x)
    
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

        """
        setattr(self,"BCtypes",types)
        setattr(self,"BCvalues",values)      
           
    def _coefficient_matrix(self):
        
        # Filling coefficient matrix

        tpore1 = self._net.throat_properties['connections'][:,0]
        tpore2 = self._net.throat_properties['connections'][:,1]
        row = tpore1
        col = tpore2
        data= self._net.throat_properties['Conductivities_Exp']
        
        pnum = self._net.pore_properties['numbering']      
        loc1 = sp.in1d(tpore2,tpore2[sp.in1d(tpore2,pnum[self.BCtypes!=1])])
        modified_tpore2 = tpore2[loc1]
        modified_tpore1 = tpore1[loc1]        
        row = sp.append(row,modified_tpore2)                
        col = sp.append(col,modified_tpore1)
        data = sp.append(data,self._net.throat_properties['Conductivities_Exp'][loc1])
        
        if (self.BCtypes==4).any():
            self.extera_Nuemann_equations = sp.unique(self.BCvalues[self.BCtypes==2])
            self.Coeff_dimension = self._net.get_num_pores()+ len(self.extera_Nuemann_equations)
            extera_neu = self.extera_Nuemann_equations
            A_dim = self.Coeff_dimension
            g_super = 1e-8

            for item in range(len(extera_neu)):
                loc_neu = sp.in1d(tpore2,pnum[self.BCvalues==extera_neu[item]])
                neu_tpore2 = tpore2[loc_neu]                
                row = sp.append(row,neu_tpore2) 
                col = sp.append(col,sp.ones([len(neu_tpore2)])*(A_dim-item))
                data = sp.append(data,sp.ones([len(neu_tpore2)])*g_super)
                 
                loc_neu_row = loc_neu * loc1
                neu_tpore2_row = tpore2[loc_neu_row]
                row = sp.append(row,sp.ones([len(neu_tpore2_row)])*(A_dim-item)) 
                col = sp.append(col,neu_tpore2_row)
                data = sp.append(data,-(self._net.throat_properties['Conductivities_Exp'][neu_tpore2_row]))
            
        else:
            A_dim = self._net.get_num_pores()
        
        # Adding positions for diagonal
        
        row = sp.append(row,range(0,A_dim))
        col = sp.append(col,range(0,A_dim))
        data = sp.append(data,sp.zeros((A_dim,1)))

        a = sprs.coo.coo_matrix((data,(row,col)),(A_dim,A_dim))
        A = a.tocsr()        
       
        for i in range(0,A_dim):                    
            if self.BCtypes[i]==1:
                A[i,i] = 1                            
            elif i>self._net.get_num_pores():
                A[i,i] = len(A[i,:][sp.nonzero(A[i,:])])*g_super
            else:
                A[i,i] = -sp.sum(A[i,:][sp.nonzero(A[i,:])])
                
        return(A)
        
        
    def _RHS_matrix(self):
        
        extera_neu = self.extera_Nuemann_equations  
        A_dim = self.Coeff_dimension
        B = sp.zeros([A_dim,1])
        Dir_pores = self._net.pore_properties['numbering'][self.BCtypes==1]
        B[Dir_pores] = sp.reshape(self.BCvalues[Dir_pores],[len(Dir_pores),1])
        if (self.BCtypes==2).any():
            for item in len(extera_neu):                
                B[self.A_dim-item,0] = extera_neu[item]
#        else:
#            A_dim = self._net.get_num_pores()             
#            B = sp.reshape(self._net.pore_properties['BCvalue'],[A_dim,1])
#
#
#        for i in range(0,self._net.get_num_pores()):
#            if self._net.pore_properties['BCvalue'][i]!=0:
#                neighbors = self._net.get_neighbor_pores(i)
#                if sp.in1d(neighbors,self.Pinvaded).all():
#                    B[i,0] = 0