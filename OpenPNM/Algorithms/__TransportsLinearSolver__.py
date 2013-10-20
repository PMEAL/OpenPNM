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

class TransportsLinearSolver(GenericAlgorithm):
    r"""   

        
    """
    
    def __init__(self,net,loglevel=10,**kwargs):
        r"""
        Initializing the class
        """
        super(TransportsLinearSolver,self).__init__(net = net,**kwargs)
        self._logger.info("Create Linear Solver Algorithm Object")
        print 'init'
            
    def _setup_for_TransportSolver(self):
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
        self._setup_for_TransportSolver(**params)       
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
        
           
    def _Coefficient_Matrix(self):
        
        # Filling coefficient matrix

        tpore1 = self._net.throat_properties['connections'][:,0]
        tpore2 = self._net.throat_properties['connections'][:,1]
        row = tpore1
        col = tpore2
        data= self._net.throat_properties['Conductivities_Exp']
        
        pnum = self._net.pore_properties['numbering']      
        loc1 = sp.in1d(tpore2,tpore2[sp.in1d(tpore2,pnum[self._net.pore_properties['BCtype']!=1])])
        modified_tpore2 = tpore2[loc1]
        modified_tpore1 = tpore1[loc1]        
        row = sp.append(row,modified_tpore2)                
        col = sp.append(col,modified_tpore1)
        data = sp.append(data,self._net.throat_properties['Conductivities_Exp'][loc1])
        
        if (self._net.pore_properties['BCtype']==2).any():
            self.extera_Nuemann_equations = sp.unique(self._net.pore_properties['BCvalue'][self._net.pore_properties['BCtype']==2])
            self.Coeff_dimension = self._net.get_num_pores()+ len(self.extera_Nuemann_equations)
            extera_neu = self.extera_Nuemann_equations
            A_dim = self.Coeff_dimension
            g_super = 1e-8

            for item in range(len(extera_neu)):
                loc_neu = sp.in1d(tpore2,pnum[self._net.pore_properties['BCvalue']==extera_neu[item]])
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
            if self._net.pore_properties['BCtype'][i]==1:
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
        Dir_pores = self._net.pore_properties['numbering'][self._net.pore_properties['BCtype']==1]
        B[Dir_pores] = sp.reshape(self._net.pore_properties['BCvalue'][Dir_pores],[len(Dir_pores),1])
        if (self._net.pore_properties['BCtype']==2).any():
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