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

import scipy as sp
import scipy.sparse as sprs
import scipy.sparse.linalg as splin
from .__GenericAlgorithm__ import GenericAlgorithm


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

        if (self.BCtypes==0).all():
            raise Exception('No boundary condition has been applied to this network.')
            self._result = sp.zeros(self._net.get_num_pores())
        else:
            A = self._build_coefficient_matrix()
            B = self._build_RHS_matrix()
            X = splin.spsolve(A,B)
            self._result = X[sp.r_[0:self._net.get_num_pores()]]
        return(self._result)

    def _boundary_conditions_setup(self,types=[],values=[]):
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
        - In Fickian algorithm, positive value for Neumann_rate or Neumann_flux for a pore with boundary condition means
        that the quantity of interest leaves the pore, but for any other algorithms, positive Neumann value  means 
        that the quantity of interest enters this pore.

        """
        self.BCtypes = sp.zeros(self._net.get_num_pores())
        self.BCvalues = sp.zeros(self._net.get_num_pores())
        for bctype in self._pore_info.keys():
            if bctype=='Dirichlet':
                self.BCtypes[self.get_pore_info(prop='Dirichlet')] = 1
                self.BCvalues[self.get_pore_info(prop='Dirichlet')] = self.get_pore_data(subdomain='Dirichlet',prop='BCval')
            elif bctype=='Neumann_flux':
                self.BCtypes[self.get_pore_info(prop='Neumann_flux')] = 2
                self.BCvalues[self.get_pore_info(prop='Neumann_flux')] = self.get_pore_data(subdomain='Neumann_flux',prop='BCval')                
            elif bctype=='Neumann_insulated':
                self.BCtypes[self.get_pore_info(prop='Neumann_insulated')] = 3
                self.BCvalues[self.get_pore_info(prop='Neumann_insulated')] = self.get_pore_data(subdomain='Neumann_insulated',prop='BCval') 
            elif bctype=='Neumann_rate':
                self.BCtypes[self.get_pore_info(prop='Neumann_rate')] = 4
                self.BCvalues[self.get_pore_info(prop='Neumann_rate')] = self.get_pore_data(subdomain='Neumann_rate',prop='BCval') 



    def _build_coefficient_matrix(self):
       
        boundaries = self._net.get_pore_data(prop='numbering')[self._net.get_pore_info(prop='boundary')]
        if (self.BCtypes[boundaries]==0).any():
            self.BCtypes[boundaries[self.BCtypes[boundaries]==0]] = 3
        
        # Filling coefficient matrix
        pnum = self._net.get_pore_data(prop='numbering')
        tpore1 = self._net.get_throat_data(prop='connections')[:,0]
        tpore2 = self._net.get_throat_data(prop='connections')[:,1]

        loc1 = sp.in1d(tpore1,pnum[self.BCtypes!=1])
        modified_tpore1 = tpore1[loc1]
        modified_tpore2 = tpore2[loc1]
        row = modified_tpore1
        col = modified_tpore2
        if sp.size(self._conductance)==1:
            self._conductance = self._conductance*sp.ones(self._net.get_num_throats())
        data_main = self._conductance
        data = data_main[loc1]

        loc2 = sp.in1d(tpore2,pnum[self.BCtypes!=1])
        modified_tpore2 = tpore2[loc2]
        modified_tpore1 = tpore1[loc2]
        row = sp.append(row,modified_tpore2)
        col = sp.append(col,modified_tpore1)
        data = sp.append(data,data_main[loc2])

        A_dim = self._net.get_num_pores()

        if (self.BCtypes==2).any():
            flux_pores = self._net.get_pore_data(prop='numbering')[self.BCtypes==2]
            flux_values = sp.unique(self.BCvalues[self.BCtypes==2])
            for i in list(range(len(flux_values))):
                f = flux_pores[sp.in1d(flux_pores,self._net.get_pore_data(prop='numbering')[self.BCvalues==flux_values[i]])]
                fn = self._net.get_neighbor_pores(f)
                fn = fn[self._net.get_pore_info(prop='internal')[fn]]
                ft = self._net.get_connecting_throat(f,fn)
                self.BCtypes[f] = 4
                self.BCvalues[f] = sp.sum(self.BCvalues[f]*(self._net.get_throat_data(prop='diameter')[ft])**2)

        if (self.BCtypes==4).any():
            self.extera_Neumann_equations = sp.unique(self.BCvalues[self.BCtypes==4])
            A_dim = A_dim + len(self.extera_Neumann_equations)
            extera_neu = self.extera_Neumann_equations
            g_super = sp.average(self._conductance[self._net.get_neighbor_throats(pnum[self.BCtypes==4])])*1e5
            for item in list(range(len(extera_neu))):
                neu_tpore2 = pnum[self.BCvalues==extera_neu[item]]

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

        dia = sp.array(list(range(0,A_dim)))
        row = sp.append(row,dia[self.BCtypes==1])
        col = sp.append(col,dia[self.BCtypes==1])
        data = sp.append(data,sp.ones_like(dia[self.BCtypes==1]))

        temp_data = sp.zeros(A_dim-len(dia[self.BCtypes==1]))
        non_Dir = dia[-sp.in1d(dia,dia[self.BCtypes==1])]
        for i in list(range(len(non_Dir))):
            temp_data[i] = -sp.sum(data[row==non_Dir[i]])
        data = sp.append(data,temp_data)
        row = sp.append(row,non_Dir)
        col = sp.append(col,non_Dir)

        a = sprs.coo.coo_matrix((data,(row,col)),(A_dim,A_dim))
        A = a.tocsr()

        return(A)


    def _build_RHS_matrix(self):

        extera_neu = self.extera_Neumann_equations
        A_dim = self.Coeff_dimension
        B = sp.zeros([A_dim,1])
        Dir_pores = self._net.get_pore_data(prop='numbering')[self.BCtypes==1]
        B[Dir_pores] = sp.reshape(self.BCvalues[Dir_pores],[len(Dir_pores),1])
        if (self.BCtypes==4).any():
            for item in list(range(len(extera_neu))):
                B[A_dim-item-1,0] = extera_neu[item]

        return(B)

    def rate(self,pores=[],throats=[]):

        if throats:
            pores1 = self._net.get_connected_pores(throats)[:,0]
            pores2 = self._net.get_connected_pores(throats)[:,1]
        else:            
            throats = self._net.get_neighbor_throats(pores,flatten=True)
            pores1 = self._net.get_connected_pores(throats)[:,0]
            pores2 = self._net.get_connected_pores(throats)[:,1]
        X1 = self._result[pores1]
        X2 = self._result[pores2]
        g = self._conductance[throats]
        R = sp.sum(sp.multiply(g,(X1-X2)))
        return(R)