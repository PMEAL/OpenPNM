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

        if (self._BCtypes==0).all():
            raise Exception('No boundary condition has been applied to this network.')
            self._result = sp.zeros(self._net.num_pores())
        else:
            A = self._build_coefficient_matrix()
            B = self._build_RHS_matrix()
            X = splin.spsolve(A,B)
            self._result = X[sp.r_[0:self._net.num_pores()]]
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
        self._BCtypes = sp.zeros(self._net.num_pores())
        self._BCvalues = sp.zeros(self._net.num_pores())
        for bctype in self._pore_info.keys():
            if bctype=='Dirichlet':
                bcpores = self.get_pore_info(label='Dirichlet')
                self._BCtypes[bcpores] = 1 
                self._BCvalues[bcpores] = self.get_pore_data(locations=bcpores,prop='BCval')
            elif bctype=='Neumann_flux':
                bcpores = self.get_pore_info(label='Neumann_flux')
                self._BCtypes[bcpores] = 2
                self._BCvalues[bcpores] = self.get_pore_data(locations=bcpores,prop='BCval')
            elif bctype=='Neumann_insulated':
                bcpores = self.get_pore_info(label='Neumann_insulated')
                self._BCtypes[bcpores] = 3
                self._BCvalues[bcpores] = self.get_pore_data(locations=bcpores,prop='BCval')
            elif bctype=='Neumann_rate':
                bcpores = self.get_pore_info(label='Neumann_rate')
                self._BCtypes[bcpores] = 4
                self._BCvalues[bcpores] = self.get_pore_data(locations=bcpores,prop='BCval')


    def _build_coefficient_matrix(self):
       
        boundaries = self._net.get_pore_indices()[-sp.in1d(self._net.get_pore_indices(),self._net.get_pore_indices('internal'))]
        if (self._BCtypes[boundaries]==0).any():
            self._BCtypes[boundaries[self._BCtypes[boundaries]==0]] = 3
        
        # Filling coefficient matrix
        pnum = self._net.get_pore_indices()
        tpore1 = self._net.get_throat_data(prop='connections')[:,0]
        tpore2 = self._net.get_throat_data(prop='connections')[:,1]

        loc1 = sp.in1d(tpore1,pnum[self._BCtypes!=1])
        modified_tpore1 = tpore1[loc1]
        modified_tpore2 = tpore2[loc1]
        row = modified_tpore1
        col = modified_tpore2
        if sp.size(self._conductance)==1:
            self._conductance = self._conductance*sp.ones(self._net.num_throats())
        data_main = self._conductance
        data = data_main[loc1]

        loc2 = sp.in1d(tpore2,pnum[self._BCtypes!=1])
        modified_tpore2 = tpore2[loc2]
        modified_tpore1 = tpore1[loc2]
        row = sp.append(row,modified_tpore2)
        col = sp.append(col,modified_tpore1)
        data = sp.append(data,data_main[loc2])

        A_dim = self._net.num_pores()

        if (self._BCtypes==2).any():
            flux_pores = self._net.get_pore_indices()[self._BCtypes==2]
            flux_values = sp.unique(self._BCvalues[self._BCtypes==2])
            for i in list(range(len(flux_values))):
                f = flux_pores[sp.in1d(flux_pores,self._net.get_pore_indices()[self._BCvalues==flux_values[i]])]
                fn = self._net.find_neighbor_pores(f)
                fn = fn[self._net.get_pore_info(label='internal')[fn]]
                ft = self._net.find_connecting_throat(f,fn)
                self._BCtypes[f] = 4
                self._BCvalues[f] = sp.sum(self._BCvalues[f]*(self._net.get_throat_data(prop='diameter')[ft])**2)

        if (self._BCtypes==4).any():
            self._extera_Neumann_equations = sp.unique(self._BCvalues[self._BCtypes==4])
            A_dim = A_dim + len(self._extera_Neumann_equations)
            extera_neu = self._extera_Neumann_equations
            g_super = sp.average(self._conductance[self._net.find_neighbor_throats(pnum[self._BCtypes==4])])*1e5
            for item in list(range(len(extera_neu))):
                neu_tpore2 = pnum[self._BCvalues==extera_neu[item]]

                row = sp.append(row,neu_tpore2)
                col = sp.append(col,len(neu_tpore2)*[A_dim-item-1])
                data = sp.append(data,len(neu_tpore2)*[g_super])

                row = sp.append(row,len(neu_tpore2)*[A_dim-item-1])
                col = sp.append(col,neu_tpore2)
                data = sp.append(data,len(neu_tpore2)*[g_super])

        else:
            self._extera_Neumann_equations = 0

        self._Coeff_dimension = A_dim

        # Adding positions for diagonal

        dia = sp.array(list(range(0,A_dim)))
        row = sp.append(row,dia[self._BCtypes==1])
        col = sp.append(col,dia[self._BCtypes==1])
        data = sp.append(data,sp.ones_like(dia[self._BCtypes==1]))

        temp_data = sp.zeros(A_dim-len(dia[self._BCtypes==1]))
        non_Dir = dia[-sp.in1d(dia,dia[self._BCtypes==1])]
        for i in list(range(len(non_Dir))):
            temp_data[i] = -sp.sum(data[row==non_Dir[i]])
        data = sp.append(data,temp_data)
        row = sp.append(row,non_Dir)
        col = sp.append(col,non_Dir)

        a = sprs.coo.coo_matrix((data,(row,col)),(A_dim,A_dim))
        A = a.tocsr()

        return(A)


    def _build_RHS_matrix(self):

        extera_neu = self._extera_Neumann_equations
        A_dim = self._Coeff_dimension
        B = sp.zeros([A_dim,1])
        Dir_pores = self._net.get_pore_indices()[self._BCtypes==1]
        B[Dir_pores] = sp.reshape(self._BCvalues[Dir_pores],[len(Dir_pores),1])
        if (self._BCtypes==4).any():
            for item in list(range(len(extera_neu))):
                B[A_dim-item-1,0] = extera_neu[item]

        return(B)

    def rate(self,pores=[],throats=[]):

        if throats:
            pores1 = self._net.find_connected_pores(throats)[:,0]
            pores2 = self._net.find_connected_pores(throats)[:,1]
        else:            
            throats = self._net.find_neighbor_throats(pores,flatten=True)
            pores1 = self._net.find_connected_pores(throats)[:,0]
            pores2 = self._net.find_connected_pores(throats)[:,1]
        X1 = self._result[pores1]
        X2 = self._result[pores2]
        g = self._conductance[throats]
        R = sp.sum(sp.multiply(g,(X1-X2)))
        return(R)
        
    def _calc_eff_prop_cubic(self,                            
                           fluid,
                           alg='',
                           face1='',
                           face2='',
                           d_term='',
                           x_term='',
                           conductance='',
                           occupancy='',
                           **params):
                
        network =self._net
        ftype1 = []
        ftype2 = []
        effective_prop = []
        if face1!='' and face2!='':             
            if type(face1)==list and type(face2)==list: 
                if len(face1)==len(face2):
                    for i in list(range(len(face1))):
                        if type(face1[i])==str and type(face2[i])==str:
                            ftype1.append(face1[i])
                            ftype2.append(face2[i])                            
                        else: self._logger.error('face1 and face2 should be string or a list of strings!')
                else: self._logger.error('face1 and face 2 should have equal length!')

            elif type(face1)==str and type(face2)==str:
                ftype1.append(face1)
                ftype2.append(face2)
            else: self._logger.error('face1 and face2 should be string or a list of strings!')

        elif face1=='' and face2=='':            
            ftype1 = ['front','right','top']
            ftype2 = ['back','left','bottom']
        else: self._logger.error('wrong input for face1 or face2')
        
        for i in list(range(len(ftype1))):
            face1 = ftype1[i] 
            face2 = ftype2[i]
            face1_pores = network.get_pore_indices(face1)
            face2_pores = network.get_pore_indices(face2)            
            ## Assign Dirichlet boundary conditions
            ## BC1
            BC1_pores = face1_pores  
            self.set_pore_info(label='Dirichlet',locations=BC1_pores)
            BC1_values = 0.8
            self.set_pore_data(prop='BCval',data=BC1_values,locations=BC1_pores)
            ## BC2
            BC2_pores = face2_pores
            self.set_pore_info(label='Dirichlet',locations=BC2_pores)
            BC2_values = 0.4
            self.set_pore_data(prop='BCval',data=BC2_values,locations=BC2_pores)        
            self.run(active_fluid=fluid) 
            x = self.get_pore_data(prop=x_term)
            if alg=='Fickian':
                X1 = sp.log(1-x[face1_pores])
                X2 = sp.log(1-x[face2_pores])
            elif alg=='Stokes':
                X1 = x[face1_pores]
                X2 = x[face2_pores]
            delta_X = sp.absolute(sp.average(X2)-sp.average(X1)) 
            d_force =sp.average(fluid.get_pore_data(prop=d_term))
            
            coordx = network.get_pore_data(prop='coords')[:,0]
            coordy = network.get_pore_data(prop='coords')[:,1]
            coordz = network.get_pore_data(prop='coords')[:,2]
            
            if sp.size(sp.unique(coordx[network.get_pore_indices(face1)]))==1:
                coord_main1 = coordx[network.get_pore_indices(face1)]
                coord_main2 = coordx[network.get_pore_indices(face2)]
                coord_temp1 = coordy[network.get_pore_indices(face1)]
                coord_temp2 = coordz[network.get_pore_indices(face1)]
             
            elif sp.size(sp.unique(coordy[network.get_pore_indices(face1)]))==1:
                coord_main1 = coordy[network.get_pore_indices(face1)]
                coord_main2 = coordy[network.get_pore_indices(face2)]
                coord_temp1 = coordx[network.get_pore_indices(face1)]
                coord_temp2 = coordz[network.get_pore_indices(face1)]
                
            elif sp.size(sp.unique(coordz[network.get_pore_indices(face1)]))==1:
                coord_main1 = coordz[network.get_pore_indices(face1)]
                coord_main2 = coordz[network.get_pore_indices(face2)]
                coord_temp1 = coordx[network.get_pore_indices(face1)]
                coord_temp2 = coordy[network.get_pore_indices(face1)]

            L = sp.absolute(sp.unique(coord_main1)[0]-sp.unique(coord_main2)[0])
            length_1 = (max(coord_temp1) - min(coord_temp1))
            length_2 = (max(coord_temp2) - min(coord_temp2))
            A = length_1*length_2            
                
            fn = network.find_neighbor_pores(face1_pores)
            fn = fn[sp.in1d(fn,network.get_pore_indices('internal'))]
            ft = network.find_connecting_throat(face1_pores,fn)
            if alg=='Fickian': X_temp = sp.log(1-x[fn])
            elif alg=='Stokes':
                X_temp = x[fn]
                d_force = 1/d_force
            g = fluid.get_throat_data(prop=conductance)
            if occupancy=='': s = sp.ones_like(g, dtype=bool)
            elif occupancy=='occupancy': s = fluid.get_throat_data(prop=occupancy)
            cond = g*s+g*(-s)/1e3
            N = sp.sum(cond[ft]*sp.absolute(X1-X_temp))
            eff = N*L/(d_force*A*delta_X)
            effective_prop.append(eff)
            del self._pore_info['Dirichlet']
            del self._pore_data['BCval']
            delattr (self,'BCtypes')
            delattr(self,'BCvalues')            
        return sp.array(effective_prop,ndmin=1)