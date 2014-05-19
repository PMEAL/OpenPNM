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
import scipy.linalg as splin
import scipy.sparse.linalg as sprslin
from .__GenericAlgorithm__ import GenericAlgorithm


class LinearSolver(GenericAlgorithm):
    r"""
    This class provides essential methods for building and solving matrices in a transport process.
    It will be inherited during Fickian diffusion, Fourier heat conduction, Hagen Poiseuille permeability and electron conduction.

    """

    def __init__(self,**kwargs):
        r"""
        Initializing the class
        """
        super(LinearSolver,self).__init__(**kwargs)

    def _do_one_inner_iteration(self):

        if (self._BCtypes==0).all():
            raise Exception('No boundary condition has been applied to this network.')
            self._result = sp.zeros(self._net.num_pores())
        else:
            self._logger.info("Creating Coefficient matrix for the algorithm")
            A = self._build_coefficient_matrix()
            self._logger.info("Creating RHS matrix for the algorithm")
            B = self._build_RHS_matrix()
            self._logger.info("Solving AX = B for the sparse matrices")
            X = sprslin.spsolve(A,B)
            self._Neumann_super_X = X[-sp.in1d(sp.r_[0:len(X)],sp.r_[0:self._net.num_pores()])]
            self._result = X[sp.r_[0:self._net.num_pores()]]        
        return(self._result)

    def set_boundary_conditions(self,bctype='',bcvalue=[],pores=[],throats=[],mode='merge'):
        r"""
        
        """
        BC_default = ['Dirichlet','Neumann_insulated','Neumann_rate_group','Neumann_rate_single']
        if pores==[] and throats==[]:  
            self._logger.error('No pore/throat has been assigned for this boundary condition!') 
            setup = 0
        else:
            elements =[]
            if pores!= []: elements.append('pore')
            if throats != []: elements.append('throat')
            for element in elements:                
                try:
                    getattr(self,'_BCtypes_'+element)
                    getattr(self,'_BCvalues_'+element)
                except: 
                    setattr(self,'_BCtypes_'+element, sp.zeros(getattr(self,'num_'+element+'s')()))
                    setattr(self,'_BCvalues_'+element,sp.zeros(getattr(self,'num_'+element+'s')()))
                existing_bc = []
                temp ='None'
                if element=='pore':
                    if pores=='all':    
                        loc = self.pores()
                        temp = 'all'
                    else:   loc = pores
                elif    element=='throat':
                    if throats=='all':    
                        loc = self.throats()
                        temp = 'all'
                    else:   loc = throats 
                for label in getattr(self,'_'+element+'_info').keys():
                    if label in BC_default and label not in existing_bc:    existing_bc.append(label)
                if mode=='remove':
                    getattr(self,'_BCtypes_'+element)[loc] = 0
                    getattr(self,'_BCvalues_'+element)[loc] = 0
                    getattr(self,'_set_data')(element=element,prop='BCval',locations=loc,mode='remove')
                    for bc_type in existing_bc:
                        if temp=='all':
                            getattr(self,'_set_info')(element=element,label=bc_type,mode='remove')
                        else:
                            getattr(self,'_set_info')(element=element,label=bc_type,locations=loc,mode='remove')
                    
                    if not (bctype=='' and bcvalue==[]):
                        self._logger.info('To remove boundary conditions from some locations, no value or type should be sent!')
                    setup = 1
                else:
                    if bctype in BC_default:
                        if bctype=='Dirichlet': bc_num = 1                    
                        elif bctype=='Neumann_insulated': bc_num = 2
                        elif bctype=='Neumann_rate_group': bc_num = 3
                        elif bctype=='Neumann_rate_single': bc_num = 4
                        if mode=='overwrite':
                            setattr(self,'_BCtypes_'+element, sp.zeros(getattr(self,'num_'+element+'s')()))
                            setattr(self,'_BCvalues_'+element,sp.zeros(getattr(self,'num_'+element+'s')()))
                            getattr(self,'_set_info')(element=element,label=bc_type,locations=loc,mode='overwrite')
                            getattr(self,'_set_data')(element=element,prop='BCval',locations=loc,mode='remove')
                            self._logger.info('Boundary conditions have been overwritten for the algorithm: '+self.name)
                            setup = 1
                        elif mode=='merge':                         
                            if (getattr(self,'_BCtypes_'+element)[loc] == 0).all():     setup = 1
                            else:
                                ind = loc[getattr(self,'_BCtypes_'+element)[loc] != 0]
                                self._logger.error('Boundary conditions have already been assigned to the '+element+'s: '+str(ind))
                                self._logger.info('To apply new bounday conditions to these locations, the existing BCs should be removed.')
                                setup = 0
                        if setup==1:
                            getattr(self,'_BCtypes_'+element)[loc] = bc_num
                            if bc_num!=2:
                                getattr(self,'_BCvalues_'+element)[loc] = bcvalue
                                getattr(self,'_set_data')(element=element,prop='BCval',locations=loc,mode='merge',data=getattr(self,'_BCvalues_'+element)[loc])
                                                  
                    else:
                        self._logger.error('The bctype: '+bctype+' has not been defined for the algorithm!')
            try:
                getattr(self,'_BCtypes_throat')
                setup = 0
                self._logger.error('The section for assigning throat boundary conditions is not implemented in this solver yet.')
            except:
                self._BCtypes = getattr(self,'_BCtypes_pore')
                self._BCvalues = getattr(self,'_BCvalues_pore')
        self.bc_setup = setup

    def _build_coefficient_matrix(self):
       
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
        if (self._BCtypes==2).any():
            insulated_pores = self._net.get_pore_indices()[self._BCtypes==2]
            insulated_throats = self._net.find_neighbor_throats(insulated_pores,flatten=True,mode='not_intersection')
            data_main[insulated_throats] = 1e-60
        data = data_main[loc1]

        loc2 = sp.in1d(tpore2,pnum[self._BCtypes!=1])
        modified_tpore2 = tpore2[loc2]
        modified_tpore1 = tpore1[loc2]
        row = sp.append(row,modified_tpore2)
        col = sp.append(col,modified_tpore1)
        data = sp.append(data,data_main[loc2])

        A_dim = self._net.num_pores()
           
        if (self._BCtypes==3).any():
            self._extera_Neumann_equations = sp.unique(self._BCvalues[self._BCtypes==3])
            A_dim = A_dim + len(self._extera_Neumann_equations)
            extera_neu = self._extera_Neumann_equations
            self._g_super = 1e-60            
            mask = self._BCtypes==3
            for item in sp.r_[0:len(extera_neu)]:
                neu_tpore2 = pnum[mask]
                neu_tpore2 = neu_tpore2[self._BCvalues[neu_tpore2]==extera_neu[item]]
                row = sp.append(row,neu_tpore2)
                col = sp.append(col,len(neu_tpore2)*[A_dim-item-1])
                data = sp.append(data,len(neu_tpore2)*[self._g_super])
                row = sp.append(row,len(neu_tpore2)*[A_dim-item-1])
                col = sp.append(col,neu_tpore2)
                data = sp.append(data,len(neu_tpore2)*[self._g_super])

        else:
            self._extera_Neumann_equations = 0

        self._Coeff_dimension = A_dim

        # Adding positions for diagonal
        dia = sp.r_[0:A_dim]
        row = sp.append(row,dia[self._BCtypes==1])
        col = sp.append(col,dia[self._BCtypes==1])
        data = sp.append(data,sp.ones_like(dia[self._BCtypes==1]))

        temp_data = sp.copy(data)
        temp_data[sp.in1d(row,dia[self._BCtypes==1])] = 0
        S_temp = sp.zeros(A_dim)
        for i in sp.r_[0:len(row)]:
            S_temp[row[i]] = S_temp[row[i]] - temp_data[i]
        non_Dir = dia[-sp.in1d(dia,dia[self._BCtypes==1])]
        data = sp.append(data,S_temp[non_Dir])
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
        individual_Neu_pores = self._net.get_pore_indices()[self._BCtypes==4]
        B[individual_Neu_pores] = sp.reshape(self._BCvalues[individual_Neu_pores],[len(individual_Neu_pores),1])
        if (self._BCtypes==3).any():
            for item in sp.r_[0:len(extera_neu)]:
                B[A_dim-item-1,0] = extera_neu[item]
            
        return(B)

    def rate(self,pores='',throats=''):
        r'''
        Send a list of pores (or throats) and recieve the cumulative rate
        of material moving into them
        '''

        if throats!='':
            p1 = self._net.find_connected_pores(throats)[:,0]
            p2 = self._net.find_connected_pores(throats)[:,1]
        elif pores!='': 
            throats = self._net.find_neighbor_throats(pores,flatten=True,mode='not_intersection')
            p1 = self._net.find_connected_pores(throats)[:,0]
            p2 = self._net.find_connected_pores(throats)[:,1]
        pores1 = sp.copy(p1)
        pores2 = sp.copy(p2)
        pores1[-sp.in1d(p1,pores)] = p2[-sp.in1d(p1,pores)]        
        pores2[-sp.in1d(p1,pores)] = p1[-sp.in1d(p1,pores)]
        X1 = self._result[pores1]
        X2 = self._result[pores2]
        g = self._conductance[throats]
        R = sp.sum(sp.multiply(g,(X1-X2)))
        return(R)
        
    def _calc_eff_prop(self,                            
                       fluid,
                       alg,
                       d_term,
                       x_term,
                       conductance,
                       occupancy,
                       direction,
                       **params):
                
        network =self._net
        ftype1 = []
        ftype2 = []
        effective_prop = []  
        result = {}
        try: fluid = self.find_object_by_name(fluid) 
        except: pass #Accept object
        if type(direction)==str and direction=='': 
            ftype1 = ['front','right','top']
            ftype2 = ['back','left','bottom']            
        elif type(direction)==str: direction = sp.array(direction,ndmin=1)              
        if type(direction)==sp.ndarray:
            ftype1 = []
            ftype2 = []
            for d in direction:
                if d=='x' or d=='X' or d=='front' or d=='back': 
                    ftype1.append('front')
                    ftype2.append('back')
                elif d=='y' or d=='Y'or d=='left' or d=='right': 
                    ftype1.append('left')
                    ftype2.append('right')
                elif d=='z' or d=='Z'or d=='top' or d=='bottom': 
                    ftype1.append('top')
                    ftype2.append('bottom') 
                else: self._logger.error('wrong input for direction!')
        
        if 'Dirichlet' in self._pore_info:
            self._dir = self.get_pore_info(label='Dirichlet')
            del self._pore_info['Dirichlet']
        if 'BCval' in self._pore_data:
            self._BCval_temp = self.get_pore_data(prop='BCval')
            del self._pore_data['BCval']
            try:
                self._BCtypes_temp = sp.copy(self._BCtypes)
                delattr (self,'_BCtypes')
                self._BCvalues_temp = sp.copy(self._BCvalues)
                delattr(self,'_BCvalues')
            except: pass
        try: self._X_temp = self.get_pore_data(prop=self._X_name)  
        except: pass          
        tensor = sp.zeros([3,3])
        for i in sp.r_[0:len(ftype1)]:
            face1 = ftype1[i] 
            face2 = ftype2[i]
            if face1=='front' or face1=='back': direct = 'X'
            elif face1=='left' or face1=='right': direct = 'Y'
            elif face1=='top' or face1=='bottom': direct = 'Z'
            if 'boundary' in self._net._pore_info:
                face1_pores = network.get_pore_indices(labels=[face1,'boundary'],mode='intersection')
                face2_pores = network.get_pore_indices(labels=[face2,'boundary'],mode='intersection')
            else:    
                face1_pores = network.get_pore_indices(face1)
                face2_pores = network.get_pore_indices(face2)            
            ## Assign Dirichlet boundary conditions
            ## BC1
            BC1_pores = face1_pores  
            self.set_pore_info(label='Dirichlet',locations=BC1_pores,mode='overwrite')
            BC1_values = 0.8
            self.set_pore_data(prop='BCval',data=BC1_values,locations=BC1_pores)
            ## BC2
            BC2_pores = face2_pores
            self.set_pore_info(label='Dirichlet',locations=BC2_pores)
            BC2_values = 0.4
            self.set_pore_data(prop='BCval',data=BC2_values,locations=BC2_pores)        
            self.run(active_fluid=fluid,
                           x_term=x_term,
                           conductance=conductance,
                           occupancy=occupancy) 
            x = self.get_pore_data(prop=x_term)
            if alg=='Fickian':
                X1 = sp.log(1-x[face1_pores])
                X2 = sp.log(1-x[face2_pores])
            elif alg=='Stokes':
                X1 = x[face1_pores]
                X2 = x[face2_pores]
            delta_X = sp.absolute(sp.average(X2)-sp.average(X1)) 
            d_force =sp.average(fluid.get_pore_data(prop=d_term))

            if  face1=='top' or face1=='bottom': 
                L = self._net.domain_size('height')
                A = self._net.domain_size('top')
            elif  face1=='left' or face1=='right':
                L = self._net.domain_size('depth')
                A = self._net.domain_size('left')
            elif  face1=='front' or face1=='back':
                L = self._net.domain_size('width')
                A = self._net.domain_size('front')
            fn = network.find_neighbor_pores(face1_pores,excl_self=True)
            fn = fn[sp.in1d(fn,network.get_pore_indices('internal'))]
            ft = network.find_connecting_throat(face1_pores,fn)
            if alg=='Fickian': X_temp = sp.log(1-x[fn])
            elif alg=='Stokes':
                X_temp = x[fn]
                d_force = 1/d_force
            cond = self._conductance
            N = sp.sum(cond[ft]*sp.absolute(X1-X_temp))
            eff = N*L/(d_force*A*delta_X)
            effective_prop.append(eff)
            del self._pore_info['Dirichlet']
            del self._pore_data['BCval']
            delattr (self,'_BCtypes')
            delattr(self,'_BCvalues')            
            result[ftype1[i]+'/'+ftype2[i]+'('+direct+')'] = sp.array(effective_prop[i],ndmin=1)
            if ftype1[i]=='top' or ftype1[i]=='bottom': tensor[2,2] = effective_prop[i]
            elif ftype1[i]=='right' or ftype1[i]=='left': tensor[1,1] = effective_prop[i]
            elif ftype1[i]=='front' or ftype1[i]=='back': tensor[0,0] = effective_prop[i]
        
        try:
            self.set_pore_data(prop=self._X_name,data=self._X_temp)
            delattr (self,'_X_temp')
        except : del self._pore_data[self._X_name]
        try:
            self.set_pore_info(label='Dirichlet',locations=self._dir,mode='overwrite')
            delattr (self,'_dir')
            self.set_pore_data(prop='BCval',data=self._BCval_temp)
            delattr (self,'_BCval_temp')
            self._BCtypes = self._BCtypes_temp
            delattr (self,'_BCtypes_temp')
            self._BCvalues = self._BCvalues_temp
            delattr (self,'_BCvalues_temp')
        except: pass        
        if len(ftype1)<3: return result
        elif len(ftype1)==3 : return tensor