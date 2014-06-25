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

        
        self._logger.info("Creating Coefficient matrix for the algorithm")
        A = self._build_coefficient_matrix()
        self._logger.info("Creating RHS matrix for the algorithm")
        B = self._build_RHS_matrix()
        self._logger.info("Solving AX = B for the sparse matrices")
        X = sprslin.spsolve(A,B)
        self._Neumann_super_X = X[-sp.in1d(sp.r_[0:len(X)],sp.r_[0:self.num_pores()])]
        self._result = X[sp.r_[0:self.num_pores()]]        
        return(self._result)

    def set_boundary_conditions(self,bctype='',bcvalue=[],pores=[],mode='merge'):
        r"""
        Apply boundary conditions to specified pores.  This does not support
        throat boundary conditions yet.
        
        Parameters
        ----------
        bctype : string
            Specifies the type of boundary condition to apply.  Can be one of:
            
            - 'Dirichlet' : Specify the quantity in each pore
            - 'Neumann' : Specify the flow rate into each pore
            - 'Neumann_group' : Specify the net flow rate into a group of pores
            
        bcvalue : array_like
            The boundary value to apply, such as concentration or rate
        pores : array_like
            The pores where the boundary conditions should be applied
        mode : string, optional
            Controls how the conditions are applied.  Options are:
            
            - 'merge': Inserts the specified values, leaving existing values elsewhere
            - 'overwrite': Inserts specified values, clearing all other values
            - 'remove': Removes boundary conditions from specified pores
            - 'clear_all': Removes ALL boundary conditions
            
        Notes
        -----
        1. At the moment is it not possible to have multiple boundary conditions 
        in the same pore, so when new conditions are applied any existing ones
        are removed from all other boundary types.
        2. It is also not yet possible to apply boundary conditions to throats.
        """
        
        BC_default = ['Dirichlet','Neumann','Neumann_group']

        #If mode is 'clear_all' then bypass checks
        if mode == 'clear_all':
            for item in self.labels():
                bcname = item.split('.')[1]
                if bcname in BC_default:
                    self._logger.debug('Removing '+bcname)
                    del self['pore.bcval_'+bcname]
                    del self['pore.'+bcname]
            return

        #Validate bctype
        if bctype.split('.')[-1] not in BC_default:
            raise Exception('Unrecognized bctype') 
            
        #Validate pores
        if pores == []:
            if mode not in ['remove','clear_all']:
                raise Exception('Pores must be specified')
        else:
            pores = sp.array(pores,ndmin=1)
                    
        #Validate bcvalue
        if bcvalue == []:
            if mode not in ['remove','clear_all']:
                raise Exception('bcvalue must be specified')
        else:                        
            bcvalue = sp.array(bcvalue,ndmin=1)
        
        #Check bcvalues are compatible with bctypes
        if bctype == 'Neumann_group':  #Only scalars are acceptable
            if sp.size(bcvalue) != 1: 
                raise Exception('When specifying Neumann_group, bcval should be a scalar')
        else: #Only scalars or Np-long are acceptable
            if sp.size(bcvalue) == 1:
                bcvalue = sp.ones(sp.shape(pores))*bcvalue
            elif sp.size(bcvalue) != sp.size(pores):
                raise Exception('The pore list and bcvalue list are different lengths')
        
        #Confirm that prop and label arrays exist
        if 'pore.bcval_'+bctype not in self.props():
            self['pore.bcval_'+bctype] = sp.ones((self.num_pores(),),dtype=float)*sp.nan
        if 'pore.'+bctype not in self.labels():
            self['pore.'+bctype] = sp.zeros((self.num_pores(),),dtype=bool)
            
        #Remove all BC from specified pores, prior to setting new ones
        for item in self.labels():
            bcname = item.split('.')[-1]
            if bcname in BC_default:
                self['pore.bcval_'+bcname][pores] = sp.nan
                self['pore.'+bcname][pores] = False
        
        #Set boundary conditions based on supplied mode
        if mode == 'merge':
            self['pore.bcval_'+bctype][pores] = bcvalue
            self['pore.'+bctype][pores] = True
        elif mode == 'overwrite':
            self['pore.bcval_'+bctype] = sp.ones((self.num_pores(),),dtype=float)*sp.nan
            self['pore.bcval_'+bctype][pores] = bcvalue
            self['pore.'+bctype] = sp.zeros((self.num_pores(),),dtype=bool)
            self['pore.'+bctype][pores] = True
        elif mode == 'remove':
            if pores == []:
                self._logger.debug('Removing '+bctype)
                del self['pore.bcval_'+bctype]
                del self['pore.'+bctype]
            else:
                self._logger.debug('Removing '+bctype+' from specified pores')
                self['pore.bcval_'+bctype][pores] = sp.nan
                self['pore.'+bctype][pores] = False

    def _build_coefficient_matrix(self):
       
        # Filling coefficient matrix
        tpore1 = self._net.get_throat_data(prop='conns')[:,0]
        tpore2 = self._net.get_throat_data(prop='conns')[:,1]

        try:            
            Dir_pores = self.pores('Dirichlet')
            non_Dir_pores = self.pores('Dirichlet',mode='difference')
            loc1 = sp.in1d(tpore1,non_Dir_pores)
            loc2 = sp.in1d(tpore2,non_Dir_pores)
        except: 
            loc1 = sp.ones(len(tpore1),dtype=bool)
            loc2 = sp.ones(len(tpore2),dtype=bool)
        
        modified_tpore1 = tpore1[loc1]
        modified_tpore2 = tpore2[loc1]
        row = modified_tpore1
        col = modified_tpore2
        if sp.size(self._conductance)==1:
            self._conductance = self._conductance*sp.ones(self.num_throats())
        data_main = self._conductance
        try:
            self.pores('Neumann_insulated')
            insulated_pores = self.pores('Neumann_insulated')
            insulated_throats = self._net.find_neighbor_throats(insulated_pores,flatten=True,mode='not_intersection')
            data_main[insulated_throats] = 1e-60
        except: pass
        
        data = data_main[loc1]

        modified_tpore2 = tpore2[loc2]
        modified_tpore1 = tpore1[loc2]
        row = sp.append(row,modified_tpore2)
        col = sp.append(col,modified_tpore1)
        data = sp.append(data,data_main[loc2])

        A_dim = self.num_pores()

        try:
            self.pores('Neumann_group')
            group_values = self.get_data(prop='bcval_Neumann_group',pores=self.pores('Neumann_group'))
            self._group_Neumann_vals = sp.unique(group_values)
            A_dim = A_dim + len(self._group_Neumann_vals)
            extera_neu = self._group_Neumann_vals
            self._g_super = 1e-60            
            for item in sp.r_[0:len(extera_neu)]:
                neu_tpore2 = self.pores('Neumann_group')
                neu_tpore2 = neu_tpore2[group_values==extera_neu[item]]
                row = sp.append(row,neu_tpore2)
                col = sp.append(col,len(neu_tpore2)*[A_dim-item-1])
                data = sp.append(data,len(neu_tpore2)*[self._g_super])
                row = sp.append(row,len(neu_tpore2)*[A_dim-item-1])
                col = sp.append(col,neu_tpore2)
                data = sp.append(data,len(neu_tpore2)*[self._g_super])

        except: pass
        # Adding positions for diagonal
        dia = sp.r_[0:A_dim]
        try:
            Dir_pores
            row = sp.append(row,dia[Dir_pores])
            col = sp.append(col,dia[Dir_pores])
            data = sp.append(data,sp.ones_like(dia[Dir_pores]))
            temp_data = sp.copy(data)
            temp_data[sp.in1d(row,dia[Dir_pores])] = 0
            non_Dir_dia = dia[-sp.in1d(dia,dia[Dir_pores])]
        except:
            temp_data = sp.copy(data)
            non_Dir_dia = dia
        S_temp = sp.zeros(A_dim)
        for i in sp.r_[0:len(row)]:
            S_temp[row[i]] = S_temp[row[i]] - temp_data[i]
        data = sp.append(data,S_temp[non_Dir_dia])
        row = sp.append(row,non_Dir_dia)
        col = sp.append(col,non_Dir_dia)
        
        self._Coeff_dimension = A_dim
        a = sprs.coo.coo_matrix((data,(row,col)),(A_dim,A_dim))
        A = a.tocsr()
        return(A)


    def _build_RHS_matrix(self):        
        
        A_dim = self._Coeff_dimension
        B = sp.zeros([A_dim,1])
        try:
            Dir_pores = self.pores('Dirichlet')
            Dir_pores_vals = self.get_data(prop='bcval_Dirichlet',pores=Dir_pores)
            B[Dir_pores] = sp.reshape(Dir_pores_vals,[len(Dir_pores),1])
        except: pass
        try:
            individual_Neu_pores = self.pores('Neumann')
            individual_Neu_pores_vals = self.get_data(prop='bcval_Neumann',pores=individual_Neu_pores)
            B[individual_Neu_pores] = sp.reshape(individual_Neu_pores_vals,[len(individual_Neu_pores),1])
        except: pass
        try:
            self.pores('Neumann_group')
            pnum = self._net.num_pores()
            B[sp.r_[pnum:(pnum+len(self._group_Neumann_vals))]] = sp.reshape(self._group_Neumann_vals[sp.r_[0:len(self._group_Neumann_vals)]],[len(self._group_Neumann_vals),1])
        except: pass
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
        
        if 'pore.Dirichlet' in self:
            self._dir = self.get_pore_info(label='Dirichlet')
            del self['pore.Dirichlet']
        if 'pore.BCval' in self:
            self._BCval_temp = self.get_pore_data(prop='BCval')
            del self['pore.BCval']
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