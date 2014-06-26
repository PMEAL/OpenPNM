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
        r'''
        Initializing the class
        '''
        super(LinearSolver,self).__init__(**kwargs)
        
    def setup(self,fluid,conductance,quantity):
        r'''
        '''
        self._fluid = fluid
        self._conductance = 'throat.'+conductance.split('.')[-1]
        self._quantity = 'pore.'+quantity.split('.')[-1]
        
        #Check health of conductance vector
        if self._fluid.check_throat_health(props=self._conductance):
            #If no nans, check for 0's
            temp = sp.where(self._fluid[self._conductance]==0)[0]
            self._fluid[self._conductance][temp] = 1e-30
        else:
            raise Exception('The provided throat conductance has problems')
            
    def update(self):
        r'''
        Send results of simulation out the the appropriate locations.
        
        This is a basic version of the update that simply sends out the main
        result (quantity). More elaborate updates should be subclassed.
        '''        
        self._fluid[self._quantity] = self[self._quantity]
        self._logger.debug('Results of '+self.name+' algorithm have been added to '+self._fluid.name)

    def set_boundary_conditions(self,bctype='',bcvalue=[],pores=[],mode='merge'):
        r'''
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
        '''
        
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
        r'''
        '''
        # Filling coefficient matrix
        tpore1 = self._net['throat.conns'][:,0]
        tpore2 = self._net['throat.conns'][:,1]
        
        #Identify Dirichlet pores, if any
        try:
            temp = self.pores('Dirichlet',mode='difference')
        except:
            temp = self.pores()
        loc1 = sp.in1d(tpore1,temp)
        loc2 = sp.in1d(tpore2,temp)
        modified_tpore1 = tpore1[loc1]
        modified_tpore2 = tpore2[loc1]
        row = modified_tpore1
        col = modified_tpore2
        
        #Expand the conductance to a vector if necessary
        g = self._fluid[self._conductance]
        if sp.size(g) == 1:
            g = g*sp.ones(self.num_throats())
        data_main = g
        data = data_main[loc1]

        modified_tpore2 = tpore2[loc2]
        modified_tpore1 = tpore1[loc2]
        row = sp.append(row,modified_tpore2)
        col = sp.append(col,modified_tpore1)
        data = sp.append(data,data_main[loc2])
        A_dim = self.num_pores()
        
        #Check for Neuman_group BCs and add superpores if necessary
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
        except: 
            pass
        
        # Adding positions for diagonal
        diag = sp.r_[0:A_dim]
        try:
            pores = self.pores('Dirichlet')
            row = sp.append(row,diag[pores])
            col = sp.append(col,diag[pores])
            data = sp.append(data,sp.ones_like(diag[pores]))
            temp_data = sp.copy(data)
            temp_data[sp.in1d(row,diag[pores])] = 0
            non_Dir_diag = diag[-sp.in1d(diag,diag[pores])]
        except:
            temp_data = sp.copy(data)
            non_Dir_diag = diag
            
        S_temp = sp.zeros(A_dim)
        for i in sp.r_[0:len(row)]:
            S_temp[row[i]] = S_temp[row[i]] - temp_data[i]
        data = sp.append(data,S_temp[non_Dir_diag])
        row = sp.append(row,non_Dir_diag)
        col = sp.append(col,non_Dir_diag)
        
        self._Coeff_dimension = A_dim
        a = sprs.coo.coo_matrix((data,(row,col)),(A_dim,A_dim))
        A = a.tocsr()
        return(A)


    def _build_RHS_matrix(self):        
        r'''
        '''
        A_dim = self._Coeff_dimension
        B = sp.zeros([A_dim,1])
        try:
            Dir_pores = self.pores('Dirichlet')
            Dir_pores_vals = self['pore.bcval_Dirichlet'][Dir_pores]
            B[Dir_pores] = sp.reshape(Dir_pores_vals,[len(Dir_pores),1])
        except: pass
        try:
            individual_Neu_pores = self.pores('Neumann')
            individual_Neu_pores_vals = self['pore.bcval_Neumann'][individual_Neu_pores]
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
        
    def _do_one_inner_iteration(self):
        r'''
        '''
        self._logger.info("Creating Coefficient matrix for the algorithm")
        A = self._build_coefficient_matrix()
        self._logger.info("Creating RHS matrix for the algorithm")
        B = self._build_RHS_matrix()
        self._logger.info("Solving AX = B for the sparse matrices")
        X = sprslin.spsolve(A,B)
        self._Neumann_super_X = X[-sp.in1d(sp.r_[0:len(X)],sp.r_[0:self.num_pores()])]
        self._result = X[sp.r_[0:self.num_pores()]]        
        self._logger.info('Writing result to '+self.__class__.__name__+'[\''+self._conductance+'\']')
        self[self._quantity] = self._result
        
        
        