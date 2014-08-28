"""
===============================================================================
module __GenericLinearTransport__: Class for solving linear transport processes
===============================================================================

"""

import scipy as sp
import scipy.sparse as sprs
import scipy.sparse.linalg as sprslin
from .__GenericAlgorithm__ import GenericAlgorithm


class GenericLinearTransport(GenericAlgorithm):
    r"""
    This class provides essential methods for building and solving matrices 
    in a transport process.  It is inherited by FickianDiffusion, 
    FourierConduction, StokesFlow and OhmicConduction.

    """

    def __init__(self,**kwargs):
        r'''
        Initializing the class
        '''
        super(GenericLinearTransport,self).__init__(**kwargs)
        
    def setup(self,phase,conductance,quantity):
        r'''
        This setup provides the initial data for the solver
        '''
        self._phase = phase
        self._conductance = 'throat.'+conductance.split('.')[-1]
        self._quantity = 'pore.'+quantity.split('.')[-1]
        
        #Check health of conductance vector
        if self._phase.check_data_health(props=self._conductance,quiet=True):
            #If no nans, check for 0's
            ind = sp.nonzero(phase[self._conductance])[0]
            gmin = sp.amin(self._phase[self._conductance][ind])
            ind = sp.where(phase[self._conductance]==0)[0]
            self['throat.conductance'] = self._phase[self._conductance]
            #To prevent singular matrix
            self['throat.conductance'][ind] = gmin/1000000
        else:
            raise Exception('The provided throat conductance has problems')
            
    def update(self):
        r'''
        Send results of simulation out the the appropriate locations.
        
        This is a basic version of the update that simply sends out the main
        result (quantity). More elaborate updates should be subclassed.
        '''        
        self._phase[self._quantity] = self[self._quantity]
        self._logger.debug('Results of '+self.name+' algorithm have been added to '+self._phase.name)

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
            
        #If mode is 'remove', also bypass checks
        if mode == 'remove':
            if pores == []:
                self._logger.debug('Cannot remove BC values unless pores are specified')
            else:
                self._logger.debug('Removing '+bctype+' from specified pores')
                self['pore.bcval_'+bctype][pores] = sp.nan
                self['pore.'+bctype][pores] = False
            return

        #Validate bctype
        if bctype.split('.')[-1] not in BC_default:
            raise Exception('Unrecognized bctype') 
            
        #Validate pores
        if pores == []:
            raise Exception('Pores must be specified')
        else:
            pores = sp.array(pores,ndmin=1)
                    
        #Validate bcvalue
        if bcvalue == []:
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

    def _build_coefficient_matrix(self):
        r'''
        This builds the sparse coefficient matrix for the linear solver.
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
        g = self['throat.conductance']
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
        #Convert the lists to the sparse matrix
        self._Coeff_dimension = A_dim
        a = sprs.coo.coo_matrix((data,(row,col)),(A_dim,A_dim))
        A = a.tocsr()
        return(A)


    def _build_RHS_matrix(self):        
        r'''
        This builds the right-hand-side matrix for the linear solver.
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

    def rate(self,pores='',mode='group'):
        r'''
        Send a list of pores and receive the net rate
        of material moving into them.
        
        Parameters
        ----------
        pores : array_like
            The pores where the net rate will be calculated
        mode : string, optional
            Controls how to return the rate.  Options are:
            
            - 'group'(default): It returns the cumulative rate moving into them
            - 'single': It calculates the rate for each pore individually.       
        
        '''
        pores = sp.array(pores,ndmin=1)
        R = []
        if mode=='group':   iteration = 1
        elif mode=='single':    iteration = sp.shape(pores)[0]
        for i in sp.r_[0:iteration]:
            if mode=='group':   P = pores
            elif mode=='single':    P = pores[i]
            throats = self._net.find_neighbor_throats(P,flatten=True,mode='not_intersection')
            p1 = self._net.find_connected_pores(throats)[:,0]
            p2 = self._net.find_connected_pores(throats)[:,1]
            pores1 = sp.copy(p1)
            pores2 = sp.copy(p2)
            #Changes to pores1 and pores2 to make them as the internal and external pores
            pores1[-sp.in1d(p1,P)] = p2[-sp.in1d(p1,P)]
            pores2[-sp.in1d(p1,P)] = p1[-sp.in1d(p1,P)]
            X1 = self[self._quantity][pores1]
            X2 = self[self._quantity][pores2]
            g = self['throat.conductance'][throats]
            R.append(sp.sum(sp.multiply(g,(X2-X1))))
        return(sp.array(R,ndmin=1))
        
    def _do_one_inner_iteration(self):
        r'''
        This method collects the A and B matrices, solves AX = B and returns the result to the corresponding algorithm.
        '''
        self._logger.info("Creating Coefficient matrix for the algorithm")
        A = self._build_coefficient_matrix()
        self._logger.info("Creating RHS matrix for the algorithm")
        B = self._build_RHS_matrix()
        self._logger.info("Solving AX = B for the sparse matrices")
        X = sprslin.spsolve(A,B)
        self._Neumann_super_X = X[-sp.in1d(sp.r_[0:len(X)],sp.r_[0:self.num_pores()])]
        #Removing the additional super pore variables from the results
        self[self._quantity] = X[sp.r_[0:self.num_pores()]]        
        self._logger.info('Writing the results to '+'[\''+self._quantity+'\'] in the '+self.name+' algorithm.')
            
    def _calc_eff_prop(self,check_health=False):
        r'''
        This returns the main parameters for calculating the effective property in a linear transport equation. 
        It also checks for the proper boundary conditions, inlets and outlets.
        
        Parameters
        ----------
        check_health : boolean(optional)
            It analyzes the inlet and outlet pores to check their spatial positions
        '''
        try:
            self[self._quantity]
        except:
            raise Exception('The algorithm has not been run yet. Cannot calculate effective property.')               
        #Determine boundary conditions by analyzing algorithm object
        Ps = self.pores(labels='pore.Dirichlet')
        BCs = sp.unique(self['pore.bcval_Dirichlet'][Ps])
        if sp.shape(BCs)[0] != 2:
            raise Exception('The supplied algorithm did not have appropriate BCs')
        inlets = sp.where(self['pore.bcval_Dirichlet']==sp.amax(BCs))[0]
        outlets = sp.where(self['pore.bcval_Dirichlet']==sp.amin(BCs))[0]        

        #Analyze input and output pores
        if check_health:
            #Check for coplanarity
            if self._net.iscoplanar(inlets) == False:
                raise Exception('The inlet pores do not define a plane. Effective property will be approximation')
            if self._net.iscoplanar(outlets) == False:
                raise Exception('The outlet pores do not define a plane. Effective property will be approximation')
            #Ensure pores are on a face of domain (only 1 non-self neighbor each)
            PnI = self._net.find_neighbor_pores(pores=inlets,mode='not_intersection',excl_self=True)
            if sp.shape(PnI) != sp.shape(inlets):
                self._logger.warning('The inlet pores have too many neighbors. Internal pores appear to be selected.')
            PnO = self._net.find_neighbor_pores(pores=outlets,mode='not_intersection',excl_self=True)
            if sp.shape(PnO) != sp.shape(outlets):
                self._logger.warning('The outlet pores have too many neighbors. Internal pores appear to be selected.')        
        
        #Fetch area and length of domain
        A = self._net.domain_area(face=inlets)
        L = self._net.domain_length(face_1=inlets,face_2=outlets)
        flow = self.rate(pores=inlets)
        D = sp.sum(flow)*L/A/(BCs[0]-BCs[1])
        return D
        
        