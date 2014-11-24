# -*- coding: utf-8 -*-
"""
===============================================================================
module __GenericLinearTransport__: Class for solving linear transport processes
===============================================================================

"""
import scipy as sp
import scipy.sparse as sprs
import scipy.sparse.linalg as sprslin
from OpenPNM.Algorithms import GenericAlgorithm
import OpenPNM.Utilities.vertexops as vo
from OpenPNM.Base import logging
logger = logging.getLogger()

class GenericLinearTransport(GenericAlgorithm):
    r"""
    This class provides essential methods for building and solving matrices
    in a transport process.  It is inherited by FickianDiffusion,
    FourierConduction, StokesFlow and OhmicConduction.

    """

    def __init__(self,phase=None,**kwargs):
        r'''
        Initializing the class
        '''
        super(GenericLinearTransport,self).__init__(**kwargs)
        if phase is None:
            self._phase = OpenPNM.Phases.GenericPhase()
        else:
            self._phase = phase  # Register phase with self
            if sp.size(phase)!=1:   self._phases = phase
            else:   self._phases.append(phase)

    def setup(self,conductance,quantity,super_pore_conductance):
        r'''
        This setup provides the initial data for the solver
        '''
        if super_pore_conductance is None:  self.super_pore_conductance = []
        else:   self.super_pore_conductance =  super_pore_conductance 
        
        if  sp.size(self._phase)==1:
            self._conductance = 'throat.'+conductance.split('.')[-1]
            self._quantity = 'pore.'+self._phase.name+'_'+quantity.split('.')[-1]

            #Check health of conductance vector
            if self._phase.check_data_health(props=self._conductance,quiet=True):
                self['throat.conductance'] = self._phase[self._conductance]
            else:
                raise Exception('The provided throat conductance has problems')
        else:
            raise Exception('The linear transport solver accepts just one phase.')

        logger.info("Creating Coefficient matrix for the algorithm")
        self.A  = self._build_coefficient_matrix()
        logger.info("Creating RHS matrix for the algorithm")
        self.b = self._build_RHS_matrix()


    def run(self,**kwargs):
        r'''
        '''
        logger.info("Setup "+self.__class__.__name__)
        self.setup(**kwargs)
        self.solve()


    def return_results(self,pores=None,throats=None,**kwargs):
        r'''
        Send results of simulation out the the appropriate locations.

        This is a basic version of the update that simply sends out the main
        result (quantity). More elaborate updates should be subclassed.
        '''
        if pores is None:
            pores = self.Ps
        if throats is None:
            throats = self.Ts

        phase_quantity = self._quantity.replace(self._phase.name+'_',"")
        if phase_quantity not in self._phase.props():
            self._phase[phase_quantity] = sp.nan
        self._phase[phase_quantity][pores] = self[self._quantity][pores]

        dx = sp.squeeze(sp.diff(self[self._quantity][self._net.find_connected_pores(self.throats())],n=1,axis=1))
        g = self['throat.conductance']
        rate = sp.absolute(g*dx)
        if 'throat.rate' not in self._phase.props():
            self._phase['throat.rate'] = sp.nan
        self._phase['throat.rate'][throats] = rate[throats]
        logger.debug('Results of '+self.name+' algorithm have been added to '+self._phase.name)



    def _build_coefficient_matrix(self):
        r'''
        This builds the sparse coefficient matrix for the linear solver.
        '''
        # Filling coefficient matrix
        tpore1 = self._net['throat.conns'][:,0]
        tpore2 = self._net['throat.conns'][:,1]

        #Identify Dirichlet pores
        try:
            temp = self.pores(self._phase.name+'_Dirichlet',mode='difference')
        except:
            raise Exception('The linear transport solver needs at least one Dirichlet boundary condition for the phase which is attached to '+self.name)
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
            self.pores(self._phase.name+'_Neumann_group')
            self._extra_Neumann_size = len(getattr(self,'_pore_'+self._phase.name+'_Neumann_group_location'))
            self._group_Neumann_vals = sp.zeros(self._extra_Neumann_size)
            
            for N in sp.r_[0:self._extra_Neumann_size]:
                neu_tpore2 = getattr(self,'_pore_'+self._phase.name+'_Neumann_group_location')[N]
                self._group_Neumann_vals[N] = sp.unique(self['pore.'+self._phase.name+'_bcval_Neumann_group'][neu_tpore2])
                neighbor_throats = self._net.find_neighbor_throats(pores=neu_tpore2)
                try:   g_super = self.super_pore_conductance[N]
                except: 
                    g_super = 1e-3*min(data_main[neighbor_throats])
                    self.super_pore_conductance.append(g_super)
                row = sp.append(row,neu_tpore2)
                col = sp.append(col,len(neu_tpore2)*[A_dim+N])
                data = sp.append(data,len(neu_tpore2)*[g_super])
                row = sp.append(row,len(neu_tpore2)*[A_dim+N])
                col = sp.append(col,neu_tpore2)
                data = sp.append(data,len(neu_tpore2)*[g_super])
            A_dim = A_dim + self._extra_Neumann_size
        except:
            pass

        # Adding positions for diagonal
        diag = sp.r_[0:A_dim]
        try:
            pores = self.pores(self._phase.name+'_Dirichlet')
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
        if hasattr(self,'_k'):
            S_temp = S_temp + self._k
        data = sp.append(data,S_temp[non_Dir_diag])
        row = sp.append(row,non_Dir_diag)
        col = sp.append(col,non_Dir_diag)
        #Convert the lists to the sparse matrix
        self._Coeff_dimension = A_dim
        a = sprs.coo.coo_matrix((data,(row,col)),(A_dim,A_dim))
        A = a.tocsr()
        A.eliminate_zeros()        
        return(A)


    def _build_RHS_matrix(self):
        r'''
        This builds the right-hand-side matrix for the linear solver.
        '''
        A_dim = self._Coeff_dimension
        b = sp.zeros([A_dim,1])
        try:
            Dir_pores = self.pores(self._phase.name+'_Dirichlet')
            Dir_pores_vals = self['pore.'+self._phase.name+'_bcval_Dirichlet'][Dir_pores]
            b[Dir_pores] = sp.reshape(Dir_pores_vals,[len(Dir_pores),1])
        except: pass
        try:
            individual_Neu_pores = self.pores(self._phase.name+'_Neumann')
            individual_Neu_pores_vals = self['pore.'+self._phase.name+'_bcval_Neumann'][individual_Neu_pores]
            b[individual_Neu_pores] = sp.reshape(individual_Neu_pores_vals,[len(individual_Neu_pores),1])
        except: pass
        try:
            self.pores(self._phase.name+'_Neumann_group')
            pnum = self._net.num_pores()
            b[sp.r_[pnum:(pnum+len(self._group_Neumann_vals))]] = sp.reshape(self._group_Neumann_vals[sp.r_[0:len(self._group_Neumann_vals)]],[len(self._group_Neumann_vals),1])
        except: pass
       
        return(b)

    def rate(self,pores=None,network=None,conductance=None,X_value=None,mode='group'):
        r'''
        Send a list of pores and receive the net rate
        of material moving into them.

        Parameters
        ----------
        pores : array_like
            The pores where the net rate will be calculated
        network : OpenPNM Network Object
            The network object to which this algorithm will apply. 
            If no network is sent, the rate will apply to the network which is attached to the algorithm.        
        conductance : array_like
            The conductance which this algorithm will use to calculate the rate. 
            If no conductance is sent, the rate will use the 'throat.conductance' which is attached to the algorithm.         
        X_value : array_like
            The values of the quantity (temperature, mole_fraction, voltage, ...), which this algorithm will use to calculate the rate. 
            If no X_value is sent, the rate will look at the '_quantity', which is attached to the algorithm.        
        mode : string, optional
            Controls how to return the rate.  Options are:
            - 'group'(default): It returns the cumulative rate moving into them
            - 'single': It calculates the rate for each pore individually.

        '''
        if network is None: network = self._net
        if conductance is None: conductance = self['throat.conductance']
        if X_value is None: X_value = self[self._quantity]
        pores = sp.array(pores,ndmin=1)
        R = []
        if mode=='group':   
            t = network.find_neighbor_throats(pores,flatten=True,mode='not_intersection')
            throat_group_num = 1
        elif mode=='single':
            t = network.find_neighbor_throats(pores,flatten=False,mode='not_intersection')
            throat_group_num = sp.size(t)
        
        for i in sp.r_[0:throat_group_num]:
            if mode=='group':   
                throats = t
                P = pores
            elif mode=='single':
                throats = t[i]
                P = pores[i]
            p1 = network.find_connected_pores(throats)[:,0]
            p2 = network.find_connected_pores(throats)[:,1]
            pores1 = sp.copy(p1)
            pores2 = sp.copy(p2)
            #Changes to pores1 and pores2 to make them as the inner and outer pores
            pores1[-sp.in1d(p1,P)] = p2[-sp.in1d(p1,P)]
            pores2[-sp.in1d(p1,P)] = p1[-sp.in1d(p1,P)]
            X1 = X_value[pores1]
            X2 = X_value[pores2]
            g = conductance[throats]
            R.append(sp.sum(sp.multiply(g,(X2-X1))))
        return(sp.array(R,ndmin=1))

    def solve(self):
        r'''
        This method collects the A and b matrices, solves AX = b and returns the result to the corresponding algorithm.
        '''
        logger.info("Solving AX = b for the sparse matrices")
        X = sprslin.spsolve(self.A,self.b)
        self.X = X
        self._Neumann_super_X = X[-sp.in1d(sp.r_[0:len(X)],sp.r_[0:self.num_pores()])]
        #Removing the additional super pore variables from the results
        self[self._quantity] = X[sp.r_[0:self.num_pores()]]
        logger.info('Writing the results to '+'[\''+self._quantity+'\'] in the '+self.name+' algorithm.')

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
        Ps = self.pores('pore.'+self._phase.name+'_Dirichlet')
        BCs = sp.unique(self['pore.'+self._phase.name+'_bcval_Dirichlet'][Ps])
        if sp.shape(BCs)[0] != 2:
            raise Exception('The supplied algorithm did not have appropriate BCs')
        inlets = sp.where(self['pore.'+self._phase.name+'_bcval_Dirichlet']==sp.amax(BCs))[0]
        outlets = sp.where(self['pore.'+self._phase.name+'_bcval_Dirichlet']==sp.amin(BCs))[0]

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
                logger.warning('The inlet pores have too many neighbors. Internal pores appear to be selected.')
                pass
            PnO = self._net.find_neighbor_pores(pores=outlets,mode='not_intersection',excl_self=True)
            if sp.shape(PnO) != sp.shape(outlets):
                logger.warning('The outlet pores have too many neighbors. Internal pores appear to be selected.')
                pass

        #Fetch area and length of domain
        if "pore.vert_index" in self._net.props():
            A = vo.vertex_dimension(network = self._net,face1=inlets, parm='area')
            L = vo.vertex_dimension(network = self._net,face1=inlets,face2=outlets,parm='length')
        else:
            A = self._net.domain_area(face=inlets)
            L = self._net.domain_length(face_1=inlets,face_2=outlets)
        flow = self.rate(pores=inlets)
        D = sp.sum(flow)*L/A/(BCs[0]-BCs[1])
        return D

