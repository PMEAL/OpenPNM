# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 10:50:28 2013

@author: Mahmoudreza Aghighi

module __FickianDiffusion__: Fick's Law Diffusion
========================================================================

.. warning:: The classes of this module should be loaded through the 'Algorithms.__init__.py' file.

"""

import OpenPNM
import numpy as np
import scipy.sparse as sprs
import scipy.sparse.linalg as splin
from __GenericAlgorithm__ import GenericAlgorithm

class FickianDiffusion(GenericAlgorithm):
    r"""   
    
    FickianDiffusion - Class to run Fick's law mass transfer diffusion on constructed networks
    
                        It returns conecentration gradient inside the network.
                        An invasion algorithm should be used before running diffusion calculations.
                                   
                            
    Parameters
    ----------
    - Alg:
        Algorithm for Non wetting phase configuration
        'OP' , 'IP' and 'None' are possible options.        
    - Pressure for 'OP':
        Applied pressure on the network which causes some pores and throats be invaded by non-wetting phase.
        It can be a single value or a list.
        The class will consider uninvaded conduits based on this pressure, as possible voids for gas diffusion.
        Every pore with Pc_invade greater than Pressure, will be considered as open pores        
    - Psequence for 'IP':
        The desired sequence of invaded pores based on IP algorithm.
        The class will consider uninvaded conduits based on this sequence number, as possible voids for gas diffusion.
        It can be a single value or a list.
        Every pore with sequence number greater than Psequence, will be considered as open pores.
    - Pinvaded for 'None':
        It is a list which contains invaded pores and based on user's definition.
            
    -loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)    
    - Total_Conc: 
        Total gas concentration
    - Diff_Coefficient:
        Diffision coefficient for diffusion of A through stagnant B(e.g. Oxygen through Nitrogen and Water vapor)
                
    Examples
    --------

    All of the variables in this class have default values, but users can define them too.
    >>>
    
    To Do:
        - Instead of sending (Pressure in OP) or (Psequence in IP) as inputs, this class should only accept 
        a list of invaded voids.
        - loglevel and documentation is not complete yet.
        - For Applying Neumann boundary condition, a function is needed to consider mass conservation, sink or
        source terms and warn user about the applied values.
        
    """
    
    def __init__(self,net=OpenPNM.Network.GenericNetwork(),loglevel=10,Alg='None',Pressure=[0],Psequence=[0],Pinvaded=[],**kwargs):
        r"""
        Initializing the class
        """
        super(FickianDiffusion,self).__init__(net = net,**kwargs)
        self._logger.info("Create Fick's Diffusion Algorithm Object")
        self.Alg = Alg
        self.Pressure = Pressure
        self.Psequence = Psequence
        self.Pinvaded = Pinvaded

             
    def _setup_for_FicksDiffusion(self):
        r"""
        Main features: Applying Boundary Conditions & Creating Transport Conductivities 

        This function executes the essential mathods for building matrices in Linear solution 
        """
        self._logger.info("Create Diffusion Conductivities")        
        self._fill_Conductivity()
        self._logger.info("Applying Boundary Conditions")
        self._Boundary_Pores_Conditions()   

    def _Boundary_Pores_Conditions(self,FaceTypes=[4,1,4,4,1,4],FaceValues=[0,0.2,0,0,0.5,0],BList=[]):
        r"""        
        - Assigning Type and Value for Boundary Pores.
        - Types of Boundary Values in FaceTypes:
           internal=0, Dirichlet=1, flux=2, periodic=3, insulated=4 mass_rate=5
        - In FaceValues, each item represents the value of that face according
          to the boundary condition type. Except for mass_rate condition which 
          the value represents the total mass rate applied to the entire face,
          not just a single pore. 
          For instance: if we have Facetypes[2]=1 and FaceValues[2]=0.75,it means
          that for face 3, concentration=0.75
          but if we have FaceTypes[4]=2 and FaceValues[4]=0.8, it means that
          for face 5, flux=0.8
        - Instead of applying boundary conditions to entire face, by using 
          BLists, it is possible to apply BC to arbitrary pores. For instance:
          BLists = [(459,2,0.097),(2043,1,0.5)]
          It means that for Pore 459: flux=0.097 & for Pore 2043: Concentration=0.5
        """
        self._net.pore_properties['BCtype'] = np.zeros(self._net.get_num_pores())
        self._net.pore_properties['BCvalue'] = np.zeros(self._net.get_num_pores())
        for i in range(1,7):
            self._net.pore_properties['BCtype'][self._net.pore_properties['type']==i] = FaceTypes[i-1]
            self._net.pore_properties['BCvalue'][self._net.pore_properties['type']==i] = FaceValues[i-1]
        if BList:
            for item in BList:
                pore = item[0]
                self._net.pore_properties['BCtype'][pore] = item[1]
                self._net.pore_properties['BCvalue'][pore] = item[2]
        
    def _Inv_Conduit(self,P_val):      
        r"""
        creates a tlist containing which is useful to indicate the conduit is
        filled with invading fluid or not.
        
        """        
     
        if self.Alg=='None':
            self._logger.info("Invaded pores have already been determined")            
        else:
            if self.Alg=='OP':
                val_name = 'Pc_invaded'
                B_condition = self._net.pore_properties[val_name]<P_val
            elif self.Alg=='IP':
                val_name = 'IP_Pseq'
                list_B=[]
                for k in np.transpose(self._net.pore_properties[val_name]<P_val):
                    list_B.append(k[0])                     
                B_condition = np.array(list_B)

            self.Pinvaded = np.array(range(self._net.get_num_pores()))[B_condition]

        self._net.throat_properties['UninvadedConduits'] = np.zeros(self._net.get_num_throats())
        for i in range(0,self._net.get_num_throats()):
            neighborPs = self._net.get_connected_pores(i)
            if -np.in1d(neighborPs,self.Pinvaded).all():
                self._net.throat_properties['UninvadedConduits'] = 1
   
#
#            self._net.throat_properties['Pconduits'] = np.zeros(self._net.get_num_throats(),np.float)
#
#            for i in range(0,self._net.get_num_throats()):
#                neighborPs = self._net.get_connected_pores(i)
#                temp = self._net.pore_properties[val_name][neighborPs]
#                self._net.throat_properties['Pconduits'][i] = min(temp)

    
    def _do_outer_iteration_stage(self):
        r"""
        Running Fick's law for different pressures or sequences.
               
        """
        self._logger.info("Outer Iteration Stage ")
        self._setup_for_FicksDiffusion()
     
        if self.Alg=='None':
            Concentration = self._do_one_inner_iteration(P_val=0)
            Concentration = np.multiply(Concentration,-np.in1d(Concentration,self.Pinvaded))
            self._net.set_pore_property(name="Concentration",ndarray=Concentration)
        else:
            if self.Alg=='OP':
                Alg_var = self.Pressure
                val_name = 'Pc_invaded'
    
            elif self.Alg=='IP':
                Alg_var = self.Psequence
                val_name = 'IP_Pseq'
    
            for P_val in Alg_var:
                self._logger.info("Applying Pressure/Sequence = "+str(P_val))
                Concentration = self._do_one_inner_iteration(P_val)
                #Store result of Concentration in each step
                if P_val!=0:
                    Concentration = np.multiply(Concentration,self._net.pore_properties[val_name]>P_val)
          
                self._net.set_pore_property(name="Concentration_at_"+str(P_val),ndarray=Concentration)
                       

    def _do_one_inner_iteration(self,P_val):
        
        if hasattr(self._net,'dry_coefficient_matrix'):
            self.logger.info('Dry matrix has already been created')     
        else:            
           setattr(self._net,"dry_coefficient_matrix",self._Coefficient_Matrix(P_val=0))
        
        if ((P_val==0)and(len(self.Pinvaded)==0)):
            A = self._net.dry_coefficient_matrix
        else:
            A = self._Coefficient_Matrix(P_val)        
          
        if (self._net.pore_properties['BCtype']==2).any():
            B = np.zeros([self.A_dim,1])
            Dir_pores = np.array(range(self._net.get_num_pores()))[self._net.pore_properties['BCtype']==1]
            B[Dir_pores] = np.reshape(self._net.pore_properties['BCvalue'][Dir_pores],[len(Dir_pores),1])
            for item in len(self.extera_neu):
                F_type = np.unique(self._net.pore_properties['BCtype'][self._net.pore_properties['BCvalue']==self.extera_neu[item]])
                if F_type==1:
                    Area = self._net.divisions[1]*self._net.divisions[2]*self._net.lattice_spacing
                elif F_type==2:
                    Area = self._net.divisions[0]*self._net.divisions[1]*self._net.lattice_spacing
                elif F_type==3:
                    Area = self._net.divisions[0]*self._net.divisions[2]*self._net.lattice_spacing
                elif F_type==4:
                    Area = self._net.divisions[0]*self._net.divisions[2]*self._net.lattice_spacing
                elif F_type==5:
                    Area = self._net.divisions[0]*self._net.divisions[1]*self._net.lattice_spacing
                elif F_type==6:
                    Area = self._net.divisions[1]*self._net.divisions[2]*self._net.lattice_spacing                 
                B[self.A_dim-item,0] = - self.extera_neu[item]*Area
        else:
            A_dim = self._net.get_num_pores()             
            B = np.reshape(self._net.pore_properties['BCvalue'],[A_dim,1])


        for i in range(0,self._net.get_num_pores()):
            if self._net.pore_properties['BCvalue'][i]!=0:
                neighbors = self._net.get_neighbor_pores(i)
                if np.in1d(neighbors,self.Pinvaded).all():
                    B[i,0] = 0

        x = splin.spsolve(A,B)        
        return(x)
      
          
    def _fill_Conductivity(self):
        r"""

        """
        self._logger.info("Calculating Diffusion Conductivity for all of the Voids ")
        C = self._net.Total_Conc
        D = self._net.Diff_Coefficient
        self._net.set_throat_property(name="Cdiff",ndarray=np.zeros(self._net.get_num_throats()))
        tcond = C*D*np.abs(( self._net.throat_properties['diameter'])**2/self._net.throat_properties['length'])
        pcond = C*D*np.abs((self._net.pore_properties['diameter'])**2/(self._net.pore_properties['diameter']/2))
        for i in range(self._net.get_num_throats()):
            PNeighbors = self._net.get_connected_pores(i)
            if (self._net.pore_properties['type'][PNeighbors[0]]!=0):
                self._net.throat_properties['Cdiff'][i] = (1/tcond[i]+1/pcond[PNeighbors[1]]+1/pcond[PNeighbors[1]])**(-1)
            elif (self._net.pore_properties['type'][PNeighbors[1]]!=0):
                self._net.throat_properties['Cdiff'][i] = (1/pcond[PNeighbors[0]]+1/pcond[PNeighbors[0]]+1/tcond[i])**(-1)
            else:
                self._net.throat_properties['Cdiff'][i] = (1/pcond[PNeighbors[0]]+1/tcond[i]+1/pcond[PNeighbors[1]])**(-1)
     
            
    def _Coefficient_Matrix(self,P_val):
        
        # Determine open and closed conduits for mass transfer 
        if hasattr(self._net,'dry_coefficient_matrix'):
            
            if ((len(self.Pinvaded)>0)or(P_val>0)):
                self._Inv_Conduit(P_val)
                list_name = 'CdiffWet'
#                if P_val>0:
#                    if self.Alg=='OP':
#                        val_name = 'Pc_invaded'    
#                    elif self.Alg=='IP':
#                        val_name = 'IP_Pseq'
#                    self.Pinvaded = np.array(range(self._net.get_num_pores()))[self._net.pore_properties[val_name]<P_val]
#                    temp_list = np.multiply(self._net.throat_properties['Cdiff'],(self._net.throat_properties['Pconduits']>P_val))           
#                else:
                temp_list = np.multiply(self._net.throat_properties['Cdiff'],self._net.throat_properties['UninvadedConduits'])           
                self._net.set_throat_property(name=list_name,ndarray=temp_list)
        else:
            list_name = 'Cdiff'
        # Filling coefficient matrix

        nodes = self._net.throat_properties[list_name]>0 
        tpore1 = self._net.throat_properties['connections'][:,0]
        tpore2 = self._net.throat_properties['connections'][:,1]
        row = np.append(tpore1[nodes],tpore1[-nodes])
        col = np.append(tpore2[nodes],tpore2[-nodes])
        data= np.append(self._net.throat_properties[list_name][nodes],np.ones([len(tpore1[-nodes])])*1e-30)
        

        pnum = np.array(range(self._net.get_num_pores()))
        pdry = tpore2[nodes]       
        loc1 = np.in1d(tpore2,pdry[np.in1d(pdry,pnum[self._net.pore_properties['BCtype']!=1])])
        modified_tpore2_dry = tpore2[loc1]
        modified_tpore1_dry = tpore1[loc1]        
        row = np.append(row,modified_tpore2_dry)                
        col = np.append(col,modified_tpore1_dry)
        data = np.append(data,self._net.throat_properties[list_name][loc1])
        
        pwet = tpore2[-nodes]
        loc2 = np.in1d(tpore2,pwet[np.in1d(pwet,pnum[self._net.pore_properties['BCtype']!=1])])
        modified_tpore2_wet = tpore2[loc2]
        modified_tpore1_wet = tpore1[loc2]        
        row = np.append(row,modified_tpore2_wet)                
        col = np.append(col,modified_tpore1_wet)
        data = np.append(data,np.ones([len(modified_tpore2_wet)])*1e-30)
        
        if (self._net.pore_properties['BCtype']==2).any():
            self.extera_neu = np.unique(self._net.pore_properties['BCvalue'][self._net.pore_properties['BCtype']==2])
            self.A_dim = self._net.get_num_pores()+ len(self.extera_neu)
            extera_neu = self.extera_neu
            A_dim = self.A_dim

            for item in len(extera_neu):
                loc_neu = np.in1d(tpore2,pnum[self._net.pore_properties['BCvalue']==extera_neu[item]])
                neu_tpore2 = tpore2[loc_neu]                
                row = np.append(row,neu_tpore2) 
                col = np.append(col,np.ones([len(neu_tpore2)])*(A_dim-item))
                data = np.append(data,np.ones([len(neu_tpore2)])*1e-18)
                        
                row = np.append(row,np.ones([len(neu_tpore2)])*(A_dim-item)) 
                col = np.append(col,neu_tpore2)
                data = np.append(data,np.ones([len(neu_tpore2)])*1e-18)
            
        else:
            A_dim = self._net.get_num_pores()
        
        # Adding positions for diagonal
        
        row = np.append(row,range(0,A_dim))
        col = np.append(col,range(0,A_dim))
        data = np.append(data,np.zeros((A_dim,1)))

        a = sprs.coo.coo_matrix((data,(row,col)),(A_dim,A_dim))
        A = a.tocsr()        
       
        for i in range(0,A_dim):                    
            if self._net.pore_properties['BCtype'][i]==1:
                A[i,i] = 1                            
            else:
                A[i,i] = -np.sum(A[i,:][np.nonzero(A[i,:])])
                
        return(A)