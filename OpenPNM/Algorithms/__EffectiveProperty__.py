"""
module __EffectiveProperty__: Base class to estimate transport properties
===============================================================================

"""
import scipy as sp
import scipy.signal as spsg
import scipy.spatial as sptl

from .__GenericAlgorithm__ import GenericAlgorithm

class EffectiveProperty(GenericAlgorithm):
    r'''
    '''
    def __init__(self,**kwargs):
        r'''
        '''
        super(EffectiveProperty,self).__init__(**kwargs)
        self._logger.info("Construct Algorithm")
        
    def run(self,algorithm,clean=False):
        r'''
        '''
        self._alg = algorithm
        self._fluid = algorithm._fluid
        self._conductance = algorithm._conductance
        self._quantity =  algorithm._quantity
        self._clean = clean
        if self._clean | ('pore.Dirichlet' not in algorithm.labels()):
            D = self._calc_eff_prop_tensor(fluid=self._fluid,alg=algorithm)
        else:
            #code that calls _execute for the algorithms preset boundaries.
            #Determine boundary conditions by analyzing algorithm object
            self._Ps = self._alg.pores(labels='pore.Dirichlet')
            self._BCs = sp.unique(self._alg['pore.bcval_Dirichlet'][self._Ps])
            if sp.shape(self._BCs)[0] != 2:
                raise Exception('The supplied algorithm did not have appropriate BCs')
            self._inlets = sp.where(self._alg['pore.bcval_Dirichlet']==sp.amax(self._BCs))[0]
            self._outlets = sp.where(self._alg['pore.bcval_Dirichlet']==sp.amin(self._BCs))[0]
            D = self._execute()
                
        return D
            
                
        
        
    def _execute(self):
        

        #Analyze input and output pores
        #Check for coplanarity
#        if self._net.iscoplanar(inlets) == False:
#            raise Exception('The inlet pores do not define a plane')
#        if self._net.iscoplanar(outlets) == False:
#            raise Exception('The outlet pores do not define a plane')
        #Ensure pores are on a face of domain (only 1 non-self neighbor each)
#        PnI = self._net.find_neighbor_pores(pores=inlets,mode='not_intersection',excl_self=True)
#        if sp.shape(PnI) != sp.shape(inlets):
#            raise Exception('The inlet pores have too many neighbors')
#        PnO = self._net.find_neighbor_pores(pores=outlets,mode='not_intersection',excl_self=True)
#        if sp.shape(PnO) != sp.shape(outlets):
#            raise Exception('The outlet pores have too many neighbors')
        inlets = self._inlets
        outlets = self._outlets
        
        
        #Fetch area and length of domain
        A = self._net.domain_area(face=self._inlets)
        L = self._net.domain_length(face_1=self._inlets,face_2=self._outlets)
        x = self._alg[self._quantity]
        #Find flow through inlet face
        Pin = []
        Pn = []
        for pore in inlets:
            pore_value = x[pore]
            neighbors = self._net.find_neighbor_pores(pore, excl_self = True)
            for neighbor in neighbors:
                neighbor_value = x[neighbor]
                if(sp.absolute(neighbor_value - pore_value) > .000001):
                    Pin.append(pore)
                    Pn.append(neighbor)
        
        Ts = self._net.find_connecting_throat(Pin,Pn)
        g = self._fluid[self._conductance][Ts]
        s = self._fluid['throat.occupancy'][Ts]
        xin = x[Pin]
        xout = x[Pn]
        flow = g*(xin - xout)
        D = sp.sum(flow)*L/A/sp.absolute(self._BCs[0]-self._BCs[1])
        return D
        
    def _calc_eff_prop_tensor(self):
        self._logger.error('Cannot compute tensor effective property. Not written yet')
#        
#        
#        ,                            
#                       fluid,
#                       alg,
#                       d_term,
#                       x_term,
#                       conductance,
#                       occupancy,
#                       direction,
#                       **params):
#                
#        network =self._net
#        ftype1 = []
#        ftype2 = []
#        effective_prop = []  
#        result = {}
#        try: fluid = self.find_object_by_name(fluid) 
#        except: pass #Accept object
#        if type(direction)==str and direction=='': 
#            ftype1 = ['front','right','top']
#            ftype2 = ['back','left','bottom']            
#        elif type(direction)==str: direction = sp.array(direction,ndmin=1)              
#        if type(direction)==sp.ndarray:
#            ftype1 = []
#            ftype2 = []
#            for d in direction:
#                if d=='x' or d=='X' or d=='front' or d=='back': 
#                    ftype1.append('front')
#                    ftype2.append('back')
#                elif d=='y' or d=='Y'or d=='left' or d=='right': 
#                    ftype1.append('left')
#                    ftype2.append('right')
#                elif d=='z' or d=='Z'or d=='top' or d=='bottom': 
#                    ftype1.append('top')
#                    ftype2.append('bottom') 
#                else: self._logger.error('wrong input for direction!')
#        
#        if 'pore.Dirichlet' in self:
#            self._dir = self.get_pore_info(label='Dirichlet')
#            del self['pore.Dirichlet']
#        if 'pore.BCval' in self:
#            self._BCval_temp = self.get_pore_data(prop='BCval')
#            del self['pore.BCval']
#            try:
#                self._BCtypes_temp = sp.copy(self._BCtypes)
#                delattr (self,'_BCtypes')
#                self._BCvalues_temp = sp.copy(self._BCvalues)
#                delattr(self,'_BCvalues')
#            except: pass
#        try: self._X_temp = self.get_pore_data(prop=self._X_name)  
#        except: pass          
#        tensor = sp.zeros([3,3])
#        for i in sp.r_[0:len(ftype1)]:
#            face1 = ftype1[i] 
#            face2 = ftype2[i]
#            if face1=='front' or face1=='back': direct = 'X'
#            elif face1=='left' or face1=='right': direct = 'Y'
#            elif face1=='top' or face1=='bottom': direct = 'Z'
#            if 'boundary' in self._net._pore_info:
#                face1_pores = network.get_pore_indices(labels=[face1,'boundary'],mode='intersection')
#                face2_pores = network.get_pore_indices(labels=[face2,'boundary'],mode='intersection')
#            else:    
#                face1_pores = network.get_pore_indices(face1)
#                face2_pores = network.get_pore_indices(face2)            
#            ## Assign Dirichlet boundary conditions
#            ## BC1
#            BC1_pores = face1_pores  
#            self.set_pore_info(label='Dirichlet',locations=BC1_pores,mode='overwrite')
#            BC1_values = 0.8
#            self.set_pore_data(prop='BCval',data=BC1_values,locations=BC1_pores)
#            ## BC2
#            BC2_pores = face2_pores
#            self.set_pore_info(label='Dirichlet',locations=BC2_pores)
#            BC2_values = 0.4
#            self.set_pore_data(prop='BCval',data=BC2_values,locations=BC2_pores)        
#            self.run(active_fluid=fluid,
#                           x_term=x_term,
#                           conductance=conductance,
#                           occupancy=occupancy) 
#            x = self.get_pore_data(prop=x_term)
#            if alg=='Fickian':
#                X1 = sp.log(1-x[face1_pores])
#                X2 = sp.log(1-x[face2_pores])
#            elif alg=='Stokes':
#                X1 = x[face1_pores]
#                X2 = x[face2_pores]
#            delta_X = sp.absolute(sp.average(X2)-sp.average(X1)) 
#            d_force =sp.average(fluid.get_pore_data(prop=d_term))
#
#            if  face1=='top' or face1=='bottom': 
#                L = self._net.domain_size('height')
#                A = self._net.domain_size('top')
#            elif  face1=='left' or face1=='right':
#                L = self._net.domain_size('depth')
#                A = self._net.domain_size('left')
#            elif  face1=='front' or face1=='back':
#                L = self._net.domain_size('width')
#                A = self._net.domain_size('front')
#            fn = network.find_neighbor_pores(face1_pores,excl_self=True)
#            fn = fn[sp.in1d(fn,network.get_pore_indices('internal'))]
#            ft = network.find_connecting_throat(face1_pores,fn)
#            if alg=='Fickian': X_temp = sp.log(1-x[fn])
#            elif alg=='Stokes':
#                X_temp = x[fn]
#                d_force = 1/d_force
#            cond = self._conductance
#            N = sp.sum(cond[ft]*sp.absolute(X1-X_temp))
#            eff = N*L/(d_force*A*delta_X)
#            effective_prop.append(eff)
#            del self._pore_info['Dirichlet']
#            del self._pore_data['BCval']
#            delattr (self,'_BCtypes')
#            delattr(self,'_BCvalues')            
#            result[ftype1[i]+'/'+ftype2[i]+'('+direct+')'] = sp.array(effective_prop[i],ndmin=1)
#            if ftype1[i]=='top' or ftype1[i]=='bottom': tensor[2,2] = effective_prop[i]
#            elif ftype1[i]=='right' or ftype1[i]=='left': tensor[1,1] = effective_prop[i]
#            elif ftype1[i]=='front' or ftype1[i]=='back': tensor[0,0] = effective_prop[i]
#        
#        try:
#            self.set_pore_data(prop=self._X_name,data=self._X_temp)
#            delattr (self,'_X_temp')
#        except : del self._pore_data[self._X_name]
#        try:
#            self.set_pore_info(label='Dirichlet',locations=self._dir,mode='overwrite')
#            delattr (self,'_dir')
#            self.set_pore_data(prop='BCval',data=self._BCval_temp)
#            delattr (self,'_BCval_temp')
#            self._BCtypes = self._BCtypes_temp
#            delattr (self,'_BCtypes_temp')
#            self._BCvalues = self._BCvalues_temp
#            delattr (self,'_BCvalues_temp')
#        except: pass        
#        if len(ftype1)<3: return result
#        elif len(ftype1)==3 : return tensor
        
        
        
        
        
        
        
        
        
        
        
        
        
        