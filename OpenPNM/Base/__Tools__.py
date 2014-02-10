"""
module __Tools__: Base class to construct pore network tools
==================================================================

.. warning:: The classes of this module should be loaded through the 'Base.__init__.py' file.

"""

import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if sys.path[1] != parent_dir:
    sys.path.insert(1, parent_dir)
import OpenPNM
import scipy as sp
from OpenPNM.Base import Utilities

class Tools(Utilities):
    r'''
    This class contains tools to read and write data in OpenPNM objects

    '''

    def __init__(self, **kwargs):
        r'''
        Initialize
        '''
        super(Tools,self).__init__(**kwargs)
        self._logger.debug("Construct Network container")
        #Initialize network properties dictionaries
        self._pore_data = {}
        self._pore_info = {}
        self._throat_data = {}
        self._throat_info = {}
        #Initialize fluid, physics, and geometry tracking lists
        self._fluids = []
        self._geometry = []
        self._logger.info("Construction of Network container complete")
        
    #--------------------------------------------------------------------------
    '''Setter and Getter Methods'''
    #--------------------------------------------------------------------------
    def set_data(self,element='',subdomain='',phase='',prop='',data='',indices=''):
        r'''
        '''
        try: subdomain = subdomain.name #allow passing of geometry objects
        except: pass #Otherwise, accept string
        try: phase = self.find_object_by_name(phase) #allow passing of fluid name by string
        except: pass #Accept object
        if phase and not subdomain: #Set fluid property
            try: getattr(phase,'_'+element+'_data')[prop]
            except: getattr(phase,'_'+element+'_data')[prop] = sp.zeros((getattr(phase,'get_num_'+element+'s')(),))
            if indices!='': getattr(phase,'_'+element+'_data')[prop][indices] = sp.array(data,ndmin=1)
            else: getattr(phase,'_'+element+'_data')[prop] = sp.array(data,ndmin=1)
            
        elif subdomain and not phase: #Set geometry property
            ind = getattr(self,'get_'+element+'_info')(subdomain)
            try: getattr(self,'_'+element+'_data')[prop] #Test existance of prop
            except: getattr(self,'_'+element+'_data')[prop] = sp.zeros((getattr(self,'get_num_'+element+'s')(),))*sp.nan
            if indices!='':
                if (sp.in1d(getattr(self,'get_'+element+'_indices')()[indices],\
                getattr(self,'get_'+element+'_indices')(subdomain))).all():
                    ind_temp = sp.zeros((getattr(self,'get_num_'+element+'s')(),),dtype=bool)
                    ind_temp[indices] = True
                    ind = ind_temp
                else: self._logger.error('Some/all of these indices do not belong to this subdomain!')
            r'''
            TODO:Shouldn't there be else here, for indices==''?  This doesn't seem to write scalar data either.
            '''
            if sp.sum(ind) == sp.shape(data)[0] or sp.shape(data)[0]==1:
                getattr(self,'_'+element+'_data')[prop][ind] = sp.array(data,ndmin=1)
            else: print('data is the wrong size!')
                
        elif phase and subdomain: #Set pore/throat scale physics property
            ind = getattr(self,'get_'+element+'_info')(subdomain)
            try: getattr(phase,'_'+element+'_data')[prop]
            except: getattr(phase,'_'+element+'_data')[prop] = sp.zeros((getattr(phase,'get_num_'+element+'s')(),))
            if indices!='':
                if (sp.in1d(getattr(self,'get_'+element+'_indices')()[indices],\
                getattr(self,'get_'+element+'_indices')(subdomain))).all():
                    ind_temp = sp.zeros((getattr(phase,'get_num_'+element+'s')(),),dtype=bool)
                    ind_temp[indices] = True
                    ind = ind_temp
                else: phase._logger.error('Some/all of these indices do not belong to this subdomain!')
            if sp.sum(ind) == sp.shape(data)[0] or sp.shape(data)[0]==1:
                getattr(phase,'_'+element+'_data')[prop][ind] = sp.array(data,ndmin=1)
            else: print('data is the wrong size!')
            
        elif not (phase or subdomain):  #Set topology property
            try: getattr(self,'_'+element+'_data')[prop]
            except: getattr(self,'_'+element+'_data')[prop] = sp.zeros_like(data)           
            if indices!='': getattr(self,'_'+element+'_data')[prop][indices] = sp.array(data,ndmin=1)
            else: getattr(self,'_'+element+'_data')[prop] = sp.array(data,ndmin=1)

    def get_data(self,element='',subdomain='',phase='',prop='',indices=''):
        r'''
        '''      
        try: subdomain = subdomain.name #allow passing of geometry objects
        except: pass #Otherwise, accept string
        try: phase = self.find_object_by_name(phase) #allow passing of fluid name by string
        except: pass #Accept object
        if phase and not subdomain:
            try: 
                getattr(phase,'_'+element+'_data')[prop]
                if indices!='':  return getattr(phase,'_'+element+'_data')[prop][indices]
                else: return getattr(phase,'_'+element+'_data')[prop] #Get fluid prop
            except: phase._logger.error(phase.name+' does not have the requested '+element+' property: '+prop)           
        elif subdomain and not phase: #Get geometry property
            ind = getattr(self,'get_'+element+'_info')(subdomain)            
            try: 
                getattr(self,'_'+element+'_data')[prop]                
                if indices!='':
                    if (sp.in1d(getattr(self,'get_'+element+'_indices')()[indices],\
                    getattr(self,'get_'+element+'_indices')(subdomain))).all():
                        ind_temp = sp.zeros((getattr(self,'get_num_'+element+'s')(),),dtype=bool)
                        ind_temp[indices] = True
                        ind = ind_temp
                    else: self._logger.error('Some/all of these indices do not belong to this subdomain!')
                return getattr(self,'_'+element+'_data')[prop][ind]
            except: self._logger.error(subdomain+' does not have the requested '+element+' property: '+prop)            
        elif phase and subdomain: #Get physics property
            ind = getattr(self,'get_'+element+'_info')(subdomain)            
            try: 
                getattr(phase,'_'+element+'_data')[prop]
                if indices!='':
                    if (sp.in1d(getattr(self,'get_'+element+'_indices')()[indices],\
                    getattr(self,'get_'+element+'_indices')(subdomain))).all():
                        ind_temp = sp.zeros((getattr(phase,'get_num_'+element+'s')(),),dtype=bool)
                        ind_temp[indices] = True
                        ind = ind_temp
                    else: phase._logger.error('Some/all of these indices do not belong to this subdomain!')                   
                return getattr(phase,'_'+element+'_data')[prop][ind]
            except: phase._logger.error(phase.name+'/'+subdomain+' does not have the requested '+element+' property: '+prop) 
        elif not (phase or subdomain): #Get topology property  
            try: 
                getattr(self,'_'+element+'_data')[prop]
                if indices!='':  return getattr(self,'_'+element+'_data')[prop][indices]
                else: return getattr(self,'_'+element+'_data')[prop] #Get fluid prop
            except: self._logger.error('Object does not have the requested '+element+' property: '+prop)      
 
    def set_pore_data(self,subdomain='',phase='',prop='',data='',indices=''):
        r'''
        Writes data to fluid or network objects according to input arguments.
        
        Parameters
        ----------
        prop : string
            Name of property to write
        subdomain : Open
        phase : OpenPNM Fluid object or fluid name string, optional
            Fluid to which data is written.  If omitted data is written to network object.
        data : array_like
            Data values to write to object
            
        See Also
        --------
        set_throat_data, set_pore_info, set_throat_info
        
        Notes
        -----
        This is wrapper method that calls set_data, which is generic for pores and throats
            
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> print(pn.name)
        test_network
        >>> pn.set_pore_data(prop='test',data=1)
        '''
        self.set_data(element='pore',subdomain=subdomain,phase=phase,prop=prop,data=data,indices=indices)
        
    def get_pore_data(self,subdomain='',phase='',prop='',indices=''):
        r'''
        Retrieves data from fluid or network objects according to input arguments.
        
        Parameters
        ----------
        prop : string
            Name of property to retrieve.  Requesting property 'all' prints a list of existing properties.
        fluid : string, optional
            Name of fluid from which to retrieve data.  If omitted data is retrieved from network object.
        
        Returns
        -------
        array_like
            An ndarray containing the requested property data from the specified object
            
        See Also
        --------
        get_throat_data, get_pore_info, get_throat_info

        Notes
        -----
        This is a wrapper method that calls get_data, which is generic for pores and throats
            
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> print(pn.name)
        test_network
        >>> pn.set_pore_data(prop='test',data=1)
        >>> pn.get_pore_data(prop='test')
        array([1])
        '''
        return self.get_data(element='pore',subdomain=subdomain,phase=phase,prop=prop,indices=indices)

    def set_throat_data(self,subdomain='',phase='',prop='',data='',indices=''):
        r'''
        Writes data to fluid or network objects according to input arguments.  
        Network topology data and pore/throat geometry data is stored on the network object.
        Fluid properties and physics properties are stored on the corresponding fluid object.
        
        Parameters
        ----------
        prop : string
            Name of property to write
        fluid : OpenPNM fluid object or fluid name string, optional
            Fluid to which data is written.  If omitted data is written to network object.
        data : array_like
            Data values to write to object
            
        See Also
        --------
        set_pore_data, set_pore_info, set_throat_info
            
        Notes
        -----
        This is wrapper method that calls set_data, which is generic for pores and throats
        
        Examples
        --------
        See set_pore_data
        '''
        self.set_data(element='throat',subdomain=subdomain,phase=phase,prop=prop,data=data,indices=indices)         

    def get_throat_data(self,subdomain='',phase='',prop='',indices=''):
        r'''
        Retrieves data from fluid or network objects according to input arguments.
        
        Parameters
        ----------
        prop : string
            Name of property to retrieve.  Requesting property 'all' prints a list of existing properties.
        fluid : string, optional
            Name of fluid from which to retrieve data.  If omitted data is retrieved from network object.
        
        Returns
        -------
        array_like
            An ndarray containing the requested property data from the specified object
            
        See Also
        --------
        get_pore_data, get_pore_info, get_throat_info

        Notes
        -----
        This is a wrapper method that calls get_data, which is generic for pores and throats
        
        Examples
        --------
        See get_pore_data
        '''
        return self.get_data(element='throat',subdomain=subdomain,phase=phase,prop=prop,indices=indices)     

    def set_info(self,element='',prop='',locations='',is_indices=False,mode='merge'):
        r'''
        This is the actual info setter method, but it should not be called directly.  
        Wrapper methods have been created.  Use set_pore_info and get_pore_info.
        
        See Also
        --------
        set_pore_info, set_throat_info
        '''
        if mode=='overwrite':
            getattr(self,'_'+element+'_info')[prop] = sp.zeros((getattr(self,'get_num_'+element+'s')(),),dtype=bool)
        if is_indices:
            try: 
                getattr(self,'_'+element+'_info')[prop]
            except: getattr(self,'_'+element+'_info')[prop] = sp.zeros((getattr(self,'get_num_'+element+'s')(),),dtype=bool)
            getattr(self,'_'+element+'_info')[prop][locations] = True
        else:
            getattr(self,'_'+element+'_info')[prop] = sp.array(locations,dtype=bool,ndmin=1)

    def get_info(self,element='',prop='',return_indices=False):
        r'''
        This is the actual info getter method, but it should not be called directly.  
        Wrapper methods have been created.  Use get_pore_info and get_throat_info
        
        See Also
        --------
        get_pore_info, get_throat_info
        
        TODO: It might be a good idea to have this return the T/F result of a list of indices given the prop type?
        '''
        if return_indices:
            return sp.where(getattr(self,'_'+element+'_info')[prop]==True)[0]
        else:
            return getattr(self,'_'+element+'_info')[prop]

    def set_pore_info(self,prop='',locations='',is_indices=False):
        r'''
        Parameters
        ----------
        prop : string
            The name of the pore label you wish to apply (e.g. 'top')
        locaitons : array_like
            An array containing the locations (pores) where the label should be applied.
            Can be either a boolean mask of Np length with True at label locations (default), 
            a list of indices where label should be applied. 
        is_indices : boolean
            This flag indicates whether locations are being sent as a boolean maks (default), 
            or a list of indices.
            
        See Also
        --------
        set_pore_data, set_throat_data, set_throat_info
            
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.set_pore_info(prop='test',locations=[0,1],is_indices=True) #Set using index notation
        >>> pn.get_pore_info(prop='test',return_indices=True) #Retrieve values as indices
        array([0, 1], dtype=int64)
        >>> loc = sp.zeros((pn.get_num_pores(),),dtype=bool)
        >>> loc[[0,1]] = True
        >>> pn.set_pore_info(prop='test2',locations=loc) #Set using boolean mask
        >>> pn.get_pore_info(prop='test',return_indices=True) #Retrieve values as indices
        array([0, 1], dtype=int64)
        '''
        self.set_info(element='pore',prop=prop,locations=locations,is_indices=is_indices)

    def get_pore_info(self,prop='',return_indices=False):
        r'''
        Retrieves locations where requested subdomain label is applies
        
        Parameters
        ----------
        prop : string
            The name of the label you wish to retrieve (e.g. 'top')
            
        return_indices : bool, optional
            This flag indicates that the returned result should be a list of indices
            
        Returns
        -------
        A boolean mask of length Np with True at all locations where subdomain label applies, 
        or a list of indices where label applies.
        
        See Also
        --------
        get_pore_data, get_throat_data, get_throat_info
        
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> result = pn.get_pore_info(prop='top',return_indices=True) #Retrieve values as indices
        >>> result[0:10]
        array([100, 101, 102, 103, 104, 105, 106, 107, 108, 109], dtype=int64)
        >>> result = pn.get_pore_info(prop='top') #Retrieve values as boolean mask
        >>> result[97:103]
        array([False, False, False,  True,  True,  True], dtype=bool)
        '''
        return self.get_info(element='pore',prop=prop,return_indices=return_indices)
        
    def set_throat_info(self,prop='',locations='',is_indices=False):
        r'''
        Parameters
        ----------
        prop : string
            The name of the pore label you wish to apply (e.g. 'top')
        locaitons : array_like
            An array containing the locations (pores) where the label should be applied.
            Can be either a boolean mask of Np length with True at label locations (default), 
            a list of indices where label should be applied. 
        is_indices : boolean
            This flag indicates whether locations are being sent as a boolean maks (default), 
            or a list of indices.
            
        See Also
        --------
        set_pore_data, set_throat_data, set_pore_info
            
        Examples
        --------
        See set_pore_info for usage
        '''
        self.set_info(element='throat',prop=prop,locations=locations,is_indices=is_indices)
        
    def get_throat_info(self,prop='',return_indices=False):
        r'''
        Retrieves locations where requested subdomain label is applies
        
        Parameters
        ----------
        prop : string
            The name of the label you wish to retrieve (e.g. 'top')
            
        return_indices : bool, optional
            This flag indicates that the returned result should be a list of indices
            
        Returns
        -------
        A boolean mask of length Np with True at all locations where subdomain label applies, 
        or a list of indices where label applies.
        
        See Also
        --------
        get_pore_data, get_throat_data, get_pore_info
        
        Examples
        --------
        See set_pore_info for usage
        '''
        return self.get_info(element='throat',prop=prop,return_indices=return_indices)
        
    def find_pore_labels(self,pnum):
        r'''
        Returns all the subdomain labels associated with the given pore
        '''
        labels = []
        for item in self._pore_info.keys():
            if self._pore_info[item][pnum]:
                labels.append(item)
        return labels
        
    def refresh_info(self):
        temp = sp.zeros_like(self.get_pore_data(prop='coords')[:,0],dtype=bool)
        self.set_pore_info(prop='all',locations=temp)
        temp = sp.zeros_like(self.get_throat_data(prop='connections')[:,0],dtype=bool)
        self.set_throat_info(prop='all',locations=temp)
    
    def find_throat_labels(self,tnum):
        r'''
        Returns all the subdomain labels associated with the given throat
        '''
        labels = []
        for item in self._throat_info.keys():
            if self._throat_info[item][tnum]:
                labels.append(item)
        return labels
        
    #--------------------------------------------------------------------------
    '''Object query methods'''
    #--------------------------------------------------------------------------
    def get_num_pores(self,subdomain=['all'],mode='union'):
        r'''
        Returns the number of pores of the specified subdomain

        Parameters
        ----------
        subdomain : list of strings, optional
            The pore subdomain labels that should be included in the count.  
            If not supplied, all pores are counted.
        mode : string, optional
            Specifies whether the count should be done as a union (default) or intersection.
            In union mode, all pores with ANY of the given labels are counted.
            In intersection mode, only pores with ALL the give labels are counted.
            
        Returns
        -------
        Np : int
            Number of pores with the specified subdomain label
            
        See Also
        --------
        get_num_throats
            
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.get_num_pores()
        125
        >>> pn.get_num_pores(subdomain=['top'])
        25
        >>> pn.get_num_pores(subdomain=['top','front'],mode='union') #'union' is default
        45
        >>> pn.get_num_pores(subdomain=['top','front'],mode='intersection')
        5
        
        '''
        #convert string to list, if necessary
        if type(subdomain) == str: subdomain = [subdomain]
        #Count number of pores of specified type
        if subdomain == ['all']: #return all pores
            return sp.shape(self.get_pore_info(prop='numbering'))[0]
        else:
            temp = self.get_pore_indices(subdomain=subdomain,mode=mode,indices=False)
            return sp.sum(temp) #return sum of Trues
            
    def get_num_throats(self,subdomain=['all'],mode='union'):
        r'''
        Return the number of throats of the specified subdomain

        Parameters
        ----------
        subdomain : list of strings, optional
            The throat subdomain labels that should be included in the count.  
            If not supplied, all throats are counted.
        mode : string, optional
            Specifies whether the count should be done as a union (default) or intersection.
            In union mode, all throats with ANY of the given labels are counted.
            In intersection mode, only throats with ALL the give labels are counted.

        Returns
        -------
        Nt : int
            Number of throats with the specified subdomain label
            
        See Also
        --------
        get_num_pores

        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.get_num_throats()
        300
        >>> pn.get_num_throats(subdomain=['top'])
        40
        >>> pn.get_num_throats(subdomain=['top','front'],mode='union') #'union' is default
        76
        >>> pn.get_num_throats(subdomain=['top','front'],mode='intersection')
        4
        
        '''
        #convert string to list, if necessary
        if type(subdomain) == str: subdomain = [subdomain]
        #Count number of pores of specified type
        if subdomain == ['all']: #return all pores
            return sp.shape(self.get_throat_info(prop='numbering'))[0]
        else:
            temp = self.get_throat_indices(subdomain=subdomain,mode=mode,indices=False)
            return sp.sum(temp) #return sum of Trues

    def get_pore_indices(self,subdomain=['all'],indices=True,mode='union'):
        r'''
        Returns pore locations where given subdomain labels exist.
        
        Parameters
        ----------
        subdomain : list of strings, optional
            The pore subdomain label(s) whose locations are requested.
            If omitted, all pore inidices are returned.
        indices : boolean, optional
            This flag specifies whether pore locations are returned a boolean mask of length Np,
            or as a list of indices (default).
        mode : string, optional
            Specifies whether the indices should be union (default) or intersection of desired labels.
            In union mode, all pores with ANY of the given labels are included.
            In intersection mode, only pores with ALL the give labels are included.
            This is ignored if no subdomain is given.
        
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.get_pore_indices(subdomain=['top','front'],mode='intersection')
        array([100, 105, 110, 115, 120], dtype=int64)
        '''
        if type(subdomain) == str: subdomain = [subdomain] #convert string to list, if necessary
        if subdomain == ['all']: #Return full index; easier than get_data(prop='nums')
            if indices:
                ind = sp.r_[0:self.get_num_pores()]
            else:
                ind = sp.ones((self.get_num_pores(),),dtype=bool)
        else:
            if mode == 'union':
                union = sp.zeros((self.get_num_pores(),),dtype=bool)
                for item in subdomain: #iterate over subdomain list and collect all indices
                    union = union + self.get_info(element='pore',prop=item)
                ind = union
            elif mode == 'intersection':
                intersect = sp.ones((self.get_num_pores(),),dtype=bool)
                for item in subdomain: #iterate over subdomain list and collect all indices
                    intersect = intersect*self.get_info(element='pore',prop=item)
                ind = intersect
            if indices: ind = sp.where(ind==True)[0]
        return ind

    def get_throat_indices(self,subdomain=['all'],indices=True,mode='union'):
        r'''
        Returns throat locations where given subdomain labels exist.
        
        Parameters
        ----------
        subdomain : list of strings, optional
            The throat subdomain label(s) whose locations are requested.
            If omitted, all throat inidices are returned.
        indices : boolean, optional
            This flag specifies whether throat locations are returned as a boolean mask of length Np,
            or as a list of indices (default).
        mode : string, optional
            Specifies whether the indices should be union (default) or intersection of desired labels.
            In union mode, all throats with ANY of the given labels are included.
            In intersection mode, only throats with ALL the give labels are included.
            This is ignored if no subdomain is given.
        
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> Tind = pn.get_throat_indices()
        >>> Tind[0:5]
        array([0, 1, 2, 3, 4])
        '''
        if type(subdomain) == str: subdomain = [subdomain] #convert string to list, if necessary
        if subdomain == ['all']: #Return full index; easier than get_data(prop='nums')
            if indices:
                ind = sp.r_[0:self.get_num_throats()]
            else:
                ind = sp.ones((self.get_num_throats(),),dtype=bool)
        else:
            if mode == 'union':
                union = sp.zeros((self.get_num_throats(),),dtype=bool)
                for item in subdomain: #iterate over subdomain list and collect all indices
                    union = union + self.get_info(element='throat',prop=item)
                ind = union
            elif mode == 'intersection':
                intersect = sp.ones((self.get_num_throats(),),dtype=bool)
                for item in subdomain: #iterate over subdomain list and collect all indices
                    intersect = intersect*self.get_info(element='throat',prop=item)
                ind = intersect
            if indices: ind = sp.where(ind==True)[0]
        return ind
        
    def find_object_by_name(self,name):
        r'''
        This is a short-cut method.  Given the string name of an 
        OpenPNM Fluid, Geometry, Physics, Algorithm, or Network object 
        this method will return that object
        
        Parameters
        ----------
        name : string
            Unique name of desired object
        
        Returns
        -------
        OpenPNM Object
            
        Notes
        -----
        If any objects are instantiated without a name (i.e. name = ''), then
        this method may start failing since the default name in many method calls
        is name = ''.
        
        '''
        for item in self._instances:
            if item.name == name:
                obj = item
        return obj

    def get_result(self,alg_obj,**kwargs):
        r'''
        This method invokes the update method on the given OpenPNM Algorithm object
        
        Parameters
        ----------
        alg_obj : OpenPNM Algorithm object
        
        Notes
        -----
        This method accepts keyword arguments which it passes on to algorithm object.
        For specific details refer to the `update` of the algorithm.
        '''
        alg_obj.update(**kwargs)        

if __name__ == '__main__':
    import doctest
    doctest.testmod(verbose=True)

