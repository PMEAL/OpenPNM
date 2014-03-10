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
from OpenPNM.Utilities import Base

class Tools(Base):
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
        self._physics = []
        self._logger.info("Construction of Object container complete")
        
    #--------------------------------------------------------------------------
    '''Setter and Getter Methods'''
    #--------------------------------------------------------------------------
    def _set_data(self,element='',phase='',prop='',data='',locations=''):
        r'''
        Documentation for this method is being updated, we are sorry for the inconvenience.
        '''
        data = sp.array(data,ndmin=1)
        if type(locations)==list:
            try: locations = getattr(self,'get_'+element+'_indices')(locations)
            except: locations = sp.array(locations,ndmin=1)
        elif type(locations)==sp.ndarray:
            try: locations = getattr(self,'get_'+element+'_indices')(locations)
            except: pass            
        elif locations!='':
            try: locations = locations.name 
            except: pass
            if type(locations)==str: locations = getattr(self,'get_'+element+'_indices')([locations])

        if phase :
            try: phase = self.find_object_by_name(phase) 
            except: pass #Accept object
            try: 
                getattr(phase,'_'+element+'_data')[prop]
                temp_word = 'updated for '
            except: temp_word = 'added to '
            if sp.shape(data)[0]==1:
                if locations!='':
                    try: getattr(phase,'_'+element+'_data')[prop]
                    except: getattr(phase,'_'+element+'_data')[prop] = sp.zeros((getattr(phase,'num_'+element+'s')(),))*sp.nan
                    getattr(phase,'_'+element+'_data')[prop][locations] = data
                    phase._logger.debug(element+' property '+prop+' has been '+temp_word+phase.name)
                else:
                    try:
                        getattr(phase,'_'+element+'_data')[prop]
                        if sp.shape(getattr(phase,'_'+element+'_data')[prop])[0]!=1:
                            print('Warning: '+prop+' '+element+' property in '+phase.name+' was an array which has been overwritten with a scalar value')
                    except: pass
                    getattr(phase,'_'+element+'_data')[prop] = data  
                    phase._logger.debug(element+' property '+prop+' has been '+temp_word+phase.name)
            else:
                if locations!='':
                    if sp.shape(locations)[0]==sp.shape(data)[0]:
                        try: getattr(phase,'_'+element+'_data')[prop]
                        except: getattr(phase,'_'+element+'_data')[prop] = sp.zeros((getattr(phase,'num_'+element+'s')(),))*sp.nan
                        getattr(phase,'_'+element+'_data')[prop][locations] = data
                        phase._logger.debug(element+' property '+prop+' has been '+temp_word+phase.name)
                    else:
                        phase._logger.error('For adding '+element+' property '+prop+' to '+phase.name+', locations and size of data do not match!')
                else:
                    try:
                        getattr(phase,'num_'+element+'s')()                        
                        if sp.shape(data)[0]==getattr(phase,'num_'+element+'s')():
                            getattr(phase,'_'+element+'_data')[prop] = data
                            phase._logger.debug(element+' property '+prop+' has been '+temp_word+phase.name)
                        else: phase._logger.error('For adding '+element+' property '+prop+' to '+phase.name+', number of '+element+'s and size of data do not match!')
                    except: phase._logger.error(element+' numbering has not been specified for '+phase.name)
                        
        else:
            try: 
                getattr(self,'_'+element+'_data')[prop]
                temp_word = 'updated for '
            except: temp_word = 'added to '            
            if sp.shape(data)[0]==1:
                if locations!='':                
                    try: getattr(self,'_'+element+'_data')[prop]
                    except: getattr(self,'_'+element+'_data')[prop] = sp.zeros((getattr(self,'num_'+element+'s')(),))*sp.nan
                    getattr(self,'_'+element+'_data')[prop][locations] = data
                    self._logger.debug(element+' property '+prop+' has been '+temp_word+self.name)
                else:
                    try: 
                        getattr(self,'_'+element+'_data')[prop]
                        if sp.shape(getattr(self,'_'+element+'_data')[prop])[0]!=1:
                            print('Warning: '+prop+' '+element+' property in '+self.name+' was an array which has been overwritten with a scalar value')
                    except: pass
                    getattr(self,'_'+element+'_data')[prop] = data   
                    self._logger.debug(element+' property '+prop+' has been '+temp_word+self.name)
            else:                
                if locations!='':
                    if sp.shape(locations)[0]==sp.shape(data)[0]:
                        try: getattr(self,'_'+element+'_data')[prop]
                        except: getattr(self,'_'+element+'_data')[prop] = sp.zeros((getattr(self,'num_'+element+'s')(),))*sp.nan
                        getattr(self,'_'+element+'_data')[prop][locations] = data
                        self._logger.debug(element+' property '+prop+' has been '+temp_word+self.name)
                    else: self._logger.error('For adding '+element+' property '+prop+' to '+self.name+', locations and size of data do not match!')
                else:
                    try: 
                        getattr(self,'num_'+element+'s')()                        
                        if sp.shape(data)[0]==getattr(self,'num_'+element+'s')():
                            getattr(self,'_'+element+'_data')[prop] = data
                            self._logger.debug(element+' property '+prop+' has been '+temp_word+self.name)
                        else: self._logger.error('For adding '+element+' property '+prop+' to '+self.name+', number of '+element+'s and size of data do not match!')
                    except: 
                        getattr(self,'_'+element+'_data')[prop] = data
                        self._logger.debug(element+' property '+prop+' has been '+temp_word+self.name)



    def _get_data(self,element='',phase='',prop='',locations=''):
        r'''
        Documentation for this method is being updated, we are sorry for the inconvenience.
        '''      
        if type(locations)==list: 
            try: locations = getattr(self,'get_'+element+'_indices')(locations)
            except: locations = sp.array(locations,ndmin=1)
        elif type(locations)==sp.ndarray:
            try: locations = getattr(self,'get_'+element+'_indices')(locations)
            except: pass            
        elif locations!='':
            try: locations = locations.name 
            except: pass
            if type(locations)==str: locations = getattr(self,'get_'+element+'_indices')([locations])
       
        if phase :
            try: phase = self.find_object_by_name(phase) 
            except: pass 
            if locations!='':                
                try: 
                    getattr(phase,'_'+element+'_data')[prop]
                    try: return getattr(phase,'_'+element+'_data')[prop][locations]
                    except: phase._logger.error('data for these locations cannot be returned')
                except: phase._logger.error(phase.name+' does not have the requested '+element+' property: '+prop) 
            else:
                try: return getattr(phase,'_'+element+'_data')[prop]
                except: phase._logger.error(phase.name+' does not have the requested '+element+' property: '+prop) 
       
        else :
            if locations!='':                
                try: 
                    getattr(self,'_'+element+'_data')[prop]
                    try: return getattr(self,'_'+element+'_data')[prop][locations]
                    except: self._logger.error('data for these locations cannot be returned')
                except: self._logger.error(self.name+' does not have the requested '+element+' property: '+prop) 
            else:
                try: return getattr(self,'_'+element+'_data')[prop]
                except: self._logger.error(self.name+' does not have the requested '+element+' property: '+prop)           
      
 
    def set_pore_data(self,phase='',prop='',data='',locations=''):
        r'''
        Writes data to fluid or network objects according to input arguments.
        
        Parameters
        ----------
        prop : string
            Name of property to write
        phase : OpenPNM Fluid object or fluid name string, optional
            Fluid to which data is written.  If omitted data is written to network object.
        data : array_like
            Data values to write to object
        locations: It can be object, location string (or a list of strings), boolean array or indices.   
        See Also
        --------
        set_throat_data, set_pore_info, set_throat_info
        
        Notes
        -----
        This is wrapper method that calls set_data, which is generic for pores and throats
            
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.set_pore_data(prop='test',data=1.1)
        >>> pn.get_pore_data(prop='test')
        array([ 1.1])
        '''
        self._set_data(element='pore',phase=phase,prop=prop,data=data,locations=locations)
        
    def get_pore_data(self,phase='',prop='',locations=''):
        r'''
        Retrieves data from fluid or network objects according to input arguments.
        
        Parameters
        ----------
        prop : string
            Name of property to retrieve.  Requesting property 'all' prints a list of existing properties.
        phase : string, optional
            Name of fluid from which to retrieve data.  If omitted data is retrieved from network object.
        locations: It can be object, location string (or a list of strings), boolean array or indices.   

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
        >>> pn.set_pore_data(prop='test',data=1.1)
        >>> pn.get_pore_data(prop='test')
        array([ 1.1])
        '''
        return self._get_data(element='pore',phase=phase,prop=prop,locations=locations)

    def set_throat_data(self,phase='',prop='',data='',locations=''):
        r'''
        Writes data to fluid or network objects according to input arguments.  
        Network topology data and pore/throat geometry data is stored on the network object.
        Fluid properties and physics properties are stored on the corresponding fluid object.
        
        Parameters
        ----------
        prop : string
            Name of property to write
        phase : OpenPNM fluid object or fluid name string, optional
            Fluid to which data is written.  If omitted data is written to network object.
        data : array_like
            Data values to write to object
        locations: It can be object, location string (or a list of strings), boolean array or indices.   

        See Also
        --------
        set_pore_data, set_pore_info, set_throat_info
            
        Notes
        -----
        This is wrapper method that calls _set_data, which is generic for pores and throats
        
        Examples
        --------
        See set_pore_data
        '''
        self._set_data(element='throat',phase=phase,prop=prop,data=data,locations=locations)         

    def get_throat_data(self,phase='',prop='',locations=''):
        r'''
        Retrieves data from fluid or network objects according to input arguments.
        
        Parameters
        ----------
        prop : string
            Name of property to retrieve.  Requesting property 'all' prints a list of existing properties.
        phase : string, optional
            Name of fluid from which to retrieve data.  If omitted data is retrieved from network object.
        locations: It can be geometry object, location string (or a list of strings), boolean array or indices.   

        Returns
        -------
        array_like
            An ndarray containing the requested property data from the specified object
            
        See Also
        --------
        get_pore_data, get_pore_info, get_throat_info

        Notes
        -----
        This is a wrapper method that calls _get_data, which is generic for pores and throats
        
        Examples
        --------
        See get_pore_data
        '''
        return self._get_data(element='throat',phase=phase,prop=prop,locations=locations)     

    def _set_info(self,element='',label='',locations='',mode='merge'):
        r'''
        This is the actual info setter method, but it should not be called directly.  
        Wrapper methods have been created.  Use set_pore_info and get_pore_info.
        
        See Also
        --------
        set_pore_info, set_throat_info
        '''
        if type(locations)==list: 
            try: locations = getattr(self,'get_'+element+'_indices')(locations)
            except: locations = sp.array(locations,ndmin=1)
        elif type(locations)==sp.ndarray:
            try: locations = getattr(self,'get_'+element+'_indices')(locations)
            except : pass 
        if locations!='':
            
            try: 
                locations = locations.name
                label = locations
            except: pass
            if type(locations)==str: locations = getattr(self,'get_'+element+'_indices')([locations])           
            locations=sp.array(locations,ndmin=1)
            if label:
                if label=='all':
                    try: 
                        old_label = getattr(self,'_'+element+'_info')[label]
                        if sp.shape(old_label)[0]<sp.shape(locations)[0]:
                            getattr(self,'_'+element+'_info')[label] = sp.ones_like(locations,dtype=bool)
                            self._logger.info('label=all for '+element+'has been updated to a bigger size!')
                            for info_labels in getattr(self,'_'+element+'_info').keys():
                                if info_labels!=label:
                                    temp = sp.zeros((getattr(self,'num_'+element+'s')(),),dtype=bool)
                                    temp[old_label] = getattr(self,'_'+element+'_info')[info_labels]
                                    getattr(self,'_'+element+'_info')[info_labels] = temp
                        elif sp.shape(old_label)[0]>sp.shape(locations)[0]: 
                            self._logger.error('To apply a new numbering label (label=all) to '+element+'s, size of the locations cannot be less than total number of '+element+'s!!')
                    except: getattr(self,'_'+element+'_info')[label] = sp.ones_like(locations,dtype=bool)
                else:    
                    try: getattr(self,'_'+element+'_info')[label]
                    except: getattr(self,'_'+element+'_info')[label] = sp.zeros((getattr(self,'num_'+element+'s')(),),dtype=bool)
                    if mode=='overwrite':
                        getattr(self,'_'+element+'_info')[label] = sp.zeros((getattr(self,'num_'+element+'s')(),),dtype=bool)
                        getattr(self,'_'+element+'_info')[label][locations] = True
                    elif mode=='remove':  getattr(self,'_'+element+'_info')[label][locations] = False                           
                    elif mode=='merge':  getattr(self,'_'+element+'_info')[label][locations] = True
            else: self._logger.error('No label has been defined for these locations')                

        elif mode=='remove':  del getattr(self,'_'+element+'_info')[label]
        else:  getattr(self,'_'+element+'_info')[label] = sp.zeros((getattr(self,'num_'+element+'s')(),),dtype=bool)

    def _get_info(self,element='',label='',return_indices=False):
        r'''
        This is the actual info getter method, but it should not be called directly.  
        Wrapper methods have been created.  Use get_pore_info and get_throat_info
        
        See Also
        --------
        get_pore_info, get_throat_info
        
        '''
        if return_indices:
            return sp.where(getattr(self,'_'+element+'_info')[label]==True)[0]
        else:
            return getattr(self,'_'+element+'_info')[label]
           
    def set_pore_info(self,label='',locations='',mode='merge'):
        r'''
        Apply a label to a selection of pores.  
        
        Parameters
        ----------
        label : string
            The name of the pore labels you wish to apply (e.g. 'top')
        locaitons : array_like
            An array containing the locations (pores) where the labels should be applied.
            Can be either a boolean mask of Np length with True at labels locations (default), 
            a list of indices where labels should be applied.
        mode : string
            Set the mode to be used for writing labels.  Options are:
            
            * 'merge' : (default) Adds label to specified locations while 
            maintaining pre-existing labels
            
            * 'overwrite' : Adds label to specified locations while 
            removing all pre-existing labels
            
            * 'remove' : Removes labels from specified locations.  If no
            locations are given then this mode will remove the entire label
            from the network.
            
        See Also
        --------
        set_pore_data, set_throat_data, set_throat_info
            
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.set_pore_info(label='test',locations=[0,1]) #Set using index notation
        >>> pn.get_pore_info(label='test',return_indices=True) #Retrieve values as indices
        array([0, 1], dtype=int64)
        >>> loc = sp.zeros((pn.num_pores(),),dtype=bool)
        >>> loc[[0,1]] = True
        >>> pn.set_pore_info(label='test',locations=loc) #Set using boolean mask
        >>> pn.get_pore_info(label='test',return_indices=True) #Retrieve values as indices
        array([0, 1], dtype=int64)
        '''
        self._set_info(element='pore',label=label,locations=locations,mode=mode)

    def get_pore_info(self,label='',return_indices=False):
        r'''
        Retrieves locations where requested label is applies
        
        Parameters
        ----------
        label : string
            The name of the labels you wish to retrieve (e.g. 'top')
            
        return_indices : bool, optional
            This flag indicates that the returned result should be a list of indices
            
        Returns
        -------
        A boolean mask of length Np with True at all locations where labels apply, 
        or a list of indices where labels apply.
        
        See Also
        --------
        get_pore_data, get_throat_data, get_throat_info
        
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> result = pn.get_pore_info(label='top',return_indices=True) #Retrieve values as indices
        >>> result[0:10]
        array([100, 101, 102, 103, 104, 105, 106, 107, 108, 109], dtype=int64)
        >>> result = pn.get_pore_info(label='top') #Retrieve values as boolean mask
        >>> result[97:103]
        array([False, False, False,  True,  True,  True], dtype=bool)
        '''
        return self._get_info(element='pore',label=label,return_indices=return_indices)
        
    def set_throat_info(self,label='',locations='',mode='merge'):
        r'''
        Apply a label to a selection of throats
        
        Parameters
        ----------
        label : string
            The name of the pore labels you wish to apply (e.g. 'top')
        mode : string
            Options are 'merge' and 'overwrite', default is 'merge'
        locaitons : array_like
            An array containing the locations (pores) where the labels should be applied.
            Can be either a boolean mask of Np length with True at labels locations (default), 
            a list of indices where labels should be applied. 
        mode : string
            Set the mode to be used for writing labels.  Options are:
            
            * 'merge' : (default) Adds label to specified locations while 
            maintaining pre-existing labels
            
            * 'overwrite' : Adds label to specified locations while 
            removing all pre-existing labels
            
            * 'remove' : Removes labels from specified locations.  If no
            locations are given then this mode will remove the entire label
            from the network.
            
        See Also
        --------
        set_pore_data, set_throat_data, set_pore_info
            
        Examples
        --------
        See set_pore_info for usage
        '''
        self._set_info(element='throat',label=label,locations=locations,mode=mode)
        
    def get_throat_info(self,label='',return_indices=False):
        r'''
        Retrieves locations where requested labels are applied
        
        Parameters
        ----------
        label : string
            The name of the labels you wish to retrieve (e.g. 'top')
            
        return_indices : bool, optional
            This flag indicates that the returned result should be a list of indices
            
        Returns
        -------
        A boolean mask of length Np with True at all locations where labels apply, 
        or a list of indices where labels apply.
        
        See Also
        --------
        get_pore_data, get_throat_data, get_pore_info
        
        Examples
        --------
        See set_pore_info for usage
        '''
        return self._get_info(element='throat',label=label,return_indices=return_indices)
        
    def list_pore_props(self):
        r'''
        Returns a list containing the names of all defined pore properties.
        This list is an iterable, so is useful for scanning through labels.

        Returns
        -------
        A an alphabetically sorted list containing the string name of all 
        pore properties currently defined.  

        See Also
        --------
        print_dicts
        '''
        temp = list(self._pore_data.keys())
        temp.sort()
        return temp

    def list_throat_props(self):
        r'''
        Returns a list containing the names of all defined throat properties. 
        This list is an iterable, so is useful for scanning through labels.

        Returns
        -------
        A an alphabetically sorted list containing the string name of all 
        throat properties currently defined.  

        See Also
        --------
        print_dicts
        '''
        temp = list(self._throat_data.keys())
        temp.sort()
        return temp
        
    def list_pore_labels(self):
        r'''
        Returns a list containing the names of all defined pore labels. This
        list is an iterable, so is useful for scanning through labels.

        Returns
        -------
        A an alphabetically sorted list containing the string name of all 
        pore labels currently defined.  

        See Also
        --------
        print_dicts
        '''
        temp = list(self._pore_info.keys())
        temp.sort()
        return sp.array(temp,ndmin=1)

    def list_throat_labels(self):
        r'''
        Returns a list containing the names of all defined throat labels. This
        list is an iterable, so is useful for scanning through labels.
        
        Returns
        -------
        A an alphabetically sorted list containing the string name of all 
        throat labels currently defined.  
        
        See Also
        --------
        print_dicts
        '''
        temp = list(self._throat_info.keys())
        temp.sort()
        return sp.array(temp,ndmin=1)
        
    def _get_labels(self,element,locations,mode='union',flatten=False):
        r'''
        This is the actual label getter method, but it should not be called directly.  
        Wrapper methods have been created.  Use get_pore_labels and get_throat_labels
        
        See Also
        --------
        get_pore_labels, get_throat_labels
        
        '''
        if element == 'pore':
            element_info = self._pore_info
            labels = self.list_pore_labels()
        elif element == 'throat':
            element_info = self._throat_info
            labels = self.list_throat_labels()
        arr = sp.ndarray((sp.shape(locations,)[0],len(labels)),dtype=bool)
        col = 0
        for item in labels:
            arr[:,col] = element_info[item][locations]
            col = col + 1
        if flatten == True:
            if mode == 'count':
                return sp.sum(arr,axis=1)
            if mode == 'union':
                return labels[sp.sum(arr,axis=0)>0]
            if mode == 'intersection':
                return labels[sp.sum(arr,axis=0)==sp.shape(locations,)[0]]
        else:
            if mode == 'raw':
                return arr
            else:
                temp = sp.ndarray((sp.shape(locations,)[0],),dtype=object)
                for i in sp.arange(0,sp.shape(locations,)[0]):
                    temp[i] = list(labels[arr[i,:]])
                return temp
            
    def get_pore_labels(self,pnums,mode='union',flatten=False):
        r'''
        Returns the labels applied to specified pore locations
        
        Parameters
        ----------
        pnums : array_like
            The pores whos labels are sought
        flatten : boolean, optional
            If False (default) the returned list is Np long, where each element
            is a list labels applied to the specified pores.  If True, the mode
            logic
        mode : string, optional
            Controls how the query should be performed
            
            * 'union' : A list of labels applied to ANY of the given locations        
            
            * 'intersection' : Label applied to ALL of the given locations
            
            * 'count' : The number of labels on each pore
            
            * 'raw' : returns an Np x Nlabels array, where each row corresponds
            to a pore location, and each column contains the truth value for 
            the existance of labels as returned from list_pore_labels().
            
        Notes
        -----
        The mode argument is ignored unless flatten is True, with the exception 
        of 'raw'.
        
        '''
        pnums = sp.array(pnums,ndmin=1)
        temp = self._get_labels(element='pore',locations=pnums, mode=mode, flatten=flatten)
        return temp

    def get_throat_labels(self,tnums,mode='union',flatten=False):
        r'''
        Returns the labels applied to specified throat locations
        
        Parameters
        ----------
        tnums : array_like
            The throats whos labels are sought
        flatten : boolean, optional
            If false (default) the returned list is Nt long, where each element
            is a list labels applied to the specified throats
        mode : string, optional
            Controls how the query should be performed
            
            * 'union' : A list of labels applied to ANY of the given locations        
            
            * 'intersection' : Label applied to ALL of the given locations
            
            * 'count' : The number of labels on each throat
            
            * 'raw' : returns an Np x Nlabels array, where each row corresponds
            to a throat location, and each column contains the truth value for 
            the existance of labels as returned from list_throat_labels().
            
        Notes
        -----
        The mode argument is ignored unless flatten is True, with the exception 
        of 'raw'.
        
        '''
        tnums = sp.array(tnums,ndmin=1)
        temp = self._get_labels(element='throat',locations=tnums,mode=mode,flatten=flatten)
        return temp
        
    def _get_indices(self,element,labels,return_indices,mode):
        r'''
        This is the actual method for getting indices, but should not be called
        directly.  
        '''
        if mode == 'union':
            union = sp.zeros_like(self._get_info(element=element,label='all'),dtype=bool)
            for item in labels: #iterate over labels list and collect all indices
                    union = union + self._get_info(element=element,label=item)
            ind = union
        elif mode == 'intersection':
            intersect = sp.ones_like(self._get_info(element=element,label='all'),dtype=bool)
            for item in labels: #iterate over labels list and collect all indices
                    intersect = intersect*self._get_info(element=element,label=item)
            ind = intersect
        elif mode == 'not_intersection':
            not_intersect = sp.zeros_like(self._get_info(element=element,label='all'),dtype=int)
            for item in labels: #iterate over labels list and collect all indices
                info = self._get_info(element=element,label=item)
                not_intersect = not_intersect + sp.int8(info)
            ind = (not_intersect == 1)
        elif mode == 'none':
            none = sp.zeros_like(self._get_info(element=element,label='all'),dtype=int)
            for item in labels: #iterate over labels list and collect all indices
                info = self._get_info(element=element,label=item)
                none = none - sp.int8(info)
            ind = (none == 0)
        if return_indices: ind = sp.where(ind==True)[0]
        return ind

    def get_pore_indices(self,labels=['all'],return_indices=True,mode='union'):
        r'''
        Returns pore locations where given labels exist.
        
        Parameters
        ----------
        labels : list of strings, optional
            The pore label(s) whose locations are requested.
            If omitted, all pore inidices are returned.
        indices : boolean, optional
            This flag specifies whether pore locations are returned a boolean mask of length Np,
            or as a list of indices (default).
        mode : string, optional
            Specifies how the query should be performed.  The options are:    

            * 'union' : (default) All pores with ANY of the given labels are 
            returned.

            * 'intersection' : Only pore with ALL the given labels are 
            returned.

            * 'not_intersection' : Only pores with exactly one of the given 
            labels are returned.
            
            * 'none' : Only pores with none of the given labels are returned.
        
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pind = pn.get_pore_indices(labels=['top','front'],mode='union')
        >>> pind[[0,1,2,-3,-2,-1]]
        array([  0,   5,  10, 122, 123, 124], dtype=int64)
        >>> pn.get_pore_indices(labels=['top','front'],mode='intersection')
        array([100, 105, 110, 115, 120], dtype=int64)
        '''
        if type(labels) == str: labels = [labels] #convert string to list, if necessary
        ind = self._get_indices(element='pore',labels=labels,return_indices=return_indices,mode=mode)
        return ind

    def get_throat_indices(self,labels=['all'],return_indices=True,mode='union'):
        r'''
        Returns throat locations where given labels exist.
        
        Parameters
        ----------
        labels : list of strings, optional
            The throat label(s) whose locations are requested.
            If omitted, all throat inidices are returned.
        return_indices : boolean, optional
            This flag specifies whether throat locations are returned as a boolean mask of length Np,
            or as a list of indices (default).
        mode : string, optional
            Specifies how the query should be performed.  The options are: 

            * 'union' : (default) All throats with ANY of the given labels are 
            returned.

            * 'intersection' : Only throats with ALL the given labels are 
            counted.

            * 'not_intersection' : Only throats with exactly one of the given 
            labels are counted.
            
            * 'none' : Only throats with none of the given labels are returned.
        
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> Tind = pn.get_throat_indices()
        >>> Tind[0:5]
        array([0, 1, 2, 3, 4], dtype=int64)
        '''
        if type(labels) == str: labels = [labels] #convert string to list, if necessary
        ind = self._get_indices(element='throat',labels=labels,return_indices=return_indices,mode=mode)
        return ind
        
    def num_pores(self,labels=['all'],mode='union'):
        r'''
        Returns the number of pores of the specified labels

        Parameters
        ----------
        labels : list of strings, optional
            The pore labels that should be included in the count.  
            If not supplied, all pores are counted.
        labels : list of strings
            Label of pores to be returned
        mode : string, optional
            Specifies how the count should be performed.  The options are: 
            
            * 'union' : (default) All pores with ANY of the given labels are counted.

            * 'intersection' : Only pores with ALL the given labels are counted.
            
            * 'not_intersection' : Only pores with exactly one of the given labels are counted.
            
            * 'none' : Only pores with none of the given labels are counted.
            
        Returns
        -------
        Np : int
            Number of pores with the specified labels 
            
        See Also
        --------
        num_throats
            
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.num_pores()
        125
        >>> pn.num_pores(labels=['top'])
        25
        >>> pn.num_pores(labels=['top','front'],mode='union') #'union' is default
        45
        >>> pn.num_pores(labels=['top','front'],mode='intersection')
        5
        >>> pn.num_pores(labels=['top','front'],mode='not_intersection')
        40
        
        '''
        #convert string to list, if necessary
        if type(labels) == str: labels = [labels]
        #Count number of pores of specified type
        temp = self.get_pore_indices(labels=labels,mode=mode,return_indices=False)
        return sp.sum(temp) #return sum of Trues
            
    def num_throats(self,labels=['all'],mode='union'):
        r'''
        Return the number of throats of the specified labels

        Parameters
        ----------
        labels : list of strings, optional
            The throat labels that should be included in the count.  
            If not supplied, all throats are counted.
        mode : string, optional
            Specifies how the count should be performed.  The options are: 

            * 'union' : (default) All throats with ANY of the given labels are counted.

            * 'intersection' : Only throats with ALL the given labels are counted.

            * 'not_intersection' : Only throats with exactly one of the given labels are counted.
            
            * 'none' : Only throats with none of the given labels are counted.

        Returns
        -------
        Nt : int
            Number of throats with the specified labels
            
        See Also
        --------
        num_pores

        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.num_throats()
        300
        >>> pn.num_throats(labels=['top'])
        40
        >>> pn.num_throats(labels=['top','front'],mode='union') #'union' is default
        76
        >>> pn.num_throats(labels=['top','front'],mode='intersection')
        4
        >>> pn.num_throats(labels=['top','front'],mode='not_intersection')
        72
        
        '''
        #convert string to list, if necessary
        if type(labels) == str: labels = [labels]
        #Count number of pores of specified type
        temp = self.get_throat_indices(labels=labels,mode=mode,return_indices=False)
        return sp.sum(temp) #return sum of Trues

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

