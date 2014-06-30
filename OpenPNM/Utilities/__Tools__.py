"""
module __Tools__: Base class to construct pore network tools
==================================================================

.. warning:: The classes of this module should be loaded through the 'Base.__init__.py' file.

"""

import sys, os,pprint
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if sys.path[1] != parent_dir:
    sys.path.insert(1, parent_dir)
import scipy as sp
import OpenPNM
from OpenPNM.Utilities import Base


class Tools(Base,dict):
    r'''
    This class contains tools to read and write data in OpenPNM objects

    '''
    def __init__(self, **kwargs):
        r'''
        Initialize
        '''
        super(Tools,self).__init__(**kwargs)
        self._logger.info("Construct Base class from Tools subclass")
        #Initialize network properties dictionaries
        self._logger.debug("Construction of Tools class complete")
        
    def __setitem__(self,key,value):
        r'''
        This is a subclass of the default __setitem__ behavior.  The main aim
        is to limit what type and shape of data can be written to protect
        the integrity of the network.
        '''
        element = key.split('.')[0]
        if type(value) == int:
            value = [value]
        value = sp.array(value,ndmin=0)
        if (element != 'pore') and (element != 'throat'):
            self._logger.error('Array name must begin with \'pore\' or \'throat\'')
            return
        if (key == 'pore.coords') or (key == 'throat.conns'):
            super(Base, self).__setitem__(key,value)
            return
        elif key.split('.')[1] == 'all':
            try: 
                self[key]
                self._logger.error(key+' is already defined.')
            except: super(Base, self).__setitem__(key,value)
            return
        else:
            if key not in self:
                if (sp.shape(value)[0] == 1):
                    self._logger.debug('Adding scalar: '+key)
                    super(Base, self).__setitem__(key,value)
                elif (sp.shape(value)[0] == self.num_pores()):
                    self._logger.debug('Adding Np length vector: ' +key)
                    super(Base, self).__setitem__(key,value)
                elif (sp.shape(value)[0] == self.num_throats()):
                    self._logger.debug('Adding Nt length vector: '+key)
                    super(Base, self).__setitem__(key,value)
                else:
                    self._logger.error('Cannot create a label or property of the wrong length')
            else:
                if (sp.shape(value)[0] == 1):
                    if sp.shape(self[key]) == 1:
                        self._logger.debug('Updating scalar: '+key)
                        super(Base, self).__setitem__(key,value)
                    else:
                        self._logger.debug('Overwriting vector with scalar: '+key)
                        super(Base, self).__setitem__(key,value)
                else:
                    if sp.shape(value)[0] == sp.shape(self[key])[0]:
                        self._logger.debug('Updating vector :'+key)
                        super(Base, self).__setitem__(key,value)
                    else:
                        if sp.shape(value)[0] == self.num_pores():
                            self._logger.debug('Updating vector: '+key)
                            super(Base, self).__setitem__(key,value)
                        elif sp.shape(value)[0] == self.num_throats():
                            self._logger.debug('Updating vector: '+key)
                            super(Base, self).__setitem__(key,value)
                        else:
                            self._logger.error('Cannot overwrite '+key+' with an array of the wrong length')
    
#    def __getitem__(self,key):
#        try:
#            temp = super(Base, self).__getitem__(key)
#        except:
#            self._logger.warning('Requested array will be created')
#            self.__setitem__(key,[])
#            temp = super(Base, self).__getitem__(key)
#        return temp
    
    #--------------------------------------------------------------------------
    '''Setter and Getter Methods'''
    #--------------------------------------------------------------------------
    def _set_data(self,element='',prop='',data='',locations='',mode='merge'):
        r'''
        Documentation for this method is being updated, we are sorry for the inconvenience.
        '''
#        if locations = 'all'
#            locations 
#        self[element+'.'+prop] = data
        if mode=='remove':
            if data=='':
                try: 
                    self[element+'.'+prop]
                    if locations!='':   
                        self[element+'.'+prop][locations] = sp.nan
                        self._logger.debug('For the '+element+' property '+prop+', the specified data have been deleted in '+self.name)                        
                    else:
                        del self[element+'.'+prop]
                        self._logger.debug(element+' property '+prop+' has been deleted from the dictionary in '+self.name)
                except: self._logger.error(element+' property '+prop+' in '+self.name+' has not been defined. Therefore, no data can be removed!')                              
            else:   self._logger.error('For the '+element+' property '+prop+' in '+self.name+': The (remove) mode, will remove the property from the dictionary or specified locations. No data should be sent!')        
        else:
            data = sp.array(data,ndmin=1)
            if data.ndim > 1: data = data.squeeze()
            if 'OpenPNM.Network' in str(self.__class__): net = self
            else: net = self._net
            if type(locations)==list:
                try: locations = getattr(net,'get_'+element+'_indices')(locations)
                except: locations = sp.array(locations,ndmin=1)
            elif type(locations)==sp.ndarray:
                try: locations = getattr(net,'get_'+element+'_indices')(locations)
                except: pass
            elif locations!='':
                try: locations = locations.name
                except: pass
                if type(locations)==str: locations = getattr(net,'get_'+element+'_indices')([locations])
            try:
                self[element+'.'+prop]
                temp_word = 'updated for '
            except: temp_word = 'added to '            
            if sp.shape(data)[0]==1:
                if locations!='':                
                    try: 
                        self[element+'.'+prop]
                        if mode=='overwrite':
                            self[element+'.'+prop] = sp.zeros((getattr(self,'num_'+element+'s')(),))*sp.nan
                        elif mode=='merge' and sp.shape(self[element+'.'+prop])[0]==1:
                            self[element+'.'+prop] = self[element+'.'+prop]*sp.ones((getattr(self,'num_'+element+'s')(),))
                    except: self[element+'.'+prop] = sp.zeros((getattr(self,'num_'+element+'s')(),))*sp.nan
                    self[element+'.'+prop][locations] = data
                    self._logger.debug(element+' property '+prop+' has been '+temp_word+self.name)
                else:
                    try: 
                        self[element+'.'+prop]
                        if mode=='overwrite':
                            if sp.shape(self[element+'.'+prop])[0]!=1:
                                self._logger.debug(element+' property '+prop+' in '+self.name+' was an array which has been overwritten with a scalar value')
                            self[element+'.'+prop] = data
                            self._logger.debug(element+' property '+prop+' has been '+temp_word+self.name)
                        if mode=='merge' and sp.shape(self[element+'.'+prop])[0]!=1:  
                            self._logger.error('a scalar data without specified locations cannot be merged into the '+element+' property '+prop+' in '+self.name+' which is (1*N) array. To do so, choose overwrite mode.')
                    except:
                        self[element+'.'+prop] = data
                        self._logger.debug(element+' property '+prop+' has been '+temp_word+self.name)
            else:                
                if locations!='':
                    if sp.shape(locations)[0]==sp.shape(data)[0]:
                        try:                                 
                            self[element+'.'+prop]
                            if mode=='overwrite':
                                self[element+'.'+prop] = sp.zeros((getattr(self,'num_'+element+'s')(),))*sp.nan
                        except: self[element+'.'+prop] = sp.zeros((getattr(self,'num_'+element+'s')(),))*sp.nan                            
                        self[element+'.'+prop][locations] = data
                        self._logger.debug(element+' property '+prop+' has been '+temp_word+self.name)
                    else: self._logger.error('For adding '+element+' property '+prop+' to '+self.name+', locations and size of data do not match!')
                else:
                    try: 
                        getattr(self,'num_'+element+'s')()                        
                        if sp.shape(data)[0]==getattr(self,'num_'+element+'s')():
                            self[element+'.'+prop] = data
                            self._logger.debug(element+' property '+prop+' has been '+temp_word+self.name)
                        else: self._logger.error('For adding '+element+' property '+prop+' to '+self.name+', no locations have been specified. Also size of the data and number of '+element+'s do not match! To add this property, specify locations or change the size of the data.')
                    except: 
                        self[element+'.'+prop] = data
                        self._logger.debug(element+' property '+prop+' has been '+temp_word+self.name)

    def _get_data(self,element='',prop='',locations='',mode=''):
        r'''
        Documentation for this method is being updated, we are sorry for the inconvenience.
        '''
        if self.__module__.split('.')[1] == 'Network': net = self
        else: net = self._net
        if type(locations)==list: 
            try: locations = getattr(net,'get_'+element+'_indices')(locations)
            except: locations = sp.array(locations,ndmin=1)
        elif type(locations)==sp.ndarray:
            try: locations = getattr(net,'get_'+element+'_indices')(locations)
            except: pass            
        elif locations!='':
            try: locations = locations.name 
            except: pass
            if type(locations)==str: locations = getattr(net,'get_'+element+'_indices')([locations])    
        if locations!='':
            if mode != '':
                if mode == 'interpolate':
                    if element == 'pore':
                        return getattr(self,'interpolate_data')(prop=prop,pores=locations)
                    else:
                        return getattr(self,'interpolate_data')(prop=prop,throats=locations)
                else:
                    self._logger.error('The requested mode '+mode+' is not valid')
            else:
                try: 
                    self[element+'.'+prop]
                    try: 
                        if  sp.shape(self[element+'.'+prop])[0]==1 and max(locations)<getattr(self,'num_'+element+'s')():
                            return self[element+'.'+prop]                        
                        else: return self[element+'.'+prop][locations]
                    except: self._logger.error('data for these locations cannot be returned')
                except: 
                    self._logger.error(self.name+' does not have the requested '+element+' property: '+prop) 
        else:
            try: return self[element+'.'+prop]
            except: self._logger.error(self.name+' does not have the requested '+element+' property: '+prop)           
 
    def get_data(self,prop='',pores=None,throats=None,mode=''):
        r'''
        Retrieves data from the object to which it is associated according to 
        input arguments.
        
        Parameters
        ----------
        prop : string
            Name of property to retrieve
        pores : array_like, or string 'all'
            List of pore indices from which to retrieve data.  If set to 'all',
            then ALL values are returned
        throats : array_like, or string 'all'
            List of throat indies from which to retrieve data.  If set to 'all'
            , then ALL values are returned.
        modes : string
            None yet

        Returns
        -------
        array_like
            An ndarray containing the requested property data from the 
            specified locations
            
        See Also
        --------
        get_info

        Notes
        -----
        This is a wrapper method that calls _get_data.  
        Only one of pores or throats should be sent.
            
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.set_data(prop='test',pores=[0],data=1.1)
        >>> pn.get_data(prop='test',pores=[0])
        array([ 1.1])
        '''
        if pores != None:
            if pores == 'all':
                pores = self.get_pore_indices(labels='all')
            return self._get_data(element='pore',prop=prop,locations=pores,mode=mode)
        if throats != None:
            if throats == 'all':
                throats = self.get_throat_indices(labels='all')
            return self._get_data(element='throat',prop=prop,locations=throats,mode=mode)
            
    def set_data(self,prop='',data='',pores=None,throats=None,mode='merge'):
        r'''
        Write data according to input arguments.
        
        Parameters
        ----------
        prop : string
            Name of property to retrieve.
        pores : array_like, or string 'all'
            List of pore indices to which data will be written.
        throats : array_like, or string 'all'
            List of throat indices to which data will be written.
            
        See Also
        --------
        set_info
        
        Notes
        -----
        This is a wrapper method that calls _set_data.  
        Only one of pores or throats should be sent.
            
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.set_data(prop='test',pores=[0],data=1.1)
        >>> pn.get_data(prop='test',pores=[0])
        array([ 1.1])
        '''
        data= sp.array(data,ndmin=1)
        if pores != None:
            if pores == 'all':
                if sp.shape(data)[0] == 1:
                    pores = ''
                else:
                    pores = self.get_pore_indices(labels='all')
            self._set_data(element='pore',prop=prop,data=data,locations=pores,mode=mode)
        if throats != None:
            if throats == 'all':
                if sp.shape(data)[0] == 1:
                    throats = ''
                else:
                    throats = self.get_throat_indices(labels='all')
            self._set_data(element='throat',prop=prop,data=data,locations=throats,mode=mode)    

    def set_pore_data(self,prop='',data='',locations='',mode='merge'):
        r'''
        THIS METHOD IS DEPRECATED, USE set_data INSTEAD
        '''
        self._set_data(element='pore',prop=prop,data=data,locations=locations,mode=mode)
        
    def get_pore_data(self,prop='',locations=''):
        r'''
        THIS METHOD IS DEPRECATED, USE get_data INSTEAD
        '''
        return self._get_data(element='pore',prop=prop,locations=locations)

    def set_throat_data(self,prop='',data='',locations='',mode='merge'):
        r'''
        THIS METHOD IS DEPRECATED, USE set_data INSTEAD        
        '''
        self._set_data(element='throat',prop=prop,data=data,locations=locations,mode=mode)         

    def get_throat_data(self,prop='',locations=''):
        r'''
        THIS METHOD IS DEPRECATED, USE get_data INSTEAD
        '''
        return self._get_data(element='throat',prop=prop,locations=locations)     

    def _set_info(self,element='',label='',locations='',mode='merge'):
        r'''
        This is the actual info setter method, but it should not be called directly.  
        Wrapper methods have been created.  Use set_info().
        '''
        locations = sp.array(locations,ndmin=1)
        if type(locations)==sp.ndarray:
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
                        old_label = self[element+'.'+label]
                        if sp.shape(old_label)[0]<sp.shape(locations)[0]:
                            self[element+'.'+label] = sp.ones_like(locations,dtype=bool)
                            self._logger.info('label=all for '+element+'has been updated to a bigger size!')
                            for info_labels in getattr(self,'_'+element+'_info').keys():
                                if info_labels!=label:
                                    temp = sp.zeros((getattr(self,'num_'+element+'s')(),),dtype=bool)
                                    temp[old_label] = getattr(self,'_'+element+'_info')[info_labels]
                                    getattr(self,'_'+element+'_info')[info_labels] = temp
                        elif sp.shape(old_label)[0]>sp.shape(locations)[0]: 
                            self._logger.error('To apply a new numbering label (label=all) to '+element+'s, size of the locations cannot be less than total number of '+element+'s!!')
                    except: self[element+'.'+label] = sp.ones_like(locations,dtype=bool)
                else:    
                    try: self[element+'.'+label]
                    except: self[element+'.'+label] = sp.zeros((getattr(self,'num_'+element+'s')(),),dtype=bool)
                    if mode=='overwrite':
                        self[element+'.'+label] = sp.zeros((getattr(self,'num_'+element+'s')(),),dtype=bool)
                        self[element+'.'+label][locations] = True
                    elif mode=='remove':  self[element+'.'+label][locations] = False
                    elif mode=='merge':  self[element+'.'+label][locations] = True
            else: self._logger.error('No label has been defined for these locations')                

        elif mode=='remove':  del self[element+'.'+label]
        else:  self[element+'.'+label] = sp.zeros((getattr(self,'num_'+element+'s')(),),dtype=bool)

    def _get_info(self,element='',label='',mode=''):
        r'''
        This is the actual info getter method, but it should not be called directly.  
        Wrapper methods have been created.  Use get_info().        
        '''
        #Clean up label argument to remove leading 'pore' or 'throat'
        if label.split('.')[0] == element:
            label = label.split('.')[1]
        return self[element+'.'+label]
            
    def set_info(self,label='',pores=None,throats=None,mode='merge'):
        r'''
        Apply a label to a selection of pores or throats
        
        Parameters
        ----------
        label : string
            The name of the pore labels you wish to apply (e.g. 'top')
        pores or throats : array_like
            An array containing the pore (or throat) indices where the labels 
            should be applied.  Can be either a boolean mask of Np length with 
            True at labels locations (default), a list of indices where labels 
            should be applied.
        mode : string
            Set the mode to be used for writing labels.  Options are:
            
            * 'merge' : (default) Adds label to specified locations while 
            maintaining pre-existing labels
            
            * 'overwrite' : Adds label to specified locations while 
            removing all pre-existing labels
            
            * 'remove' : Removes labels from specified locations.  If no
            locations are given then this mode will remove the entire label
            from the network.
        '''
        if pores != None:
            if pores == 'all':
                pores = self.pores(labels='all')
            self._set_info(element='pore',label=label,locations=pores,mode=mode)
        if throats != None:
            if throats == 'all':
                throats = self.get_throat_indices(labels='all')
            self._set_info(element='throat',label=label,locations=throats,mode=mode)
            
    def get_info(self,label='',pores=None,throats=None):
        r'''
        Retrieves the locations where the specified label is applied
        
        Parameters
        ----------
        label : string
            Label of interest
        pores (or throats) : array_like
            List of pores or throats
            
        See Also
        --------
        get_labels
            
        '''
        if pores != None:
            if pores == 'all':
                pores = self.pores(labels='all')
            temp = self._get_info(element='pore',label=label)
            return temp[sp.in1d(temp,pores)]
        if throats != None:
            if throats == 'all':
                throats = self.throats(labels='all')
            temp = self._get_info(element='throat',label=label)
            return temp[sp.in1d(temp,throats)]
           
    def set_pore_info(self,label='',locations='',mode='merge'):
        r'''
        THIS METHOD IS DEPRECATED, USE set_info INSTEAD
        '''
        self._set_info(element='pore',label=label,locations=locations,mode=mode)

    def get_pore_info(self,label=''):
        r'''
        THIS METHOD IS DEPRECATED, USE get_info INSTEAD
        '''
        return self._get_info(element='pore',label=label)
        
    def set_throat_info(self,label='',locations='',mode='merge'):
        r'''
        THIS METHOD IS DEPRECATED, USE set_info INSTEAD
        '''
        self._set_info(element='throat',label=label,locations=locations,mode=mode)
        
    def get_throat_info(self,label=''):
        r'''
        THIS METHOD IS DEPRECATED, USE get_info INSTEAD
        '''
        return self._get_info(element='throat',label=label)
        
    def amalgamate_data(self,objs=[]):
        r"""
        Returns a dictionary containing ALL pore data from all netowrk and/or
        fluid objects received as arguments
        """
        if type(objs) != list:
            objs = list(objs)
        data_amalgamated = {}
        for item in objs:
            try:
                for key in item.keys():
                    if sp.amax(item[key]) < sp.inf:
                        dict_name = item.name+'.'+key
                        data_amalgamated.update({dict_name : item[key]})
            except: 
                self._logger.error('Only network and fluid items contain data')
                
        #Add geometry labels as pore values for Paraview plotting
        geoms = self.find_object(obj_type='Geometry')
        data_amalgamated[self.name+'.pore.geometry'] = sp.ones((self.num_pores(),))*sp.nan
        index = 0;
        for item in geoms:
            index = index + 1
            data_amalgamated[self.name+'.pore.geometry'][item.pores()] = index
        return data_amalgamated
        
    def _get_props(self,mode='all'):
        r'''
        This is the actual prop list getter method, but it should not be
        called directly.  Wrapper methods have been created, use props().
        '''
        props = []
        temp = []
        for item in self.keys():
            if self[item].dtype != bool:
                temp.append(item)
        if mode == 'all':
            props = temp
        if mode == 'vectors':
            for item in temp:
                try: 
                    self[item][1]
                    props.append(item)
                except: pass
        if mode == 'scalars':
            for item in temp:
                try: self[item][1]
                except: props.append(item)
        props.sort()
        return props
            
    def props(self,element='',pores=[],throats=[],mode='all'):
        r'''
        Returns a list containing the names of all defined pore or throat
        properties. 
        
        Parameters
        ----------
        pores or throats : array_like
            hmmm
        mode : string, optional
            Set the mode to be used for retrieving props.  Options are:
            
            * 'all' : (default) Returns all pore or throat properties
            
            * 'scalars' : Only return properties that are stored as scalars
            
            * 'vectors' : Only return properties that are stored as vectors

        Returns
        -------
        A an alphabetically sorted list containing the string name of all 
        pore or throat properties currently defined.  This list is an iterable,
        so is useful for scanning through properties.

        See Also
        --------
        labels
        
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.props()
        ['pore.coords', 'throat.conns']
        >>> pn.props('pore')
        ['pore.coords']
        '''
        
        props = self._get_props(mode=mode)
        if (pores == []) and (throats == []):
            if element == '':
                return props            
            elif element == 'pore':
                temp = [item for item in props if item.split('.')[0]=='pore']
            elif element == 'throat':
                temp = [item for item in props if item.split('.')[0]=='throat']
            else:
                self._logger.error('Unrecognized element')
                return
            return temp
        elif pores != []:
            temp = {}
            for item in props:
                if item.split('.')[0] == 'pore':
                    if sp.shape(self[item])[0] == 1: #is scalar
                        vals = sp.ones((sp.shape(pores)[0],))*self[item]
                    else:
                        vals = self[item][pores]
                    temp.update({item:vals})
            return temp
        elif throats != []:
            temp = {}
            for item in props:
                if item.split('.')[0] == 'throat':
                    if sp.shape(self[item])[0] == 1: #is scalar
                        vals = sp.ones((sp.shape(throats)[0],))*self[item]
                    else:
                        vals = self[item][throats]
                    temp.update({item:vals})
            return temp
            
    def _get_labels(self,element='',locations=[],mode='union'):
        r'''
        This is the actual label getter method, but it should not be called directly.  
        Wrapper methods have been created, use get_labels().
        '''
        labels = []
        for item in self.keys():
            if item.split('.')[0] == element:
                if self[item].dtype == bool:
                    labels.append(item)
        labels.sort()
        if locations == []:
            return labels
        else:
            labels = sp.array(labels)
            arr = sp.zeros((sp.shape(locations)[0],len(labels)),dtype=bool)
            col = 0
            for item in labels:
                arr[:,col] = self[item][locations]
                col = col + 1
            if mode == 'count':
                return sp.sum(arr,axis=1)
            if mode == 'union':
                temp = labels[sp.sum(arr,axis=0)>0]
                return temp.tolist()
            if mode == 'intersection':
                temp = labels[sp.sum(arr,axis=0)==sp.shape(locations,)[0]]
                return temp.tolist()
            if mode == 'difference':
                temp = labels[sp.sum(arr,axis=0)!=sp.shape(locations,)[0]]
                return temp.tolist()
            if mode == 'mask':
                return arr
            if mode == 'none':
                temp = sp.ndarray((sp.shape(locations,)[0],),dtype=object)
                for i in sp.arange(0,sp.shape(locations,)[0]):
                    temp[i] = list(labels[arr[i,:]])
                return temp
            else:
                print('unrecognized mode')
                
    def labels(self,element='',pores=[],throats=[],mode='union'):
        r'''
        Returns the labels applied to specified pore locations
        
        Parameters
        ----------
        pores (or throats) : array_like
            The pores (or throats) whose labels are sought.  If left empty a 
            dictionary containing all pore and throat labels is returned.
        mode : string, optional
            Controls how the query should be performed
            
            * 'none' : An N x Li list of all labels applied to each input pore (or throats). Li can vary betwen pores (and throats)
            
            * 'union' : A list of labels applied to ANY of the given pores (or throats)
            
            * 'intersection' : Label applied to ALL of the given pores (or throats)
            
            * 'difference' : Labels NOT applied to ALL pores (or throats)
            
            * 'count' : The number of labels on each pores (or throats)
            
            * 'mask' : returns an N x Lt array, where each row corresponds to a pore (or throat) location, and each column contains the truth value for the existance of labels as returned from labels(pores='all',mode='union')).
            
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.labels(pores=[0,1,5,6])
        ['pore.all', 'pore.bottom', 'pore.front', 'pore.internal', 'pore.left']
        >>> pn.labels(pores=[0,1,5,6],mode='intersection')
        ['pore.all', 'pore.bottom', 'pore.internal']
        '''
        if (pores == []) and (throats == []):
            if element == '':
                temp = []
                temp = self._get_labels(element='pore')
                temp = temp + self._get_labels(element='throat')
            elif element == 'pore':
                temp = self._get_labels(element='pore',locations=[], mode=mode)
            elif element == 'throat':
                temp = self._get_labels(element='throat',locations=[], mode=mode)
            else:
                self._logger.error('Unrecognized element')
                return
            return temp
        elif pores != []:
            if pores == 'all':
                pores = self.pores()
            pores = sp.array(pores,ndmin=1)
            temp = self._get_labels(element='pore',locations=pores, mode=mode)
            return temp
        elif throats != []:
            if throats == 'all':
                throats = self.throats()
            throats = sp.array(throats,ndmin=1)
            temp = self._get_labels(element='throat',locations=throats,mode=mode)
            return temp
            
    def filter_by_label(self,pores=[],throats=[],label=''):
        r'''
        Returns which of the supplied pores (or throats) has the specified label
        
        Parameters
        ----------
        pores, or throats : array_like
            List of pores or throats to be filtered
            
        label : string
            The label to apply as a filter
            
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.filter_by_label(pores=[0,1,5,6],label='left')
        array([0,1])
        '''
        if pores != []:
            label = 'pore.'+label.split('.')[-1]
            all_labels = self.labels('pore')
            mask = self.labels(pores=pores,mode='mask')
            ind = all_labels.index(label)
            temp = mask[:,ind]
            pores = sp.array(pores,ndmin=1)
            return pores[temp]
        elif throats != []:
            label = 'throat.'+label.split('.')[-1]
            all_labels = self.labels('throat')
            mask = self.labels(throats=throats,mode='mask')
            ind = all_labels.index(label)
            temp = mask[:,ind]
            throats = sp.array(throats,ndmin=1)
            return throats[temp]            
        
    def _get_indices(self,element,labels,mode):
        r'''
        This is the actual method for getting indices, but should not be called
        directly.  Use pores or throats instead.
        '''
        try: labels = [labels.name]  # Check if object was sent
        except: pass
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
        elif mode == 'difference':
            none = sp.zeros_like(self._get_info(element=element,label='all'),dtype=int)
            for item in labels: #iterate over labels list and collect all indices
                info = self._get_info(element=element,label=item)
                none = none - sp.int8(info)
            ind = (none == 0)
        #Extract indices from boolean mask
        ind = sp.where(ind==True)[0]
        return ind
        
    def pores(self,labels='all',mode='union'):
        r'''
        Returns pore locations where given labels exist.
        
        Parameters
        ----------
        labels : list of strings, optional
            The pore label(s) whose locations are requested.
            If omitted, all pore inidices are returned.
        mode : string, optional
            Specifies how the query should be performed.  The options are:    

            * 'union' : (default) All pores with ANY of the given labels are 
            returned.

            * 'intersection' : Only pore with ALL the given labels are 
            returned.

            * 'not_intersection' : Only pores with exactly one of the given 
            labels are returned.
            
            * 'difference' : Only pores with none of the given labels are returned.
        
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
        ind = self._get_indices(element='pore',labels=labels,mode=mode)
        return ind
        
    def throats(self,labels='all',mode='union'):
        r'''
        Returns throat locations where given labels exist.
        
        Parameters
        ----------
        labels : list of strings, optional
            The throat label(s) whose locations are requested.
            If omitted, 'all' throat inidices are returned.
        mode : string, optional
            Specifies how the query should be performed.  The options are: 

            * 'union' : (default) All throats with ANY of the given labels are 
            returned.

            * 'intersection' : Only throats with ALL the given labels are 
            counted.

            * 'not_intersection' : Only throats with exactly one of the given 
            labels are counted.
            
            * 'difference' : Only throats with none of the given labels are returned.
        
        Notes
        -----
        This method replaces get_throat_indices
        
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> Tind = pn.throats()
        >>> Tind[0:5]
        array([0, 1, 2, 3, 4], dtype=int64)
        
        '''
        if type(labels) == str: labels = [labels] #convert string to list, if necessary
        ind = self._get_indices(element='throat',labels=labels,mode=mode)
        return ind
        
    def get_pore_indices(self,labels=['all'],mode='union'):
        r'''
        THIS METHOD IS DEPRECATED, USE pores() INSTEAD
        '''
        return self.pores(labels=labels,mode=mode)

    def get_throat_indices(self,labels=['all'],mode='union'):
        r'''
        THIS METHOD IS DEPRECATED, USE throats() INSTEAD
        '''
        return self.throats(labels=labels,mode=mode)
        
    def tomask(self,pores=None,throats=None):
        r'''
        Convert a list of pore or throat indices into a boolean mask
        
        Parameters
        ----------
        pores or throats : array_like
            List of pore or throat indices
            
        Returns
        -------
        mask : array_like
            A boolean mask of length Np or Nt with True in the locations of
            pores or throats received.  
        
        '''
        if pores != None:
            Np = sp.shape(self['pore.all'])[0]
            pores = sp.array(pores,ndmin=1)
            mask = sp.zeros((Np,),dtype=bool)
            mask[pores] = True
            return mask
        if throats != None:
            Nt = sp.shape(self['throat.all'])[0]
            throats = sp.array(throats,ndmin=1)
            mask = sp.zeros((Nt,),dtype=bool)
            mask[throats] = True
            return mask

    def interpolate_data(self,data=[],prop='',throats=[],pores=[]):
        r"""
        Determines a pore (or throat) property as the average of it's neighboring 
        throats (or pores)

        Parameters
        ----------
        pores : array_like
            The pores for which values are desired
        throats : array_like
            The throats for which values are desired
        data : array_like
            A list of specific values to be interploated.  List MUST be either
            Np or Nt long
        
        Returns
        -------
        values : array_like
            An array containing interpolated pore (or throat) data

        Notes
        -----
        - This uses an unweighted average, without attempting to account for 
        distances or sizes of pores and throats.
        - Only one of pores, throats OR data are accepted

        """
        if self.__module__.split('.')[1] == 'Network': net = self
        else: net = self._net
        if sp.shape(data)[0] > 0:
            prop='tmp'
            values = sp.array(data,ndmin=1)
            if sp.shape(values)[0] == net.num_pores():
                throats = net.get_throat_indices('all')
                self.set_data(prop=prop,pores='all',data=values)
            elif sp.shape(values)[0] == net.num_throats():
                pores = self.get_pore_indices('all')
                self.set_data(prop=prop,throats='all',data=values)
            elif sp.shape(values)[0] == 1: #if scalar was sent
                return values
            else:
                self._logger.error('Received data of an ambiguous length')
                return
        if sp.shape(pores)[0] > 0:
            throats = net.find_neighbor_throats(pores,flatten=False)
            throat_data = self.get_data(prop=prop,throats='all')
            values = sp.ones((sp.shape(pores)[0],))*sp.nan
            if sp.shape(throat_data)[0] == 1:
                values = throat_data
            else:
                for i in pores:
                    values[i] = sp.mean(throat_data[throats[i]])
        elif sp.shape(throats)[0] > 0:
            pores = net.find_connected_pores(throats,flatten=False)
            pore_data = self.get_data(prop=prop,pores='all')
            values = sp.ones((sp.shape(throats)[0],))*sp.nan
            if sp.shape(pore_data)[0] == 1:
                values = pore_data
            else:
                values = sp.mean(pore_data[pores],axis=1)
        return values
        
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
            
            * 'difference' : Only pores with none of the given labels are counted.
            
        Returns
        -------
        Np : int
            Number of pores with the specified labels 
            
        See Also
        --------
        num_throats, count
            
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
        Np = self.pores(labels=labels,mode=mode)
        Np = self.tomask(pores=Np)
        return sp.sum(Np) #return sum of Trues
            
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
            
            * 'difference' : Only throats with none of the given labels are counted.

        Returns
        -------
        Nt : int
            Number of throats with the specified labels
            
        See Also
        --------
        num_pores, count

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
        Nt = self.throats(labels=labels,mode=mode)
        Nt = self.tomask(throats=Nt)
        return sp.sum(Nt) #return sum of Trues
        
    def count(self,element=None):
        r'''
        Returns a dictionary containing the number of pores and throats in 
        the network, stored under the keys 'pore' or 'throat'
        
        Parameters
        ----------
        element : string, optional
            Can be either 'pore' or 'throat', which specifies which count to return.
            
        Returns
        -------
        A dictionary containing the number of pores and throats under the 
        'pore' and 'throat' key respectively.  
        
        See Also
        --------
        num_pores, num_throats
            
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.count()
        {'pore': 125, 'throat': 300}
        >>> pn.count('pore')
        125
        '''
        temp = {}
        temp['pore'] = self.num_pores()
        temp['throat'] = self.num_throats()
        if element != None:
            temp = temp[element]
        return temp
        
    def data_health(self,element='',props=[],quiet=False):
        r'''
        Check the health of pore and throat data arrays.  
        
        Parameters
        ----------
        element : string, optional
            Can be either 'pore' or 'throat', which will limit the checks to 
            only those data arrays.
        props : list of pore (or throat) properties, optional
            If given, will limit the health checks to only the specfied
            properties.  Also useful for checking existance.
        quiet : bool, optional
            By default this method will output a summary of the health check.
            This can be disabled by setting quiet to False.
            
        Returns
        -------
        Returns a True if all check pass, and False if any checks fail.  This
        is ideal for programatically checking data integrity prior to running
        an algorithm.
        
        '''
        health = {}
        flag = True
        if props == []:
            props = self.props(element)
        else:
            if type(props) == str:
                props = [props]
            if props[0].split('.')[0] not in ['pore','throat']:
                self._logger.error('Properties must be either pore or throat')
        for item in props:
            try: 
                if sp.sum(sp.isnan(self[item])) > 0:
                    health[item] = 'Has NaNs'
                    flag = False
                elif sp.shape(self[item])[0] == 1:
                    health[item] = 'Healthy Scalar'
                elif sp.shape(self[item])[0] == self.count(item.split('.')[0]):
                    health[item] = 'Healthy Vector'
                else:
                    health[item] = 'Wrong Length'
                    flag = False
            except: 
                health[item] = 'Does not exist'
                flag = False
        if quiet == False:
            pprint.pprint(health)
        return flag
            
    def check_pore_health(self,props=[],quiet=False):
        r'''
        This method is deprecated, use data_health instead
        '''
        if props != []:
            if type(props) == str:
                props = [props]
            temp = []
            for item in props:
                item = item.split('.')[-1]
                temp.append('pore.' + item)
            props = temp
        return self.data_health(element='pore',props=props,quiet=quiet)            
        
    def check_throat_health(self,props=[],quiet=False):
        r'''
        This method is deprecated, use data_health instead
        '''
        if props != []:
            if type(props) == str:
                props = [props]
            temp = []
            for item in props:
                item = item.split('.')[-1]
                temp.append('throat.' + item)
            props = temp
        return self.data_health(element='throat',props=props,quiet=quiet)
        

if __name__ == '__main__':
    import doctest
    doctest.testmod(verbose=True)

