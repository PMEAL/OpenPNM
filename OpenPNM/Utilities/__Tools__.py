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
        self._logger.info("Construction of Object container complete")
        
    #--------------------------------------------------------------------------
    '''Setter and Getter Methods'''
    #--------------------------------------------------------------------------
    def _set_data(self,element='',phase='',prop='',data='',locations=''):
        r'''
        Documentation for this method is being updated, we are sorry for the inconvenience.
        '''
        if type(data)!=sp.ndarray: data = sp.array(data,ndmin=1)
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
              
        try: phase = self.find_object_by_name(phase) 
        except: pass #Accept object

        if phase :

            if sp.shape(data)[0]==1:
                if locations!='':                
                    try: getattr(phase,'_'+element+'_data')[prop]
                    except: getattr(phase,'_'+element+'_data')[prop] = sp.zeros((getattr(phase,'num_'+element+'s')(),))*sp.nan
                    getattr(phase,'_'+element+'_data')[prop][locations] = data
                else:
                    try: 
                        getattr(phase,'_'+element+'_data')[prop]
                        if sp.shape(getattr(phase,'_'+element+'_data')[prop])[0]!=1:
                            print('Warning: '+prop+' '+element+' property in '+phase.name+' was an array which has been overwritten with a scalar value')
                    except: pass
                    getattr(phase,'_'+element+'_data')[prop] = data  
            else:                
                if locations!='':
                    if sp.shape(locations)[0]==sp.shape(data)[0]:
                        try: getattr(phase,'_'+element+'_data')[prop]
                        except: getattr(phase,'_'+element+'_data')[prop] = sp.zeros((getattr(phase,'num_'+element+'s')(),))*sp.nan
                        getattr(phase,'_'+element+'_data')[prop][locations] = data
                    else: 
                        phase._logger.error('For adding '+element+' property '+prop+' to '+phase.name+', locations and size of data do not match!')
                else:
                    try: 
                        getattr(phase,'num_'+element+'s')()                        
                        if sp.shape(data)[0]==getattr(phase,'num_'+element+'s')():
                            getattr(phase,'_'+element+'_data')[prop] = data
                        else: phase._logger.error('For adding '+element+' property '+prop+' to '+phase.name+', number of '+element+'s and size of data do not match!')
                    except: phase._logger.error(element+' numbering has not been specified for '+phase.name)
            phase._logger.debug(element+' property '+prop+' has been added to '+phase.name)
                        
        else:
            
            if sp.shape(data)[0]==1:
                if locations!='':                
                    try: getattr(self,'_'+element+'_data')[prop]
                    except: getattr(self,'_'+element+'_data')[prop] = sp.zeros((getattr(self,'num_'+element+'s')(),))*sp.nan
                    getattr(self,'_'+element+'_data')[prop][locations] = data
                else:
                    try: 
                        getattr(self,'_'+element+'_data')[prop]
                        if sp.shape(getattr(self,'_'+element+'_data')[prop])[0]!=1:
                            print('Warning: '+prop+' '+element+' property in '+self.name+' was an array which has been overwritten with a scalar value')
                    except: pass
                    getattr(self,'_'+element+'_data')[prop] = data            
            else:                
                if locations!='':
                    if sp.shape(locations)[0]==sp.shape(data)[0]:
                        try: getattr(self,'_'+element+'_data')[prop]
                        except: getattr(self,'_'+element+'_data')[prop] = sp.zeros((getattr(self,'num_'+element+'s')(),))*sp.nan
                        getattr(self,'_'+element+'_data')[prop][locations] = data
                    else: self._logger.error('For adding '+element+' property '+prop+' to '+self.name+', locations and size of data do not match!')
                else:
                    try: 
                        getattr(self,'num_'+element+'s')()                        
                        if sp.shape(data)[0]==getattr(self,'num_'+element+'s')():
                            getattr(self,'_'+element+'_data')[prop] = data
                        else: self._logger.error('For adding '+element+' property '+prop+' to '+self.name+', number of '+element+'s and size of data do not match!')
                    except: getattr(self,'_'+element+'_data')[prop] = data
            self._logger.debug(element+' property '+prop+' has been added to '+self.name)



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
        try: phase = self.find_object_by_name(phase) 
        except: pass        
        
        if phase :
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
        This is wrapper method that calls set_data, which is generic for pores and throats
        
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
        This is a wrapper method that calls get_data, which is generic for pores and throats
        
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
          
            
            if label:
                if label=='all': getattr(self,'_'+element+'_info')[label] = sp.ones_like(locations,dtype=bool)
                else:    
                    try: getattr(self,'_'+element+'_info')[label]
                    except: getattr(self,'_'+element+'_info')[label] = sp.zeros((getattr(self,'num_'+element+'s')(),),dtype=bool)
                    if mode=='overwrite': getattr(self,'_'+element+'_info')[label] = sp.zeros((getattr(self,'num_'+element+'s')(),),dtype=bool)
                    getattr(self,'_'+element+'_info')[label][locations] = True
            else: self._logger.error('No label has been defined for these locations')                

        else: getattr(self,'_'+element+'_info')[label] = sp.zeros((getattr(self,'num_'+element+'s')(),),dtype=bool)

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
        Parameters
        ----------
        label : string
            The name of the pore labels you wish to apply (e.g. 'top')
        locaitons : array_like
            An array containing the locations (pores) where the labels should be applied.
            Can be either a boolean mask of Np length with True at labels locations (default), 
            a list of indices where labels should be applied. 
        mode : string
            Options are 'merge' and 'overwrite', default is 'merge'
        is_indices : boolean
            This flag indicates whether locations are being sent as a boolean maks (default), 
            or a list of indices.
            
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
        self._set_info(element='pore',label=label,locations=locations)

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
        self._set_info(element='throat',label=label,locations=locations)
        
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
        Returns a list containing the names of all defined pore properties
        '''
        return list(self._pore_data.keys())

    def list_throat_props(self):
        r'''
        Returns a list containing the names of all defined throat properties
        '''
        return list(self._throat_data.keys())
        
    def list_pore_labels(self):
        r'''
        Returns a list containing the names of all defined pore labels
        '''
        return list(self._pore_info.keys())

    def list_throat_labels(self):
        r'''
        Returns a list containing the names of all defined throat labels
        '''
        return list(self._throat_info.keys())
        
    def find_labels(self,pnum='',tnum=''):
        r'''
        Returns a list of all labels that have been applied to the given pore or throat.

        Parameters
        ----------
        pnum : int
            The pore who's labels are sought
        tnum : int
            The throat who's labels are sought
            
        Returns
        -------
        A list of labels that have been applied to the given pore or throat
        
        Notes
        -----
        This only accepts a single pore or throat value for now, and works with
        brute force appraoch.  Vectorization would be nice.
        
        Also, the logic for allowing pnum/tnum to be either int or list is clunky (but works)
        
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.find_labels(pnum=0)
        ['all', 'bottom', 'front', 'internal', 'left']
        >>> pn.find_labels(tnum=124)
        ['all', 'internal', 'right']
        '''
        if pnum != '' and tnum == '':
            element = 'pore'
            num = pnum
        elif tnum != '' and pnum == '':
            element = 'throat'
            num = tnum
        else: self._logger.error(sys._getframe().f_code.co_name+' can only accept one of tnum or pnum')
        labels = []
        for item in getattr(self,'_'+element+'_info').keys():
            if getattr(self,'_'+element+'_info')[item][num]:
                labels.append(item)
        labels.sort()
        return labels
        
    def has_labels(self,pnums='',tnums='',labels='all',mode='union'):
        r'''
        This method accepts a list of pores (or throats) and a list of labels, 
        and returns True if pore or throat has any (or all) of the labels.
        
        Parameters
        ----------
        pnums : int or list of ints
            Pore numbers (locations) of interest
        tnum : int or list of ints
            Throat numbers (locations) of interest
        labels : str or list of strings
            Labels of interest, defaults to 'all'
        mode : str, optional
            Flag to control logic, options are 'union' (default) or 'intersection'
            
        Returns
        -------
        A boolean list the same shape as pnums (or tnums) containing truth values
        indicating whether the corresponding pore (or throat) has the specfied labels
        
        Notes
        -----
        The logic for allowing pnum/tnum to be either int or list is clunky (but works)
        
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.has_labels(pnums=[0,1,2,122,123,124],labels=['top','front'],mode='union') #Default mode is 'union'
        array([ True, False, False,  True,  True,  True], dtype=bool)
        >>> pn.has_labels(pnums=[0,1,2,122,123,124],labels=['top','right'],mode='intersection')
        array([False, False, False,  True,  True,  True], dtype=bool)
        '''
        #Parse input arguments
        if pnums != '' and tnums == '':
            element = 'pore'
            nums = pnums
        elif tnums != '' and pnums == '':
            element = 'throat'
            nums = tnums
        else: self._logger.error(sys._getframe().f_code.co_name+' can only accept one of tnum or pnum')
        if type(labels) == str: labels = [labels] #convert string to list, if necessary
        mask = getattr(self,'get_'+element+'_indices')(labels=labels,indices=False,mode=mode)
        return mask[nums]

    def check_info(self):
        r'''
        Documentation for this method is being updated, we are sorry for the inconvenience.
        '''
        temp = sp.zeros_like(self.get_pore_data(prop='coords')[:,0],dtype=bool)
        self.set_pore_info(label='all',locations=temp)
        for item in self._pore_info.keys():
            if sp.shape(self._pore_info[item])[0] != sp.shape(self._pore_info['all'])[0]:
                print('warning, info arrays are wrong size!')
        temp = sp.zeros_like(self.get_throat_data(prop='connections')[:,0],dtype=bool)
        self.set_throat_info(label='all',locations=temp)
        for item in self._throat_info.keys():
            if sp.shape(self._throat_info[item])[0] != sp.shape(self._throat_info['all'])[0]:
                print('warning, info arrays are wrong size!')

    #--------------------------------------------------------------------------
    '''Object query methods'''
    #--------------------------------------------------------------------------
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
            Specifies whether the count should be done as a union (default) or intersection.
            In union mode, all pores with ANY of the given labels are counted.
            In intersection mode, only pores with ALL the give labels are counted.
            
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
        
        '''
        #convert string to list, if necessary
        if type(labels) == str: labels = [labels]
        #Count number of pores of specified type
        temp = self.get_pore_indices(labels=labels,mode=mode,indices=False)
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
            Specifies whether the count should be done as a union (default) or intersection.
            In union mode, all throats with ANY of the given labels are counted.
            In intersection mode, only throats with ALL the give labels are counted.

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
        
        '''
        #convert string to list, if necessary
        if type(labels) == str: labels = [labels]
        #Count number of pores of specified type
        temp = self.get_throat_indices(labels=labels,mode=mode,indices=False)
        return sp.sum(temp) #return sum of Trues

    def get_pore_indices(self,labels=['all'],indices=True,mode='union'):
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
            Specifies whether the indices should be union (default) or intersection of desired labels.
            In union mode, all pores with ANY of the given labels are included.
            In intersection mode, only pores with ALL the give labels are included.
            This is ignored if no label is given.
        
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.get_pore_indices(labels=['top','front'],mode='intersection')
        array([100, 105, 110, 115, 120], dtype=int64)
        '''
        if type(labels) == str: labels = [labels] #convert string to list, if necessary
        if mode == 'union':
            union = sp.zeros_like(self.get_pore_info(label='all'),dtype=bool)
            for item in labels: #iterate over labels list and collect all indices
                    union = union + self._get_info(element='pore',label=item)
            ind = union
        elif mode == 'intersection':
            intersect = sp.ones((self.num_pores(),),dtype=bool)
            for item in labels: #iterate over labels list and collect all indices
                    intersect = intersect*self._get_info(element='pore',label=item)
            ind = intersect
        if indices: ind = sp.where(ind==True)[0]
        return ind

    def get_throat_indices(self,labels=['all'],indices=True,mode='union'):
        r'''
        Returns throat locations where given labels exist.
        
        Parameters
        ----------
        labels : list of strings, optional
            The throat label(s) whose locations are requested.
            If omitted, all throat inidices are returned.
        indices : boolean, optional
            This flag specifies whether throat locations are returned as a boolean mask of length Np,
            or as a list of indices (default).
        mode : string, optional
            Specifies whether the indices should be union (default) or intersection of desired labels.
            In union mode, all throats with ANY of the given labels are included.
            In intersection mode, only throats with ALL the give labels are included.
            This is ignored if no labels are given.
        
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> Tind = pn.get_throat_indices()
        >>> Tind[0:5]
        array([0, 1, 2, 3, 4], dtype=int64)
        '''
        if type(labels) == str: labels = [labels] #convert string to list, if necessary
        if mode == 'union':
            union = sp.zeros_like(self.get_throat_info(label='all'),dtype=bool)
            for item in labels: #iterate over labels list and collect all indices
                    union = union + self._get_info(element='throat',label=item)
            ind = union
        elif mode == 'intersection':
            intersect = sp.ones((self.num_throats(),),dtype=bool)
            for item in labels: #iterate over labels list and collect all indices
                    intersect = intersect*self._get_info(element='throat',label=item)
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

