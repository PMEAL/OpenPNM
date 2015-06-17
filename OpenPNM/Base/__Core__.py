"""
###############################################################################
Core:  Core Data Class
###############################################################################
"""
import string, random
import scipy as sp
import scipy.constants
from OpenPNM.Base import logging, Tools
from OpenPNM.Base import ModelsDict
logger = logging.getLogger()
from OpenPNM.Base import Controller
ctrl = Controller()

class Core(dict):
    r"""
    Contains OpenPNM specificmethods for working with the data in the dictionaries
    """

    def __new__(typ, *args, **kwargs):
        obj = dict.__new__(typ, *args, **kwargs)
        obj.update({'pore.all': sp.array([], ndmin=1, dtype=bool)})
        obj.update({'throat.all': sp.array([], ndmin=1, dtype=bool)})
        # Initialize phase, physics, and geometry tracking lists
        obj._name = None
        obj._ctrl = {}
        obj._phases = []
        obj._geometries = []
        obj._physics = []
        obj._net = None
        obj._parent = None
        # Initialize ordered dict for storing property models
        obj.models = ModelsDict()
        return obj

    def __init__(self, name=None, **kwargs):
        super().__init__()
        logger.debug('Initializing Core class')
        self.name = name
        self.controller = ctrl

    def __repr__(self):
        return '<%s.%s object at %s>' % (
            self.__class__.__module__,
            self.__class__.__name__,
            hex(id(self)))

    def __eq__(self,other):
        if hex(id(self)) == hex(id(other)):
            return True
        else:
            return False

    def __setitem__(self,key,value):
        r"""
        This is a subclass of the default __setitem__ behavior.  The main aim
        is to limit what type and shape of data can be written to protect
        the integrity of the network.


        Example
        -------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn['pore.example_property'] = 100
        >>> pn['pore.example_property'][0]
        100

        """
        # Enforce correct dict naming
        element = key.split('.')[0]
        if (element != 'pore') and (element != 'throat'):
            logger.error('Array name \''+key+'\' does not begin with \'pore\' or \'throat\'')
            return
        # Convert value to an ndarray
        value = sp.array(value,ndmin=1)
        #Skip checks for 'coords', 'conns'
        if (key == 'pore.coords') or (key == 'throat.conns'):
            super(Core, self).__setitem__(key,value)
            return
        # Skip checks for protected props, and prevent changes if defined
        if key.split('.')[1] in ['all']:
            if key in self.keys():
                if sp.shape(self[key]) == (0,):
                    logger.debug(key+' is being defined.')
                    super(Core, self).__setitem__(key,value)
                else:
                    logger.warning(key+' is already defined.')
                    return
            else:
                logger.debug(key+' is being defined.')
                super(Core, self).__setitem__(key,value)
            return
        # Write value to dictionary
        if sp.shape(value)[0] == 1:  # If value is scalar
            logger.debug('Broadcasting scalar value into vector: '+key)
            value = sp.ones((self._count(element),),dtype=value.dtype)*value
            super(Core, self).__setitem__(key,value)
        elif sp.shape(value)[0] == self._count(element):
            logger.debug('Updating vector: '+key)
            super(Core, self).__setitem__(key,value)
        else:
            if self._count(element) == 0:
                self.update({key:value})
            else:
                logger.warning('Cannot write vector with an array of the wrong length: '+key)
                pass

    def _set_ctrl(self,controller):
        if self.name in controller.keys():
            raise Exception('An object named '+self.name+' is already present in simulation')
        self._ctrl = controller
        controller.update({self.name: self})

    def _get_ctrl(self):
        return self._ctrl

    controller = property(_get_ctrl,_set_ctrl)

    def _set_name(self,name):
        if name in self.controller.keys():
            raise Exception('An object named '+name+' already exists')
        elif name is None:
            name = ''.join(random.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for _ in range(5))
            name = self.__module__.split('.')[-1].strip('__') + '_' + name
        elif self._name is not None:
            logger.info('Changing the name of '+self.name+' to '+name)
            # Check if name collides with any arrays in the simulation
            objs = self._simulation()
            for item in objs:
                keys = [key.split('.')[-1] for key in item.keys()]
                if name in keys:
                    raise Exception(name+' is already in use as an array name')
            for item in objs:
                if 'pore.'+self.name in item.keys():
                    item['pore.'+name] = item.pop('pore.'+self.name)
                if 'throat.'+self.name in item.keys():
                    item['throat.'+name] = item.pop('throat.'+self.name)
            self._ctrl[name] = self._ctrl.pop(self.name)
        self._name = name

    def _get_name(self):
        return self._name

    name = property(_get_name,_set_name)

    def _simulation(self):
        temp = []
        temp += [self._net]
        temp += self._net._phases
        temp += self._net._geometries
        temp += self._net._physics
        return temp

    def clear(self):
        r"""
        A subclassed version of the standard dict's clear method
        """
        pall = self['pore.all']
        tall = self['throat.all']
        super(Core,self).clear()
        self.update({'throat.all':tall})
        self.update({'pore.all':pall})

    #--------------------------------------------------------------------------
    """Model Manipulation Methods"""
    #--------------------------------------------------------------------------
    #Note: These methods have been moved to the ModelsDict class but are left
    #here for backward compatibility
    def add_model(self,propname,model,regen_mode='normal',**kwargs):
        self.models.add(propname=propname,model=model,regen_mode=regen_mode,**kwargs)

    add_model.__doc__ = ModelsDict.add.__doc__

    def regenerate(self,props='',mode='inclusive'):
        self.models.regenerate(props=props,mode=mode)

    regenerate.__doc__ = ModelsDict.regenerate.__doc__

    #--------------------------------------------------------------------------
    'Object lookup methods'
    #--------------------------------------------------------------------------

    def _find_object(self,obj_name='',obj_type=''):
        r"""
        Find objects associated with a given network model by name or type

        Parameters
        ----------
        obj_name : string
           Name of sought object

        obj_type : string
            The type of object beign sought.  Options are:

            1. 'Network' or 'Networks'
            2. 'Geometry' or 'Geometries'
            3. 'Phase' or 'Phases'
            4. 'Physics'

        Returns
        -------
        OpenPNM object or list of objects

        """
        if obj_name != '':
            obj = []
            if obj_name in ctrl.keys():
                obj = ctrl[obj_name]
            return obj
        elif obj_type != '':
            if obj_type in ['Geometry','Geometries','geometry','geometries']:
                objs = ctrl.geometries()
            elif obj_type in ['Phase','Phases','phase','phases']:
                objs = ctrl.phases()
            elif obj_type in ['Physics','physics']:
                objs = ctrl.physics()
            elif obj_type in ['Network','Networks','network','networks']:
                objs = ctrl.networks()
            return objs

    def physics(self,phys_name=[]):
        r"""
        Retrieves Physics associated with the object

        Parameters
        ----------
        name : string or list of strings, optional
            The name(s) of the Physics object to retrieve
        Returns
        -------
            If name is NOT provided, then a list of Physics names is returned.
            If a name or list of names IS provided, then the Physics object(s)
            with those name(s) is returned.
        """
        # If arg given as string, convert to list
        if type(phys_name) == str:
            phys_name = [phys_name]
        if phys_name == []:  # If default argument received
            phys = [item.name for item in self._physics]
        else:  # If list of names received
            phys = []
            for item in self._physics:
                if item.name in phys_name:
                    phys.append(item)
        return phys

    def phases(self,phase_name=[]):
        r"""
        Retrieves Phases associated with the object

        Parameters
        ----------
        name : string or list of strings, optional
            The name(s) of the Phase object(s) to retrieve.
        Returns
        -------
            If name is NOT provided, then a list of phase names is returned. If
            a name are provided, then a list containing the requested objects
            is returned.
        """
        # If arg given as string, convert to list
        if type(phase_name) == str:
            phase_name = [phase_name]
        if phase_name == []:  # If default argument received
            phase = [item.name for item in self._phases]
        else:  # If list of names received
            phase = []
            for item in self._phases:
                if item.name in phase_name:
                    phase.append(item)
        return phase

    def geometries(self,geom_name=[]):
        r"""
        Retrieves Geometry object(s) associated with the object

        Parameters
        ----------
        name : string or list of strings, optional
            The name(s) of the Geometry object to retrieve.
        Returns
        -------
            If name is NOT provided, then a list of Geometry names is returned.
            If a name IS provided, then the Geometry object of that name is
            returned.
        """
        # If arg given as string, convert to list
        if type(geom_name) == str:
            geom_name = [geom_name]
        if geom_name == []:  # If default argument received
            geom = [item.name for item in self._geometries]
        else:  # If list of names received
            geom = []
            for item in self._geometries:
                if item.name in geom_name:
                    geom.append(item)
        return geom

    def network(self,name=''):
        r"""
        Retrieves the network associated with the object.  If the object is
        a network, then it returns a handle to itself.

        Parameters
        ----------
        name : string, optional
            The name of the Network object to retrieve.

        Returns
        -------
            If a name IS provided, then the parent netowrk object is returned.

        Notes
        -----
        This doesn't quite work yet...we have to decide how to treat sub-nets first
        """
        if name == '':
            if self._net is None:
                net = [self]
            else:
                net = [self._net]
        else:
            net = []
            temp = self._find_object(obj_name=name)
            if hasattr(temp,'_isa'):
                if temp._isa('Network'):
                    net = temp
        return net

    #--------------------------------------------------------------------------
    """Data Query Methods"""
    #--------------------------------------------------------------------------
    def props(self,element='',mode='all'):
        r"""
        Returns a list containing the names of all defined pore or throat
        properties.

        Parameters
        ----------
        element : string, optional
            Can be either 'pore' or 'throat' to specify what properties are
            returned.  If no element is given, both are returned

        mode : string, optional
            Controls what type of properties are returned.  Options are:

            - 'all' : Returns all properties on the object
            - 'models' : Returns only properties that are associated with a model
            - 'constants' : Returns only properties that are set as constant values

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
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.props('pore')
        ['pore.coords']
        >>> pn.props('throat')
        ['throat.conns']
        >>> #pn.props() # this lists both, but in random order, which breaks
        >>> #           # our automatic document testing so it's commented here
        """

        props = []
        for item in self.keys():
            if self[item].dtype != bool:
                props.append(item)

        all_models = list(self.models.keys())
        constants = [item for item in props if item not in all_models]
        models = [item for item in props if item in all_models]

        if element in ['pore','pores']:
            element = 'pore'
        elif element in ['throat','throats']:
            element = 'throat'

        temp = []
        if mode == 'all':
            if element == '': temp = props
            else: temp = [item for item in props if item.split('.')[0]==element]
        elif mode == 'models':
            if element == '': temp = models
            else: temp = [item for item in models if item.split('.')[0]==element]
        elif mode == 'constants':
            if element == '': temp = constants
            else: temp = [item for item in constants if item.split('.')[0]==element]
        a = Tools.PrintableList(temp)
        return a


    def _get_labels(self,element='',locations=[],mode='union'):
        r"""
        This is the actual label getter method, but it should not be called directly.
        Wrapper methods have been created, use labels().
        """
        # Collect list of all pore OR throat labels
        labels = []
        for item in self.keys():
            if item.split('.')[0] == element:
                if self[item].dtype in ['bool']:
                    labels.append(item)
        labels.sort()
        if sp.size(locations) == 0:
            return Tools.PrintableList(labels)
        else:
            labels = sp.array(labels)
            locations = sp.array(locations,ndmin=1)
            if locations.dtype in ['bool']:
                locations = self._get_indices(element=element)[locations]
            else:
                locations = sp.array(locations,dtype=int)
            arr = sp.zeros((sp.shape(locations)[0],len(labels)),dtype=bool)
            col = 0
            for item in labels:
                arr[:,col] = self[item][locations]
                col = col + 1
            if mode == 'count':
                return sp.sum(arr,axis=1)
            if mode == 'union':
                temp = labels[sp.sum(arr,axis=0)>0]
                temp.tolist()
                return Tools.PrintableList(temp)
            if mode == 'intersection':
                temp = labels[sp.sum(arr,axis=0)==sp.shape(locations,)[0]]
                temp.tolist()
                return Tools.PrintableList(temp)
            if mode in ['difference', 'not']:
                temp = labels[sp.sum(arr,axis=0)!=sp.shape(locations,)[0]]
                temp.tolist()
                return Tools.PrintableList(temp)
            if mode == 'mask':
                return arr
            if mode == 'none':
                temp = sp.ndarray((sp.shape(locations,)[0],),dtype=object)
                for i in sp.arange(0,sp.shape(locations,)[0]):
                    temp[i] = list(labels[arr[i,:]])
                return temp
            else:
                logger.error('unrecognized mode:'+mode)

    def labels(self,element='',pores=[],throats=[],mode='union'):
        r"""
        Returns the labels applied to specified pore or throat locations

        Parameters
        ----------
        pores (or throats) : array_like
            The pores (or throats) whose labels are sought.  If left empty a
            list containing all pore and throat labels is returned.

        element : string
            Controls whether pore or throat labels are returned.  If empty then
            both are returned.

        mode : string, optional
            Controls how the query should be performed

            * 'none' : An N x Li list of all labels applied to each input pore (or throats). Li can vary betwen pores (and throats)

            * 'union' : A list of labels applied to ANY of the given pores (or throats)

            * 'intersection' : Label applied to ALL of the given pores (or throats)

            * 'not' : Labels NOT applied to ALL pores (or throats)

            * 'count' : The number of labels on each pores (or throats)

            * 'mask' : returns an N x Lt array, where each row corresponds to a pore (or throat) location, and each column contains the truth value for the existance of labels as returned from labels(pores='all',mode='union')).

        Examples
        --------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.labels(pores=[0,1,5,6])
        ['pore.all', 'pore.bottom', 'pore.front', 'pore.left']
        >>> pn.labels(pores=[0,1,5,6],mode='intersection')
        ['pore.all', 'pore.bottom']
        """
        if (sp.size(pores) == 0) and (sp.size(throats) == 0):
            if element == '':
                temp = []
                temp = self._get_labels(element='pore')
                temp.extend(self._get_labels(element='throat'))
            elif element in ['pore','pores']:
                temp = self._get_labels(element='pore',locations=[], mode=mode)
            elif element in ['throat','throats']:
                temp = self._get_labels(element='throat',locations=[], mode=mode)
            else:
                logger.error('Unrecognized element')
                return
        elif sp.size(pores) != 0:
            if pores == 'all':
                pores = self.pores()
            pores = sp.array(pores,ndmin=1)
            temp = self._get_labels(element='pore',locations=pores, mode=mode)
        elif sp.size(throats) != 0:
            if throats == 'all':
                throats = self.throats()
            throats = sp.array(throats,ndmin=1)
            temp = self._get_labels(element='throat',locations=throats,mode=mode)
        return temp

    def filter_by_label(self,pores=[],throats=[],labels='',mode='union'):
        r"""
        Returns which of the supplied pores (or throats) has the specified label

        Parameters
        ----------
        pores, or throats : array_like
            List of pores or throats to be filtered

        labels : list of strings
            The labels to apply as a filter

        mode : string
            Controls how the filter is applied.  Options include:

            * 'union' : (default) All locations with ANY of the given labels are kept.

            * 'intersection' : Only locations with ALL the given labels are kept.

            * 'not_intersection' : Only locations with exactly one of the given labels are kept.

            * 'not' : Only locations with none of the given labels are kept.

        See Also
        --------
        pores
        throats

        Examples
        --------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.filter_by_label(pores=[0,1,5,6], labels='left')
        array([0, 1])
        >>> Ps = pn.pores(['top', 'bottom', 'front'], mode='union')
        >>> pn.filter_by_label(pores=Ps, labels=['top', 'front'], mode='intersection')
        array([100, 105, 110, 115, 120])
        """
        if labels == '':  # Handle empty labels
            labels = 'all'
        if type(labels) == str:  # Convert input to list
            labels = [labels]
        # Convert inputs to locations and element
        if sp.size(pores) > 0:
            element = 'pore'
            locations = sp.array(pores)
        if sp.size(throats) > 0:
            element = 'throat'
            locations = sp.array(throats)
        # Do it
        labels = [element+'.'+item.split('.')[-1] for item in labels]
        all_locs = self._get_indices(element=element, labels=labels, mode=mode)
        mask = self._tomask(locations=all_locs, element=element)
        ind = mask[locations]
        return locations[ind]

    def _get_indices(self,element,labels=['all'],mode='union'):
        r"""
        This is the actual method for getting indices, but should not be called
        directly.  Use pores or throats instead.
        """
        element.rstrip('s')  # Correct plural form of element keyword
        if element+'.all' not in self.keys():
            raise Exception('Cannot proceed without {}.all'.format(element))
        if type(labels) == str:  # Convert string to list, if necessary
            labels = [labels]
        for label in labels:  # Parse the labels list for wildcards "*"
            if label.startswith('*'):
                labels.remove(label)
                temp = [item for item in self.labels() if item.split('.')[-1].endswith(label.strip('*'))]
                if temp == []:
                    temp = [label.strip('*')]
                labels.extend(temp)
            if label.endswith('*'):
                labels.remove(label)
                temp = [item for item in self.labels() if item.split('.')[-1].startswith(label.strip('*'))]
                if temp == []:
                    temp = [label.strip('*')]
                labels.extend(temp)
        # Begin computing label array
        if mode == 'union':
            union = sp.zeros_like(self[element+'.all'],dtype=bool)
            for item in labels: #iterate over labels list and collect all indices
                    union = union + self[element+'.'+item.split('.')[-1]]
            ind = union
        elif mode == 'intersection':
            intersect = sp.ones_like(self[element+'.all'],dtype=bool)
            for item in labels: #iterate over labels list and collect all indices
                    intersect = intersect*self[element+'.'+item.split('.')[-1]]
            ind = intersect
        elif mode == 'not_intersection':
            not_intersect = sp.zeros_like(self[element+'.all'],dtype=int)
            for item in labels: #iterate over labels list and collect all indices
                info = self[element+'.'+item.split('.')[-1]]
                not_intersect = not_intersect + sp.int8(info)
            ind = (not_intersect == 1)
        elif mode in ['difference','not']:
            none = sp.zeros_like(self[element+'.all'],dtype=int)
            for item in labels: #iterate over labels list and collect all indices
                info = self[element+'.'+item.split('.')[-1]]
                none = none - sp.int8(info)
            ind = (none == 0)
        #Extract indices from boolean mask
        ind = sp.where(ind==True)[0]
        ind = ind.astype(dtype=int)
        return ind

    def pores(self,labels='all',mode='union'):
        r"""
        Returns pore locations where given labels exist.

        Parameters
        ----------
        labels : list of strings, optional
            The pore label(s) whose locations are requested.  If omitted, all
            pore inidices are returned. This argument also accepts '*' for
            wildcard searches.
        mode : string, optional
            Specifies how the query should be performed.  The options are:

            * 'union' : (default) All pores with ANY of the given labels are returned.

            * 'intersection' : Only pore with ALL the given labels are returned.

            * 'not_intersection' : Only pores with exactly one of the given labels are returned.

            * 'not' : Only pores with none of the given labels are returned.

        Examples
        --------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> pind = pn.pores(labels=['top','front'],mode='union')
        >>> pind[[0,1,2,-3,-2,-1]]
        array([  0,   5,  10, 122, 123, 124])
        >>> pn.pores(labels=['top','front'],mode='intersection')
        array([100, 105, 110, 115, 120])
        """
        if labels == 'all':
            Np = sp.shape(self['pore.all'])[0]
            ind = sp.arange(0,Np)
        else:
            ind = self._get_indices(element='pore',labels=labels,mode=mode)
        return ind

    @property
    def Ps(self):
        r"""
        A shortcut to get a list of all pores on the object
        """
        return self.pores()

    def throats(self,labels='all',mode='union'):
        r"""
        Returns throat locations where given labels exist.

        Parameters
        ----------
        labels : list of strings, optional
            The throat label(s) whose locations are requested.  If omitted,
            'all' throat inidices are returned.  This argument also accepts
            '*' for wildcard searches.
        mode : string, optional
            Specifies how the query should be performed.  The options are:

            * 'union' : (default) All throats with ANY of the given labels are returned.

            * 'intersection' : Only throats with ALL the given labels are counted.

            * 'not_intersection' : Only throats with exactly one of the given labels are counted.

            * 'not' : Only throats with none of the given labels are returned.

        Examples
        --------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> Tind = pn.throats()
        >>> Tind[0:5]
        array([0, 1, 2, 3, 4])

        """
        if labels == 'all':
            Nt = sp.shape(self['throat.all'])[0]
            ind = sp.arange(0,Nt)
        else:
            ind = self._get_indices(element='throat',labels=labels,mode=mode)
        return ind

    @property
    def Ts(self):
        r"""
        A shortcut to get a list of all throats on the object
        """
        return self.throats()

    def _tomask(self,locations,element):
        r"""
        This is a generalized version of tomask that accepts a string of
        'pore' or 'throat' for programmatic access.
        """
        if sp.shape(locations)[0] == 0:
            return sp.zeros_like(self._get_indices(element=element),dtype=bool)
        if element in ['pore','pores']:
            Np = sp.shape(self['pore.all'])[0]
            pores = sp.array(locations,ndmin=1)
            mask = sp.zeros((Np,),dtype=bool)
            mask[pores] = True
        if element in ['throat','throats']:
            Nt = sp.shape(self['throat.all'])[0]
            throats = sp.array(locations,ndmin=1)
            mask = sp.zeros((Nt,),dtype=bool)
            mask[throats] = True
        return mask

    def tomask(self,pores=None,throats=None):
        r"""
        Convert a list of pore or throat indices into a boolean mask of the
        correct length

        Parameters
        ----------
        pores or throats : array_like
            List of pore or throat indices

        Returns
        -------
        mask : array_like
            A boolean mask of length Np or Nt with True in the locations of
            pores or throats received.

        """
        if pores is not None:
            mask = self._tomask(element='pore',locations=pores)
        if throats is not None:
            mask = self._tomask(element='throat',locations=throats)
        return mask

    def toindices(self,mask):
        r"""
        Convert a boolean mask a list of pore or throat indices

        Parameters
        ----------
        mask : array_like booleans
            A boolean array with True at locations where indices are desired.
            The appropriate indices are returned based an the length of mask,
            which must be either Np or Nt long.

        Returns
        -------
        indices : array_like
            A list of pore or throat indices corresponding the locations where
            the received mask was True.

        Notes
        -----
        This behavior could just as easily be accomplished by using the mask
        in pn.pores()[mask] or pn.throats()[mask].  This method is just a thin
        convenience function and is a compliment to tomask().

        """
        mask = sp.array(mask,ndmin=1)
        if sp.shape(mask)[0] == self.num_pores():
            indices = self.pores()[mask]
        elif sp.shape(mask)[0] == self.num_throats():
            indices = self.throats()[mask]
        else:
            raise Exception('Mask received was neither Np nor Nt long')
        return indices

    def interpolate_data(self,data):
        r"""
        Determines a pore (or throat) property as the average of it's neighboring
        throats (or pores)

        Parameters
        ----------
        data : array_like
            A list of specific values to be interpolated.  List MUST be either
            Np or Nt long

        Returns
        -------
        An array containing interpolated pore (or throat) data

        Notes
        -----
        - This uses an unweighted average, without attempting to account for distances or sizes of pores and throats.
        - Only one of pores, throats OR data are accepted

        """
        mro = [module.__name__ for module in self.__class__.__mro__]
        if 'GenericNetwork' in mro:
            net = self
            Ts = net.throats()
            Ps = net.pores()
            label = 'all'
        elif ('GenericPhase' in mro) or ('GenericAlgorithm' in mro):
            net = self._net
            Ts = net.throats()
            Ps = net.pores()
            label = 'all'
        elif ('GenericGeometry' in mro) or ('GenericPhysics' in mro):
            net = self._net
            Ts = net.throats(self.name)
            Ps = net.pores(self.name)
            label = self.name
        if sp.shape(data)[0] == self.Nt:
            #Upcast data to full network size
            temp = sp.ones((net.Nt,))*sp.nan
            temp[Ts] = data
            data = temp
            temp = sp.ones((net.Np,))*sp.nan
            for pore in Ps:
                neighborTs = net.find_neighbor_throats(pore)
                neighborTs = net.filter_by_label(throats=neighborTs,labels=label)
                temp[pore] = sp.mean(data[neighborTs])
            values = temp[Ps]
        elif sp.shape(data)[0] == self.Np:
            #Upcast data to full network size
            temp = sp.ones((net.Np,))*sp.nan
            temp[Ps] = data
            data = temp
            Ps12 = net.find_connected_pores(throats=Ts,flatten=False)
            values = sp.mean(data[Ps12],axis=1)
        else:
            logger.error('Received data was an ambiguous length')
            raise Exception()
        return values

    def _interleave_data(self,prop,sources):
        r"""
        Retrieves requested property from associated objects, to produce a full
        Np or Nt length array.

        Parameters
        ----------
        prop : string
            The property name to be retrieved
        sources : list
            List of object names OR objects from which data is retrieved

        Returns
        -------
        A full length (Np or Nt) array of requested property values.

        Notes
        -----
        This makes an effort to maintain the data 'type' when possible; however
        when data is missing this can be tricky.  Float and boolean data is
        fine, but missing ints are converted to float when nans are inserted.

        Examples
        --------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> Ps = pn.pores('top',mode='not')
        >>> Ts = pn.find_neighbor_throats(pores=Ps,mode='intersection',flatten=True)
        >>> geom = OpenPNM.Geometry.TestGeometry(network=pn,pores=Ps,throats=Ts)
        >>> Ps = pn.pores('top')
        >>> Ts = pn.find_neighbor_throats(pores=Ps,mode='not_intersection')
        >>> boun = OpenPNM.Geometry.Boundary(network=pn,pores=Ps,throats=Ts)
        >>> geom['pore.test_int'] = sp.random.randint(0, 100, geom.Np)
        >>> print(pn['pore.test_int'].dtype)
        float64
        >>> boun['pore.test_int'] = sp.ones(boun.Np).astype(int)
        >>> boun['pore.test_int'] = sp.rand(boun.Np)<0.5
        >>> print(pn['pore.test_int'].dtype)
        bool
        >>> geom['pore.test_bool'] = sp.rand(geom.Np)<0.5
        >>> print(pn['pore.test_bool'].dtype)
        bool
        >>> boun['pore.test_bool'] = sp.ones(boun.Np).astype(int)
        >>> print(pn['pore.test_bool'].dtype)
        bool
        >>> boun['pore.test_bool'] = sp.rand(boun.Np)<0.5
        >>> print(pn['pore.test_bool'].dtype)
        bool
        """
        element = prop.split('.')[0]
        temp = sp.ndarray((self._count(element)))
        nan_locs = sp.ndarray((self._count(element)), dtype='bool')
        nan_locs.fill(False)
        bool_locs = sp.ndarray((self._count(element)), dtype='bool')
        bool_locs.fill(False)
        dtypes = []
        dtypenames = []
        prop_found = False  #Flag to indicate if prop was found on a sub-object
        values_dim=0
        for item in sources:
            #Check if sources were given as list of objects OR names
            try: item.name
            except: item = self._find_object(obj_name=item)
            locations = self._get_indices(element=element,labels=item.name,mode='union')
            if prop not in item.keys():
                values = sp.ones_like(temp[locations])*sp.nan
                dtypenames.append('nan')
                dtypes.append(sp.dtype(bool))
                nan_locs[locations]=True
            else:
                prop_found = True
                values = item[prop]
                dtypenames.append(values.dtype.name)
                dtypes.append(values.dtype)
                if values.dtype == 'bool':
                    bool_locs[locations]=True
                try: values_dim = sp.shape(values)[1]
                except: pass
            if values_dim > 0:
                try:
                    temp_dim = sp.shape(temp)[1]
                    if temp_dim != values_dim:
                        logger.warning(prop+' data has different dimensions, consider revising data in object '+str(item.name))
                except:
                    temp = sp.ndarray([self._count(element),values_dim])
            if values.dtype == 'object' and temp.dtype != 'object':
                temp = temp.astype('object')
            temp[locations] = values  #Assign values
        #Check if requested prop was found on any sub-objects
        if prop_found == False:
            raise KeyError(prop)
        #Analyze and assign data type
        if sp.all([t in ['bool','nan'] for t in dtypenames]):  # If all entries are 'bool' (or 'nan')
            temp = sp.array(temp,dtype='bool')
            if sp.sum(nan_locs)>0:
                temp[nan_locs]=False
        elif sp.all([t == dtypenames[0] for t in dtypenames]) :  # If all entries are same type
            temp = sp.array(temp,dtype=dtypes[0])
        elif sp.all([t in ['int','nan','float','int32','int64','float32','float64','bool'] for t in dtypenames]):  # If all entries are 'bool' (or 'nan')
            if 'bool' in dtypenames:
                temp = sp.array(temp,dtype='bool')
                temp[~bool_locs]=False
                logger.info(prop+' has been converted to bool, some data may be lost')
            else:
                temp = sp.array(temp,dtype='float')
                logger.info(prop+' has been converted to float.')
        elif sp.all([t in ['object','nan'] for t in dtypenames]):  # If all entries are 'bool' (or 'nan')
            pass
        else:
            temp = sp.array(temp,dtype=max(dtypes))
            logger.info('Data type of '+prop+' differs between sub-objects...converting to larger data type')
        return temp

    def num_pores(self,labels='all',mode='union'):
        r"""
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
        num_throats
        count

        Examples
        --------
        >>> import OpenPNM
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

        """
        if labels == 'all':
            Np = sp.shape(self.get('pore.all'))[0]
        else:
            #convert string to list, if necessary
            if type(labels) == str:
                labels = [labels]
            #Count number of pores of specified type
            Ps = self.pores(labels=labels,mode=mode)
            Np = sp.shape(Ps)[0]
        return Np

    @property
    def Np(self):
        r"""
        A shortcut to query the total number of pores on the object'
        """
        return self.num_pores()

    def num_throats(self,labels='all',mode='union'):
        r"""
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
        num_pores
        count

        Examples
        --------
        >>> import OpenPNM
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

        """
        if labels == 'all':
            Nt = sp.shape(self.get('throat.all'))[0]
        else:
            #convert string to list, if necessary
            if type(labels) == str: labels = [labels]
            #Count number of pores of specified type
            Ts = self.throats(labels=labels,mode=mode)
            Nt = sp.shape(Ts)[0]
        return Nt

    @property
    def Nt(self):
        r"""
        A shortcut to query the total number of throats on the object'
        """
        return self.num_throats()

    def _count(self,element=None):
        r"""
        Returns a dictionary containing the number of pores and throats in
        the network, stored under the keys 'pore' or 'throat'

        Parameters
        ----------
        element : string, optional
            Can be either 'pore' , 'pores', 'throat' or 'throats', which
            specifies which count to return.

        Returns
        -------
        A dictionary containing the number of pores and throats under the
        'pore' and 'throat' key respectively.

        See Also
        --------
        num_pores
        num_throats

        Notes
        -----
        The ability to send plurals is useful for some types of 'programmatic'
        access.  For instance, the standard argument for locations is pores
        or throats.  If these are bundled up in a **kwargs dict then you can
        just use the dict key in count() without removing the 's'.

        Examples
        --------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn._count('pore')
        125
        >>> pn._count('throat')
        300
        """
        if element in ['pore','pores']:
            temp = self.num_pores()
        elif element in ['throat','throats']:
            temp = self.num_throats()
        elif element is None:
            temp = {}
            temp['pore'] = self.num_pores()
            temp['throat'] = self.num_throats()
        return temp

    def _set_locations(self,element,locations,mode='add'):
        r"""
        Private method used for assigning Geometry and Physics objects to
        specified locations

        Parameters
        ----------
        element : string
            Either 'pore' or 'throat' indicating which type of element is being
            work upon
        locations : array_like
            The pore or throat locations in terms of Network numbering to add
            (or remove) from the object
        mode : string
            Either 'add' or 'remove', the default is add.

        Examples
        --------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.Np
        125
        >>> geom = OpenPNM.Geometry.GenericGeometry(network=pn,pores=sp.arange(5,125),throats=pn.Ts)
        >>> [geom.Np, geom.Nt]
        [120, 300]
        >>> geom['pore.dummy'] = True
        >>> health = pn.check_geometry_health()
        >>> pores = health['undefined_pores']
        >>> geom.set_locations(pores=pores)
        >>> [geom.Np, geom.Nt]
        [125, 300]
        >>> geom.pores(labels='dummy',mode='not')  # Dummy as assigned BEFORE these pores were added
        array([0, 1, 2, 3, 4])
        >>> geom.set_locations(pores=pores,mode='remove')
        >>> [geom.Np, geom.Nt]
        [120, 300]
        >>> geom.num_pores(labels='dummy',mode='not')  # All pores without 'dummy' label are gone
        0
        """
        net = self._net
        if self._isa('Geometry'):
            boss_obj = self._net
            co_objs = boss_obj.geometries()
        elif self._isa('Physics'):
            boss_obj = self._phases[0]
            co_objs = boss_obj.physics()
        else:
            logger.warning('Setting locations only applies to Geometry or Physics objects')
            return

        if mode == 'add':
            # Check if any constant values exist on the object
            for item in self.props():
                if (item not in self.models.keys()) or \
                   (self.models[item]['regen_mode'] == 'constant'):
                    logger.critical('Constant models found on object,' +
                                    'models must be rerun manually')
            # Ensure locations are not already assigned to another object
            temp = sp.zeros((net._count(element), ), dtype=bool)
            for key in co_objs:
                temp += net[element+'.'+key]
            overlaps = sp.sum(temp*net._tomask(locations=locations,
                                               element=element))
            if overlaps > 0:
                self.controller.purge_object(self)
                raise Exception('Some of the given '+element+'s overlap with an existing object')

            # Store original Network indices for later use
            old_inds = sp.copy(net[element+'.'+self.name])

            # Create new 'all' label for new size
            new_len = self._count(element=element) + sp.size(locations)
            # Initialize new 'all' array
            self.update({element+'.all': sp.ones((new_len, ), dtype=bool)})

            # Set locations in Network (and Phase) dictionary
            if element+'.'+self.name not in net.keys():
                net[element+'.'+self.name] = False
            net[element+'.'+self.name][locations] = True
            if element+'.'+self.name not in boss_obj.keys():
                boss_obj[element+'.'+self.name] = False
            boss_obj[element+'.'+self.name][locations] = True

            # Increase size of labels (add False at new locations)
            blank = ~sp.copy(self[element+'.all'])
            labels = self.labels()
            labels.remove(element+'.all')
            for item in labels:
                if item.split('.')[0] == element:
                    blank[old_inds] = self[item]
                    self.update({item: blank[net[element+'.all']]})

        if mode == 'remove':
            self_inds = boss_obj._map(element=element,
                                      locations=locations,
                                      target=self)
            keep = ~self._tomask(locations=self_inds, element=element)
            for item in self.keys():
                if item.split('.')[0] == element:
                    temp = self[item][keep]
                    self.update({item: temp})
            # Set locations in Network dictionary
            net[element+'.'+self.name][locations] = False
            boss_obj[element+'.'+self.name][locations] = False

        # Finally, regenerate models to correct the length of all prop array
        self.models.regenerate()

    def _map(self, element, locations, target, return_mapping=False):
        r"""
        """
        # Initialize things
        locations = sp.array(locations, ndmin=1)
        mapping = {}

        # Analyze input object's relationship
        if self._net == target._net:  # Objects are siblings...easy
            maskS = self._net[element+'.'+self.name]
            maskT = target._net[element+'.'+target.name]
        else:  # One or more of the objects is a clone
            if self._parent is None:  # Self is parent object
                maskS = self._net[element+'.'+self.name]
                maskT = ~self._net[element+'.all']
                tempT = target._net[element+'.'+target.name]
                inds = target._net[element+'.'+self._net.name][tempT]
                maskT[inds] = True
            if target._parent is None:  # Target is parent object
                maskT = target._net[element+'.'+target.name]
                maskS = ~target._net[element+'.all']
                tempS = self._net[element+'.'+self.name]
                inds = self._net[element+'.'+target._net.name][tempS]
                maskS[inds] = True

        # Convert source locations to Network indices
        temp = sp.zeros(sp.shape(maskS), dtype=int)-1
        temp[maskS] = self._get_indices(element=element)
        locsS = sp.where(sp.in1d(temp, locations))[0]
        mapping['source'] = locations

        # Find locations in target
        temp = sp.zeros(sp.shape(maskT), dtype=int)-1
        temp[maskT] = target._get_indices(element=element)
        locsT = temp[locsS]
        mapping['target'] = locsT

        # Find overlapping locations in source and target to define mapping
        keep = (locsS >= 0)*(locsT >= 0)
        mapping['source'] = mapping['source'][keep]
        mapping['target'] = mapping['target'][keep]

        # Return results as an arrary or one-to-one mapping if requested
        if return_mapping is True:
            return mapping
        else:
            if sp.sum(locsS >= 0) < sp.shape(sp.unique(locations))[0]:
                raise Exception('Some locations not found on Source object')
            if sp.sum(locsT >= 0) < sp.shape(sp.unique(locations))[0]:
                raise Exception('Some locations not found on Target object')
            return mapping['target']

    def map_pores(self, target=None, pores=None, return_mapping=False):
        r"""
        Accepts a list of pores from the caller object and maps them onto the
        given target object

        Parameters
        ----------
        pores : array_like
            The list of pores on the caller object.  If no pores are supplied
            then all the pores of the calling object are used.

        target : OpenPNM object, optional
            The object for which a list of pores is desired.  If no object is
            supplied then the object's associated Network is used.

        return_mapping : boolean (default is False)
            If True, a dictionary containing 'source' locations, and 'target'
            locations is returned.  Any 'source' locations not found in the
            'target' object are removed from the list.

        Returns
        -------
        pores : array_like
            A list of pores mapped onto the target object

        Examples
        --------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> Ps = pn.pores(labels=['top','left'],mode='intersection')
        >>> Ps
        array([100, 101, 102, 103, 104])
        >>> geom = OpenPNM.Geometry.GenericGeometry(network=pn,pores=Ps)
        >>> geom.Ps
        array([0, 1, 2, 3, 4])
        >>> geom.map_pores(target=pn,pores=geom.Ps)
        array([100, 101, 102, 103, 104])
        >>> pn.map_pores(target=geom,pores=Ps)
        array([0, 1, 2, 3, 4])
        """
        if pores is None:
            pores = self.Ps
        if target is None:
            if self._net is None:
                target = self
            else:
                target = self._net
        Ps = self._map(element='pore',
                       locations=pores,
                       target=target,
                       return_mapping=return_mapping)
        return Ps

    def map_throats(self,
                    target=None,
                    throats=None,
                    return_mapping=False):
        r"""
        Accepts a list of throats from the caller object and maps them onto the
        given target object

        Parameters
        ----------
        throats : array_like
            The list of throats on the caller object.  If no throats are
            supplied then all the throats of the calling object are used.

        target : OpenPNM object, optional
            The object for which a list of pores is desired.  If no object is
            supplied then the object's associated Network is used.

        return_mapping : boolean (default is False)
            If True, a dictionary containing 'source' locations, and 'target'
            locations is returned.  Any 'source' locations not found in the
            'target' object are removed from the list.

        Returns
        -------
        throats : array_like
            A list of throats mapped onto the target object

        Examples
        --------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> Ts = pn.throats(labels=['top','left'],mode='intersection')
        >>> Ts
        array([260, 262, 264, 266])
        >>> geom = OpenPNM.Geometry.GenericGeometry(network=pn,throats=Ts)
        >>> geom.Ts
        array([0, 1, 2, 3])
        >>> geom.map_throats(target=pn,throats=geom.Ts)
        array([260, 262, 264, 266])
        >>> pn.map_throats(target=geom,throats=Ts)
        array([0, 1, 2, 3])
        """
        if throats is None:
            throats = self.Ts
        if target is None:
            if self._net is None:
                target = self
            else:
                target = self._net
        Ts = self._map(element='throat',
                       locations=throats,
                       target=target,
                       return_mapping=return_mapping)
        return Ts

    Tnet = property(fget=map_throats)
    Pnet = property(fget=map_pores)

    def _isa(self, keyword=None, obj=None):
        r"""
        """
        if keyword is None:
            mro = [item.__name__ for item in self.__class__.__mro__]
        if obj is None:
            query = False
            mro = [item.__name__ for item in self.__class__.__mro__]
            if keyword in ['net', 'Network', 'GenericNetwork']:
                if 'GenericNetwork' in mro:
                    query = True
            elif keyword in ['geom', 'Geometry', 'GenericGeometry']:
                if 'GenericGeometry' in mro:
                    query = True
            elif keyword in ['phase', 'Phase', 'GenericPhase']:
                if 'GenericPhase' in mro:
                    query = True
            elif keyword in ['phys', 'Physics', 'GenericPhysics']:
                if 'GenericPhysics' in mro:
                    query = True
            elif keyword in ['alg', 'Algorithm', 'GenericAlgorithm']:
                if 'GenericAlgorithm' in mro:
                    query = True
            elif keyword in ['clone']:
                if self._net is None:
                    if self._parent is not None:
                        query = True
                else:
                    if self._net._parent is not None:
                        query = True
            return query
        else:
            query = False
            if keyword in ['sibling']:
                if (self._isa('net')) and (obj._net is self):
                    query = True
                elif (obj._isa('net')) and (self._net is obj):
                    query = True
                elif self._net is obj._net:
                    query = True
            return query

    def check_data_health(self, props=[], element=''):
        r"""
        Check the health of pore and throat data arrays.

        Parameters
        ----------
        element : string, optional
            Can be either 'pore' or 'throat', which will limit the checks to
            only those data arrays.
        props : list of pore (or throat) properties, optional
            If given, will limit the health checks to only the specfied
            properties.  Also useful for checking existance.

        Returns
        -------
        Returns a HealthDict object which a basic dictionary with an added
        ``health`` attribute that is True is all entries in the dict are
        deemed healthy (empty lists), or False otherwise.

        Examples
        --------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> health_check = pn.check_data_health()
        >>> health_check.health
        True
        """
        health = Tools.HealthDict()
        if props == []:
            props = self.props(element)
        else:
            if type(props) == str:
                props = [props]
        for item in props:
            health[item] = []
            try:
                if sp.sum(sp.isnan(self[item])) > 0:
                    health[item] = 'Has NaNs'
                elif sp.shape(self[item])[0] != self._count(item.split('.')[0]):
                    health[item] = 'Wrong Length'
            except:
                health[item] = 'Does not exist'
        return health

    def __str__(self):
        horizonal_rule = '-' * 60
        lines = [horizonal_rule]
        lines.append(self.__module__.replace('__', '') + ': \t' + self.name)
        lines.append(horizonal_rule)
        lines.append("{0:<5s} {1:<35s} {2:<10s}".format('#',
                                                        'Properties',
                                                        'Valid Values'))
        lines.append(horizonal_rule)
        props = self.props()
        props.sort()
        for i, item in enumerate(props):
            if self[item].dtype != object:
                prop = item
                if len(prop) > 35:
                    prop = prop[0:32] + '...'
                required = self._count(item.split('.')[0])
                a = sp.isnan(self[item])
                defined = sp.shape(self[item])[0] - a.sum(axis=0,
                                                          keepdims=(a.ndim-1)==0)[0]
                lines.append("{0:<5d} {1:<35s} {2:>5d} / {3:<5d}".format(i + 1,
                                                                         prop,
                                                                         defined,
                                                                         required))
        lines.append(horizonal_rule)
        lines.append("{0:<5s} {1:<35s} {2:<10s}".format('#',
                                                        'Labels',
                                                        'Assigned Locations'))
        lines.append(horizonal_rule)
        labels = self.labels()
        labels.sort()
        for i, item in enumerate(labels):
            prop = item
            if len(prop) > 35:
                prop = prop[0:32] + '...'
            lines.append("{0:<5d} {1:<35s} {2:<10d}".format(i + 1,
                                                            prop,
                                                            sp.sum(self[item])))
        lines.append(horizonal_rule)
        return '\n'.join(lines)
