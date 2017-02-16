"""
###############################################################################
Core:  Core Data Class
###############################################################################
"""
from OpenPNM.Base import Workspace
import string
import random
import scipy as sp
import scipy.constants
from OpenPNM.Base import logging, Tools
from OpenPNM.Base import ModelsDict
logger = logging.getLogger()
mgr = Workspace()


class Core(dict):
    r"""
    Contains methods for working with the data in the OpenPNM dictionaries
    """

    def __new__(typ, *args, **kwargs):
        obj = dict.__new__(typ, *args, **kwargs)
        obj.update({'pore.all': sp.array([], ndmin=1, dtype=bool)})
        obj.update({'throat.all': sp.array([], ndmin=1, dtype=bool)})
        # Initialize phase, physics, and geometry tracking lists
        obj._name = None
        obj.phases = Tools.ObjectContainer()
        obj.geometries = Tools.ObjectContainer()
        obj.physics = Tools.ObjectContainer()
        obj.network = Tools.ObjectContainer()
        obj._parent = None
        # Initialize ordered dict for storing property models
        obj.models = ModelsDict()
        return obj

    def __init__(self, name=None, **kwargs):
        super().__init__()
        logger.debug('Initializing Core class')
        self.name = name

    def __repr__(self):
        return '<%s.%s object at %s>' % (
            self.__class__.__module__,
            self.__class__.__name__,
            hex(id(self)))

    def __eq__(self, other):
        if hex(id(self)) == hex(id(other)):
            return True
        else:
            return False

    def __setitem__(self, key, value):
        r"""
        This is a subclass of the default __setitem__ behavior.  The main aim
        is to limit what type and shape of data can be written to protect
        the integrity of the network.  Specifically, this means only Np or Nt
        long arrays can be written, and they must be called 'pore.___' or
        'throat.___'.  Also, any scalars are cast into full length vectors.


        Example
        -------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn['pore.example_property'] = 100
        >>> pn['pore.example_property'][0]
        100

        """
        # Enforce correct dict naming
        element = self._parse_element(key.split('.')[0], single=True)
        # Convert value to an ndarray
        value = sp.array(value, ndmin=1)
        # Skip checks for 'coords', 'conns'
        if key in ['pore.coords', 'throat.conns']:
            super(Core, self).__setitem__(key, value)
            return
        # Skip checks for protected props, and prevent changes if defined
        protected_keys = ['all']
        if key.split('.')[1] in protected_keys:
            if key in self.keys():
                if sp.shape(self[key]) == (0,):
                    logger.debug(key+' is being defined.')
                    super(Core, self).__setitem__(key, value)
                else:
                    logger.warning(key+' is already defined.')
            else:
                logger.debug(key+' is being defined.')
                super(Core, self).__setitem__(key, value)
            return
        # Write value to dictionary
        if sp.shape(value)[0] == 1:  # If value is scalar
            logger.debug('Broadcasting scalar value into vector: '+key)
            value = sp.ones((self._count(element), ), dtype=value.dtype)*value
            super(Core, self).__setitem__(key, value)
        elif sp.shape(value)[0] == self._count(element):
            logger.debug('Updating vector: '+key)
            super(Core, self).__setitem__(key, value)
        else:
            if self._count(element) == 0:
                self.update({key: value})
            else:
                logger.warning('Cannot write vector with an array of the ' +
                               'wrong length: '+key)
                # TODO: This should probably raise the following exception
                # raise Exception('Cannot write vector of the wrong length')

    def _get_mgr(self):
        if self in mgr.values():
            return mgr
        else:
            return {}

    def _set_mgr(self, mgr):
        if self not in mgr.values():
            mgr.update({self.name: self})

    workspace = property(fget=_get_mgr, fset=_set_mgr)

    def _set_name(self, name):
        if name in mgr.keys():
            raise Exception('An object named '+name+' already exists')
        elif name is None:
            name = ''.join(random.choice(string.ascii_uppercase +
                                         string.ascii_lowercase +
                                         string.digits) for _ in range(5))
            name = self.__class__.__name__ + '_' + name
        elif self._name is not None:
            logger.info('Changing the name of '+self.name+' to '+name)
            # Check if name collides with any arrays in the simulation
            if mgr._validate_name(name):
                # Rename any label arrays
                for item in self._simulation():
                    if 'pore.'+self.name in item.keys():
                        item['pore.'+name] = item.pop('pore.'+self.name)
                    if 'throat.'+self.name in item.keys():
                        item['throat.'+name] = item.pop('throat.'+self.name)
            else:
                raise Exception('The provided name is already in use')
        # Remove reference to object under old name, if present
        for item in list(mgr.items()):
            if item[1] is self:
                mgr.pop(item[0])
        # Add object to workspace under new name
        mgr.update({name: self})
        self._name = name

    def _get_name(self):
        return self._name

    name = property(_get_name, _set_name)

    def _simulation(self):
        temp = []
        temp += [self._net]
        temp += self._net._phases
        temp += self._net._geometries
        temp += self._net._physics
        return temp

    def clear(self, mode='complete'):
        r"""
        A subclassed version of the standard dict's clear method.  This can be
        used to selectively clear certain aspects of the object, including
        properties, labels and/or models.  It can also clear everything,
        except for the 'pore.all' and 'throat.all' labels which are required
        for object to remain functional.

        Parameters
        ----------
        mode : string of list of strings
            This controls what is cleared from the object.  Options are:

            **'props'** : Removes all numerical property values from the object
            dictionary.

            **'labels'** : Removes all labels from the object dictionary

            **'models'** : Removes all pore scale models from the object's
            models dictionary (object.models)

            **'complete'** : Removes all of the above AND sets the \'pore.all\'
            and \'throat.all\' labels to zero length.  This also removes any
            pore and throat locations that were previously set.  This mode
            should be used carefully since it can break some subtle aspects
            of the framework; it is meant for advanced users and developers.

        Notes
        -----
        The first three modes listed can be combined by sending a list
        containing all desired modes.  The \'complete\' mode essentially calls
        all three so need not be combined with any other modes.
        """
        allowed = ['props', 'labels', 'models', 'complete']
        mode = self._parse_mode(mode=mode, allowed=allowed)
        if 'complete' in mode:
            if self._isa('Geometry') or self._isa('Physics'):
                self.set_locations(pores=self.Pnet,
                                   throats=self.Tnet,
                                   mode='remove')
            super().clear()
            self.models.clear()
            self.update({'throat.all': sp.array([], ndmin=1, dtype=bool)})
            self.update({'pore.all': sp.array([], ndmin=1, dtype=bool)})
        if 'props' in mode:
            for item in self.props():
                del self[item]
        if 'labels' in mode:
            for item in self.labels():
                if item not in ['pore.all', 'throat.all']:
                    del self[item]
        if 'models' in mode:
            self.models.clear()

    def _find_object(self, obj_name='', obj_type=''):
        all_dicts = {}
        all_dicts.update(self.geometries)
        all_dicts.update(self.physics)
        all_dicts.update(self.phases)
        all_dicts.update({self._net.name: self._net})
        if obj_name != '':
            return all_dicts.get(obj_name)
        if obj_type != '':
            objs = []
            for item in all_dicts.values():
                if item._isa(obj_type):
                    objs.append(item.name)
            return objs

    @property
    def _geometries(self):
        return list(self.geometries.values())

    @property
    def _phases(self):
        return list(self.phases.values())

    @property
    def _physics(self):
        return list(self.physics.values())

    @property
    def _net(self):
        return list(self.network.values())[0]

    # -------------------------------------------------------------------------
    """Model Manipulation Methods"""
    # -------------------------------------------------------------------------
    # Note: These methods have been moved to the ModelsDict class but are left
    # here for backward compatibility
    def add_model(self, propname, model, regen_mode='normal', **kwargs):
        self.models.add(propname=propname,
                        model=model,
                        regen_mode=regen_mode,
                        **kwargs)

    add_model.__doc__ = ModelsDict.add.__doc__

    def regenerate(self, props='', mode='inclusive'):
        self.models.regenerate(props=props, mode=mode)

    regenerate.__doc__ = ModelsDict.regenerate.__doc__

    # -------------------------------------------------------------------------
    """Data Query Methods"""
    # -------------------------------------------------------------------------
    def props(self, element=None, mode='all'):
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

            **'all'** : Returns all properties on the object

            **'models'** : Returns only properties that are associated with a
            model

            **'constants'** : Returns only properties that are set as constant
            values

        deep : Boolean (default if False)
            If True, all properties on the object and all sub-objects are
            returned. For instance, all Geometry properties will be returned
            along with all Network properties, and all Physics properties will
            be returned with all Phase properties.  This arg has no effect
            when this query is called from a Geometry or Phase object.

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
        >>> #pn.props()
        ['pore.coords', 'pore.index', 'throat.conns']
        """
        # Parse Inputs
        allowed = ['all', 'models', 'constants', 'deep']
        mode = self._parse_mode(mode=mode, allowed=allowed)
        element = self._parse_element(element=element)
        # Prepare lists of each type of array
        props = [item for item in self.keys() if self[item].dtype != bool]
        models = list(self.models.keys())
        constants = [item for item in props if item not in models]
        # Execute desired array lookup
        vals = Tools.PrintableList()
        if 'deep' in mode:
            logger.warning('The \'deep\' mode only works when called from a ' +
                           'Network or Phase object')
        if 'all' in mode:
            temp = [item for item in props if item.split('.')[0] in element]
            vals.extend(temp)
            return vals
        if 'models' in mode:
            temp = [item for item in models if item.split('.')[0] in element]
            vals.extend(temp)
        if 'constants' in mode:
            temp = [item for item in constants
                    if item.split('.')[0] in element]
            vals.extend(temp)
        return vals

    def _get_labels(self, element, locations, mode):
        r"""
        This is the actual label getter method, but it should not be called
        directly.  Wrapper methods have been created, use ``labels``.
        """
        # Parse inputs
        locations = self._parse_locations(locations)
        allowed = ['none', 'union', 'intersection', 'not', 'count', 'mask',
                   'difference']
        mode = self._parse_mode(mode=mode, allowed=allowed, single=True)
        element = self._parse_element(element=element)
        # Collect list of all pore OR throat labels
        a = set([k for k in self.keys() if k.split('.')[0] == element[0]])
        b = set([k for k in self.keys() if self[k].dtype == bool])
        labels = list(a.intersection(b))
        labels.sort()
        labels = sp.array(labels)  # Convert to ND-array for following checks
        arr = sp.zeros((sp.shape(locations)[0], len(labels)), dtype=bool)
        col = 0
        for item in labels:
            arr[:, col] = self[item][locations]
            col = col + 1
        if mode in ['count']:
            return sp.sum(arr, axis=1)
        if mode in ['union']:
            temp = labels[sp.sum(arr, axis=0) > 0]
            temp.tolist()
            return Tools.PrintableList(temp)
        if mode in ['intersection']:
            temp = labels[sp.sum(arr, axis=0) == sp.shape(locations, )[0]]
            temp.tolist()
            return Tools.PrintableList(temp)
        if mode in ['not', 'difference']:
            temp = labels[sp.sum(arr, axis=0) != sp.shape(locations, )[0]]
            temp.tolist()
            return Tools.PrintableList(temp)
        if mode in ['mask']:
            return arr
        if mode in ['none']:
            temp = sp.ndarray((sp.shape(locations, )[0], ), dtype=object)
            for i in sp.arange(0, sp.shape(locations, )[0]):
                temp[i] = list(labels[arr[i, :]])
            return temp
        else:
            logger.error('unrecognized mode:'+str(mode))

    def labels(self, pores=[], throats=[], element=None, mode='union'):
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

            **'none'** : An N x Li list of all labels applied to each input
            pore (or throats). Li can vary betwen pores (and throats)

            **'union'** : A list of labels applied to ANY of the given pores
            (or throats)

            **'intersection'** : Label applied to ALL of the given pores
            (or throats)

            **'not'** : Labels NOT applied to ALL pores (or throats)

            **'count'** : The number of labels on each pores (or throats)

            **'mask'**: returns an N x Lt array, where each row corresponds to
            a pore (or throat) location, and each column contains the truth
            value for the existance of labels as returned from
            ``labels(element='pores')``.

        Returns
        -------
        A list containing the dictionary keys on the object limited by the
        specified ``mode``.

        Examples
        --------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.labels(pores=[0, 1, 5, 6])
        ['pore.all', 'pore.bottom', 'pore.front', 'pore.left']
        >>> pn.labels(pores=[0, 1, 5, 6], mode='intersection')
        ['pore.all', 'pore.bottom']
        """
        labels = Tools.PrintableList()
        # Short-circuit query when no pores or throats are given
        if (sp.size(pores) == 0) and (sp.size(throats) == 0):
            element = self._parse_element(element=element)
            for item in element:
                a = set([key for key in self.keys()
                         if key.split('.')[0] == item])
                b = set([key for key in self.keys()
                         if self[key].dtype == bool])
                labels.extend(list(a.intersection(b)))
        elif (sp.size(pores) > 0) and (sp.size(throats) > 0):
            raise Exception('Cannot perform label query on pores and ' +
                            'throats simultaneously')
        elif sp.size(pores) > 0:
            labels = self._get_labels(element='pore',
                                      locations=pores,
                                      mode=mode)
        elif sp.size(throats) > 0:
            labels = self._get_labels(element='throat',
                                      locations=throats,
                                      mode=mode)
        return labels

    def filter_by_label(self, pores=[], throats=[], labels=None, mode='union'):
        r"""
        Returns which of the supplied pores (or throats) has the specified
        label

        Parameters
        ----------
        pores, or throats : array_like
            List of pores or throats to be filtered

        labels : list of strings
            The labels to apply as a filter

        mode : string
            Controls how the filter is applied.  Options include:

            **'union'** : (default) All locations with ANY of the given labels
            are kept.

            **'intersection'** : Only locations with ALL the given labels are
            kept.

            **'not_intersection'** : Only locations with exactly one of the
            given labels are kept.

            **'not'** : Only locations with none of the given labels are kept.

        See Also
        --------
        pores
        throats

        Returns
        -------
        A list of pores (or throats) that have been filtered according the
        given criteria.  The returned list is a subset of the received list of
        pores (or throats).

        Examples
        --------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.filter_by_label(pores=[0,1,5,6], labels='left')
        array([0, 1])
        >>> Ps = pn.pores(['top', 'bottom', 'front'], mode='union')
        >>> pn.filter_by_label(pores=Ps, labels=['top', 'front'],
        ...                    mode='intersection')
        array([100, 105, 110, 115, 120])
        """
        allowed = ['union', 'intersection', 'not_intersection', 'not']
        mode = self._parse_mode(mode=mode, allowed=allowed, single=True)
        # Convert inputs to locations and element
        if (sp.size(throats) > 0) and (sp.size(pores) > 0):
            raise Exception('Can only filter either pores OR labels per call')
        if sp.size(pores) > 0:
            element = 'pore'
            locations = self._parse_locations(pores)
        elif sp.size(throats) > 0:
            element = 'throat'
            locations = self._parse_locations(throats)
        else:
            return(sp.array([], dtype=int))
        # Do it
        labels = self._parse_labels(labels=labels, element=element)
        labels = [element+'.'+item.split('.')[-1] for item in labels]
        all_locs = self._get_indices(element=element, labels=labels, mode=mode)
        mask = self._tomask(locations=all_locs, element=element)
        ind = mask[locations]
        return locations[ind]

    def _get_indices(self, element, labels='all', mode='union'):
        r"""
        This is the actual method for getting indices, but should not be called
        directly.  Use pores or throats instead.
        """
        # Parse and validate all input values
        allowed = ['union', 'intersection', 'not_intersection', 'not',
                   'difference']
        mode = self._parse_mode(mode=mode, allowed=allowed, single=True)
        element = self._parse_element(element, single=True)
        labels = self._parse_labels(labels=labels, element=element)
        if element+'.all' not in self.keys():
            raise Exception('Cannot proceed without {}.all'.format(element))

        # Begin computing label array
        if mode in ['union']:
            union = sp.zeros_like(self[element+'.all'], dtype=bool)
            for item in labels:  # Iterate over labels and collect all indices
                union = union + self[element+'.'+item.split('.')[-1]]
            ind = union
        elif mode in ['intersection']:
            intersect = sp.ones_like(self[element+'.all'], dtype=bool)
            for item in labels:  # Iterate over labels and collect all indices
                intersect = intersect*self[element+'.'+item.split('.')[-1]]
            ind = intersect
        elif mode in ['not_intersection']:
            not_intersect = sp.zeros_like(self[element+'.all'], dtype=int)
            for item in labels:  # Iterate over labels and collect all indices
                info = self[element+'.'+item.split('.')[-1]]
                not_intersect = not_intersect + sp.int8(info)
            ind = (not_intersect == 1)
        elif mode in ['not', 'difference']:
            none = sp.zeros_like(self[element+'.all'], dtype=int)
            for item in labels:  # Iterate over labels and collect all indices
                info = self[element+'.'+item.split('.')[-1]]
                none = none - sp.int8(info)
            ind = (none == 0)
        # Extract indices from boolean mask
        ind = sp.where(ind)[0]
        ind = ind.astype(dtype=int)
        return ind

    def pores(self, labels='all', mode='union'):
        r"""
        Returns pore locations where given labels exist, according to the logic
        specified by the mode argument.

        Parameters
        ----------
        labels : string or list of strings
            The label(s) whose pores locations are requested.  If omitted, all
            pore inidices are returned. This argument also accepts '*' for
            wildcard searches.

        mode : string
            Specifies how the query should be performed.  The options are:

            **'union'** : (default) All pores with ANY of the given labels are
            returned.

            **'intersection'** : Only pore with ALL the given labels are
            returned.

            **'not_intersection'** : Only pores with exactly one of the given
            labels are returned.

            **'not'** : Only pores with none of the given labels are returned.

        Returns
        -------
        A Numpy array containing pore indices where the specified label(s)
        exist.

        Examples
        --------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> pind = pn.pores(labels=['top', 'front'], mode='union')
        >>> pind[[0, 1, 2, -3, -2, -1]]
        array([  0,   5,  10, 122, 123, 124])
        >>> pn.pores(labels=['top', 'front'], mode='intersection')
        array([100, 105, 110, 115, 120])
        """
        ind = self._get_indices(element='pore', labels=labels, mode=mode)
        return ind

    @property
    def Ps(self):
        r"""
        A shortcut to get a list of all pores on the object
        """
        return sp.arange(0, self.Np)

    def throats(self, labels='all', mode='union'):
        r"""
        Returns throat locations where given labels exist, according to the
        logic specified by the mode argument.

        Parameters
        ----------
        labels : string or list of strings
            The throat label(s) whose locations are requested.  If omitted,
            'all' throat inidices are returned.  This argument also accepts
            '*' for wildcard searches.

        mode : string
            Specifies how the query should be performed.  The options are:

            **'union'** : (default) All throats with ANY of the given labels
            are returned.

            **'intersection'** : Only throats with ALL the given labels are
            counted.

            **'not_intersection'** : Only throats with exactly one of the
            given labels are counted.

            **'not'** : Only throats with none of the given labels are
            returned.

        Returns
        -------
        A Numpy array containing the throat indices where the specified
        label(s) exist.

        Examples
        --------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> Tind = pn.throats()
        >>> Tind[0:5]
        array([0, 1, 2, 3, 4])

        """
        ind = self._get_indices(element='throat', labels=labels, mode=mode)
        return ind

    @property
    def Ts(self):
        r"""
        A shortcut to get a list of all throats on the object
        """
        return sp.arange(0, self.Nt)

    def _tomask(self, locations, element):
        r"""
        This is a generalized version of tomask that accepts a string of
        'pore' or 'throat' for programmatic access.
        """
        element = self._parse_element(element, single=True)
        locations = self._parse_locations(locations)
        N = sp.shape(self[element + '.all'])[0]
        ind = sp.array(locations, ndmin=1)
        mask = sp.zeros((N, ), dtype=bool)
        mask[ind] = True
        return mask

    def tomask(self, pores=None, throats=None):
        r"""
        Convert a list of pore or throat indices into a boolean mask of the
        correct length

        Parameters
        ----------
        pores or throats : array_like
            List of pore or throat indices.  Only one of these can be specified
            at a time, and the returned result will be of the corresponding
            length.

        Returns
        -------
        A boolean mask of length Np or Nt with True in the specified pore or
        throat locations.

        Examples
        --------
        >>> import OpenPNM as op
        >>> pn = op.Network.Cubic(shape=[3, 3, 3])
        >>> mask = pn.tomask(pores=[0, 10, 20])
        >>> sum(mask)  # 3 non-zero elements exist in the mask (0, 10 and 20)
        3
        >>> len(mask)  # Mask size is equal to the number of pores in network
        27
        >>> mask = pn.tomask(throats=[0, 10, 20])
        >>> len(mask)  # Mask is now equal to number of throats in network
        54

        """
        if (pores is not None) and (throats is None):
            mask = self._tomask(element='pore', locations=pores)
        elif (throats is not None) and (pores is None):
            mask = self._tomask(element='throat', locations=throats)
        else:
            raise Exception('Cannot specify both pores and throats')
        return mask

    def toindices(self, mask):
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
        A list of pore or throat indices corresponding the locations where
        the received mask was True.

        Notes
        -----
        This behavior could just as easily be accomplished by using the mask
        in ``pn.pores()[mask]`` or ``pn.throats()[mask]``.  This method is just
        a thin convenience function and is a compliment to ``tomask``.

        """
        indices = self._parse_locations(mask)
        return indices

    def interpolate_data(self, data):
        r"""
        Determines a pore (or throat) property as the average of it's
        neighboring throats (or pores)

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
        - This uses an unweighted average, without attempting to account for
        distances or sizes of pores and throats.
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
            # Upcast data to full network size
            temp = sp.ones((net.Nt,))*sp.nan
            temp[Ts] = data
            data = temp
            temp = sp.ones((net.Np,))*sp.nan
            for pore in Ps:
                neighborTs = net.find_neighbor_throats(pore)
                neighborTs = net.filter_by_label(throats=neighborTs,
                                                 labels=label)
                temp[pore] = sp.mean(data[neighborTs])
            values = temp[Ps]
        elif sp.shape(data)[0] == self.Np:
            # Upcast data to full network size
            temp = sp.ones((net.Np, ))*sp.nan
            temp[Ps] = data
            data = temp
            Ps12 = net.find_connected_pores(throats=Ts, flatten=False)
            values = sp.mean(data[Ps12], axis=1)
        else:
            logger.error('Received data was an ambiguous length')
            raise Exception()
        return values

    def _interleave_data(self, prop, sources):
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
        when data is missing this can be tricky.  Data can be missing in two
        different ways: A set of pores is not assisgned to a geometry or the
        network contains multiple geometries and data does not exist on all.
        Float and boolean data is fine, but missing ints are converted to float
        when nans are inserted.

        Examples
        --------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> Ps = pn.pores('top',mode='not')
        >>> Ts = pn.find_neighbor_throats(pores=Ps,
        ...                               mode='intersection',
        ...                               flatten=True)
        >>> geom = OpenPNM.Geometry.TestGeometry(network=pn,
        ...                                      pores=Ps,
        ...                                      throats=Ts)
        >>> Ps = pn.pores('top')
        >>> Ts = pn.find_neighbor_throats(pores=Ps,
        ...                               mode='not_intersection')
        >>> boun = OpenPNM.Geometry.Boundary(network=pn, pores=Ps, throats=Ts)
        >>> geom['pore.test_float'] = sp.random.random(geom.Np)
        >>> print(sp.sum(~sp.isnan(pn['pore.test_float'])) == geom.Np)
        True
        >>> boun['pore.test_float'] = sp.random.random(boun.Np)
        >>> print(sp.sum(~sp.isnan(pn['pore.test_float'])) == pn.Np)
        True
        >>> geom['pore.test_int'] = sp.random.randint(0, 100, geom.Np)
        >>> print(pn['pore.test_int'].dtype.name.startswith('float'))
        True
        >>> boun['pore.test_int'] = sp.ones(boun.Np).astype(int)
        >>> print(pn['pore.test_int'].dtype.name.startswith('int'))
        True
        >>> geom['pore.test_bool'] = True
        >>> print(sp.sum(pn['pore.test_bool']) == geom.Np)
        True
        >>> boun['pore.test_bool'] = True
        >>> print(sp.sum(pn['pore.test_bool']) == pn.Np)
        True
        """
        element = self._parse_element(prop.split('.')[0], single=True)
        N = self._net._count(element)

        # Make sure sources contains objects, not just names of objects
        temp_sources = []
        for item in sources:
            # Check if sources were given as list of objects OR names
            try:
                item.name
            except:
                item = self._find_object(obj_name=item)
            temp_sources.append(item)
        sources = temp_sources

        # Attempt to fetch the requested prop array from each object
        arrs = [item.get(prop) for item in sources]
        locs = [item._net._get_indices(element, item.name) for item in sources]
        sizes = [sp.size(a) for a in arrs]
        if all([item is None for item in arrs]):  # prop not found anywhere
            raise KeyError(prop)
        if sp.any([i is None for i in arrs]):  # prop not found everywhere
            logger.warning('\''+prop+'\' not found on at least one object')

        # Check the general type of each array
        atype = []
        for a in arrs:
            if a is not None:
                t = a.dtype.name
                if t.startswith('int') or t.startswith('float'):
                    atype.append('numeric')
                elif t.startswith('bool'):
                    atype.append('boolean')
                else:
                    atype.append('other')
        if not all([item == atype[0] for item in atype]):
            raise Exception('The array types are not compatible')
        else:
            dummy_val = {'numeric': sp.nan, 'boolean': False, 'other': None}

        # Create an empty array of the right type and shape
        for item in arrs:
            if item is not None:
                if len(item.shape) == 1:
                    temp_arr = sp.zeros((N, ), dtype=item.dtype)
                else:
                    temp_arr = sp.zeros((N, item.shape[1]), dtype=item.dtype)
                temp_arr.fill(dummy_val[atype[0]])

        # Convrert int arrays to float IF NaNs are expected
        if (temp_arr.dtype.name.startswith('int') and
            (sp.any([i is None for i in arrs]) or
             sp.sum(sizes) != N)):
            temp_arr = temp_arr.astype(float)
            temp_arr.fill(sp.nan)

        # Fill new array with values in the corresponding locations
        for vals, inds in zip(arrs, locs):
            if vals is not None:
                temp_arr[inds] = vals
            else:
                temp_arr[inds] = dummy_val[atype[0]]
        return temp_arr

    def num_pores(self, labels='all', mode='union'):
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

            **'union'** : (default) All pores with ANY of the given labels are
            counted.

            **'intersection'** : Only pores with ALL the given labels are
            counted.

            **'not_intersection'** : Only pores with exactly one of the given
            labels are counted.

            **'difference'** : Only pores with none of the given labels are
            counted.

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
        >>> pn.num_pores(labels=['top', 'front'], mode='union')
        45
        >>> pn.num_pores(labels=['top', 'front'], mode='intersection')
        5
        >>> pn.num_pores(labels=['top', 'front'], mode='not_intersection')
        40

        """
        # Count number of pores of specified type
        Ps = self._get_indices(labels=labels, mode=mode, element='pore')
        Np = sp.shape(Ps)[0]
        return Np

    @property
    def Np(self):
        r"""
        A shortcut to query the total number of pores on the object'
        """
        return sp.shape(self.get('pore.all'))[0]

    def num_throats(self, labels='all', mode='union'):
        r"""
        Return the number of throats of the specified labels

        Parameters
        ----------
        labels : list of strings, optional
            The throat labels that should be included in the count.
            If not supplied, all throats are counted.

        mode : string, optional
            Specifies how the count should be performed.  The options are:

            **'union'** : (default) All throats with ANY of the given labels
            are counted.

            **'intersection'** : Only throats with ALL the given labels are
            counted.

            **'not_intersection'** : Only throats with exactly one of the given
            labels are counted.

            **'difference'** : Only throats with none of the given labels are
            counted.

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
        >>> pn.num_throats(labels=['top', 'front'], mode='union')
        76
        >>> pn.num_throats(labels=['top', 'front'], mode='intersection')
        4
        >>> pn.num_throats(labels=['top', 'front'], mode='not_intersection')
        72

        """
        # Count number of pores of specified type
        Ts = self._get_indices(labels=labels, mode=mode, element='throat')
        Nt = sp.shape(Ts)[0]
        return Nt

    @property
    def Nt(self):
        r"""
        A shortcut to query the total number of throats on the object'
        """
        return sp.shape(self.get('throat.all'))[0]

    def _count(self, element=None):
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
        element = self._parse_element(element=element)
        temp = {elem: sp.size(self[elem+'.all']) for elem in element}
        # TODO: In a future version this should always just return a dict
        if len(temp) == 1:
            temp = list(temp.values())[0]
        return temp

    def _map(self, element, locations, target, return_mapping=False):
        r"""
        """
        # Initialize things
        locations = self._parse_locations(locations=locations)
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
        >>> Ps = pn.pores(labels=['top', 'left'], mode='intersection')
        >>> Ps
        array([100, 101, 102, 103, 104])
        >>> geom = OpenPNM.Geometry.GenericGeometry(network=pn, pores=Ps)
        >>> geom.Ps
        array([0, 1, 2, 3, 4])
        >>> geom.map_pores(target=pn, pores=geom.Ps)
        array([100, 101, 102, 103, 104])
        >>> pn.map_pores(target=geom, pores=Ps)
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
        >>> Ts = pn.throats(labels=['top', 'left'], mode='intersection')
        >>> Ts
        array([260, 262, 264, 266])
        >>> geom = OpenPNM.Geometry.GenericGeometry(network=pn, throats=Ts)
        >>> geom.Ts
        array([0, 1, 2, 3])
        >>> geom.map_throats(target=pn, throats=geom.Ts)
        array([260, 262, 264, 266])
        >>> pn.map_throats(target=geom, throats=Ts)
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

    @property
    def Tnet(self):
        r"""
        A shortcut to retrieve the mapping of the current object's throats onto
        the network.
        """
        return self.map_throats()

    @property
    def Pnet(self):
        r"""
        A shortcut to retrieve the mapping of the current object's pores onto
        the network.
        """
        return self.map_pores()

    def _parse_locations(self, locations):
        r"""
        This private method accepts a list of pores or throats and returns a
        properly structured Numpy array of indices.

        Parameters
        ----------
        locations : multiple options
            This argument can accept numerous different data types including
            boolean masks, integers and arrays.

        Returns
        -------
        A Numpy array of locations.

        Notes
        -----
        This method should only be called by the method that is actually using
        the locations, to avoid calling it multiple times.
        """
        if locations is None:
            locations = sp.array([], ndmin=1, dtype=int)
        locs = sp.array(locations, ndmin=1)
        if locs.dtype not in [bool, sp.float16, sp.float32, sp.float64,
                              sp.int8, sp.int16, sp.int32, sp.int64]:
            raise Exception('Invalid data type received as locations, ' +
                            ' must be boolean or numeric')
        if locs.dtype == bool:
            if sp.size(locs) == self.Np:
                locs = self.Ps[locs]
            elif sp.size(locs) == self.Nt:
                locs = self.Ts[locs]
            else:
                raise Exception('Boolean list of locations must be either ' +
                                'Np nor Nt long')
        locs = locs.astype(dtype=int)
        return locs

    def _parse_element(self, element, single=False):
        r"""
        This private method is used to parse the keyword \'element\' in many
        of the above methods.

        Parameters
        ----------
        element : string or list of strings
            The element argument to check.  If is None is recieved, then a list
            containing both \'pore\' and \'throat\' is returned.

        single : boolean (default is False)
            When set to True only a single element is allowed and it will also
            return a string containing the element.

        Returns
        -------
        When ``single`` is False (default) a list contain the element(s) is
        returned.  When ``single`` is True a bare string containing the element
        is returned.
        """
        if element is None:
            element = ['pore', 'throat']
        # Convert element to a list for subsequent processing
        if type(element) is str:
            element = [element]
        # Convert 'pore.prop' and 'throat.prop' into just 'pore' and 'throat'
        element = [item.split('.')[0] for item in element]
        # Make sure all are lowercase
        element = [item.lower() for item in element]
        # Deal with an plurals
        element = [item.rsplit('s', maxsplit=1)[0] for item in element]
        for item in element:
            if item not in ['pore', 'throat']:
                raise Exception('Invalid element received')
        # Remove duplicates if any
        [element.remove(L) for L in element if element.count(L) > 1]
        if single:
            if len(element) > 1:
                raise Exception('Both elements recieved when single element ' +
                                'allowed')
            else:
                element = element[0]
        return element

    def _parse_labels(self, labels, element):
        r"""
        This private method is used for converting \'labels\' to a proper
        format, including dealing with wildcards (\*).

        Parameters
        ----------
        labels : string or list of strings
            The label or list of labels to be parsed. Note that the \* can be
            used as a wildcard.

        Returns
        -------
        A list of label strings, with all wildcard matches included if
        applicable.
        """
        if labels is None:
            raise Exception('Labels cannot be None')
        if type(labels) is str:
            labels = [labels]
        # Parse the labels list
        parsed_labels = []
        for label in labels:
            # Remove element from label, if present
            if element in label:
                label = label.split('.')[-1]
            # Deal with wildcards
            if '*' in label:
                Ls = [L.split('.')[-1] for L in self.labels(element=element)]
                if label.startswith('*'):
                    temp = [L for L in Ls if L.endswith(label.strip('*'))]
                if label.endswith('*'):
                    temp = [L for L in Ls if L.startswith(label.strip('*'))]
                temp = [element+'.'+L for L in temp]
            elif element+'.'+label in self.keys():
                temp = [element+'.'+label]
            else:
                # TODO: The following Error should/could be raised but it
                # breaks the net-geom and phase-phys look-up logic
                # raise KeyError('\''+element+'.'+label+'\''+' not found')
                logger.warning('\''+element+'.'+label+'\''+' not found')
                temp = [element+'.'+label]
            parsed_labels.extend(temp)
            # Remove duplicates if any
            [parsed_labels.remove(L) for L in parsed_labels
             if parsed_labels.count(L) > 1]
        return parsed_labels

    def _parse_mode(self, mode, allowed=None, single=False):
        r"""
        This private method is for checking the \'mode\' used in the calling
        method.

        Parameters
        ----------
        mode : string or list of strings
            The mode(s) to be parsed

        allowed : list of strings
            A list containing the allowed modes.  This list is defined by the
            calling method.  If any of the received modes are not in the
            allowed list an exception is raised.

        single : boolean (default is False)
            Indicates if only a single mode is allowed.  If this argument is
            True than a string is returned rather than a list of strings, which
            makes it easier to work with in the caller method.

        Returns
        -------
        A list containing the received modes as strings, checked to ensure they
        are all within the allowed set (if provoided).  Also, if the ``single``
        argument was True, then a string is returned.
        """
        if type(mode) is str:
            mode = [mode]
        for item in mode:
            if (allowed is not None) and (item not in allowed):
                raise Exception('\'mode\' must be one of the following: ' +
                                allowed.__str__())
        # Remove duplicates, if any
        [mode.remove(L) for L in mode if mode.count(L) > 1]
        if single:
            if len(mode) > 1:
                raise Exception('Multiple modes received when only one mode ' +
                                'allowed')
            else:
                mode = mode[0]
        return mode

    def _isa(self, keyword=None, obj=None):
        r"""
        """
        if keyword is None:
            mro = [item.__name__ for item in self.__class__.__mro__]
        if obj is None:
            query = False
            mro = [item.__name__ for item in self.__class__.__mro__]
            if 'net' in keyword.lower():
                if 'GenericNetwork' in mro:
                    query = True
            elif 'geo' in keyword.lower():
                if 'GenericGeometry' in mro:
                    query = True
            elif 'phase' in keyword.lower():
                if 'GenericPhase' in mro:
                    query = True
            elif 'phys' in keyword.lower():
                if 'GenericPhysics' in mro:
                    query = True
            elif 'alg' in keyword.lower():
                if 'GenericAlgorithm' in mro:
                    query = True
            elif 'clone' in keyword.lower():
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

    def check_data_health(self, props=[], element=None):
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
            prop = item
            required = self._count(item.split('.')[0])
            if len(prop) > 35:  # Trim overly long prop names
                prop = prop[0:32] + '...'
            if self[item].dtype == object:  # Print objects differently
                invalid = [i for i in self[item] if i is None]
                defined = sp.size(self[item]) - len(invalid)
                lines.append("{0:<5d} {1:<35s} {2:>5d} / {3:<5d}".format(i + 1,
                                                                         prop,
                                                                         defined,
                                                                         required))
            elif '._' not in prop:
                a = sp.isnan(self[item])
                defined = sp.shape(self[item])[0] \
                    - a.sum(axis=0, keepdims=(a.ndim-1) == 0)[0]
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
            if '._' not in prop:
                lines.append("{0:<5d} {1:<35s} {2:<10d}".format(i + 1,
                                                                prop,
                                                                sp.sum(self[item])))
        lines.append(horizonal_rule)
        return '\n'.join(lines)
