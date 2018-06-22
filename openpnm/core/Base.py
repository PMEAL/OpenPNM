from collections import namedtuple
import matplotlib.pyplot as plt
from openpnm.core import Workspace, logging
from openpnm.utils.misc import PrintableList, SettingsDict, HealthDict
import scipy as sp
logger = logging.getLogger(__name__)
ws = Workspace()


class Base(dict):
    r"""
    Contains methods for working with the data in the OpenPNM dict objects
    """

    def __new__(cls, *args, **kwargs):
        instance = super(Base, cls).__new__(cls, *args, **kwargs)
        # The SettingsDict implements the __missing__ magic method, which
        # returns None instead of KeyError.  This is useful for checking the
        # value of a settings without first ensuring it exists.
        instance.settings = SettingsDict()
        return instance

    def __init__(self, Np=0, Nt=0, name=None, project=None):
        self.settings.setdefault('prefix', 'base')
        super().__init__()
        if project is None:
            project = ws.new_project()
        project.extend(self)
        self.name = name
        self.update({'pore.all': sp.ones(shape=(Np, ), dtype=bool)})
        self.update({'throat.all': sp.ones(shape=(Nt, ), dtype=bool)})

    def __repr__(self):
        return '<%s object at %s>' % (self.__class__.__module__, hex(id(self)))

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
        long arrays can be written, and they must be called 'pore.***' or
        'throat.***'.  Also, any scalars are cast into full length vectors.

        """
        # If value is a dictionary, then break it up into constituent arrays
        if hasattr(value, 'keys'):
            for item in value.keys():
                prop = item.replace('pore.', '').replace('throat.', '')
                self.__setitem__(key+'.'+prop, value[item])
            return

        value = sp.array(value, ndmin=1)  # Convert value to an ndarray

        # Enforce correct dict naming
        element = key.split('.')[0]
        if element == 'constant':
            super(Base, self).__setitem__(key, value)
            return
        else:
            element = self._parse_element(element, single=True)

        # Skip checks for 'coords', 'conns'
        if key in ['pore.coords', 'throat.conns']:
            super(Base, self).__setitem__(key, value)
            return

        # Skip checks for protected props, and prevent changes if defined
        protected_keys = ['all']
        if key.split('.')[1] in protected_keys:
            if key in self.keys():
                if sp.shape(self[key]) == (0, ):
                    super(Base, self).__setitem__(key, value)
                else:
                    logger.warning(key+' is already defined.')
            else:
                super(Base, self).__setitem__(key, value)
            return

        # Write value to dictionary
        if sp.shape(value)[0] == 1:  # If value is scalar
            value = sp.ones((self._count(element), ), dtype=value.dtype)*value
            super(Base, self).__setitem__(key, value)
        elif sp.shape(value)[0] == self._count(element):
            super(Base, self).__setitem__(key, value)
        else:
            if self._count(element) == 0:
                self.update({key: value})
            else:
                raise Exception('Cannot write array, wrong length: '+key)

    def _set_name(self, name, validate=True):
        if not hasattr(self, '_name'):
            self._name = None
        if name is None:
            name = self.project._generate_name(self)
        if self.name == name:
            return
        if validate:
            self.project._validate_name(name)
        if self._name is not None:
            # Rename any label arrays in other objects
            for item in self.project:
                if 'pore.'+self.name in item.keys():
                    item['pore.'+name] = item.pop('pore.'+self.name)
                if 'throat.'+self.name in item.keys():
                    item['throat.'+name] = item.pop('throat.'+self.name)
        self._name = name

    def _get_name(self):
        if not hasattr(self, '_name'):
            self._name = None
        return self._name

    name = property(_get_name, _set_name)

    def _get_project(self):
        for proj in ws.values():
            if self in proj:
                return proj

    project = property(fget=_get_project)

    def clear(self, element=None, mode='all'):
        r"""
        A subclassed version of the standard dict's clear method.  This can be
        used to selectively clear certain data from the object, including
        properties and/or labels.  Importantly, it does NOT clear items that
        are required to maintain the integrity of the simulation.  These are
        arrays that define the topology (ie. 'pore.all', 'pore.coords',
        'throat.all', 'throat.conns'), as well as arrays that indicate

        Parameters
        ----------
        element : string or list of strings
            Can be either 'pore' or 'throat', which specifies whether 'pore'
            and/or 'throat' arrays should be cleared.  The default is both.

        mode : string or list of strings
            This controls what is cleared from the object.  Options are:

            **'props'** : Removes all numerical property values from the object
            dictionary.

            **'model_data'** : Removes only numerical data that were produced
            by an associated model.

            **'labels'** : Removes all labels from the object dictionary,
            except those relating to the pore and throat locations of
            associated objects.

            **'all'** : Removes both 'props' and 'labels'

        Notes
        -----
        If you wish to selectively remove some properties but not others, use
        something like ``del object['pore.blah']`` at the python prompt.
        Obviously this can be done in a for-loop to remove a list of items.

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> len(pn.labels())  # There are 10 total labels on the network
        10
        >>> pn.clear(mode='labels')
        >>> len(pn.labels())  # Kept only 'pore.all' and 'throat.all'
        2
        >>> geom = op.geometry.GenericGeometry(network=pn, pores=pn.Ps,
        ...                                    throats=pn.Ts, name='geo1')
        >>> len(pn.labels())  # 2 new labels were added for geometry locations
        4
        >>> pn.clear(mode='labels')
        >>> 'pore.'+geom.name in pn.keys()  # The geometry labels were kept
        True
        >>> len(pn.props())  # The network has two properties
        2
        >>> pn.clear(element='pore', mode='props')
        >>> 'pore.coords' in pn.keys()  # The pore property was removed
        True
        >>> pn.clear()  # Remove everything except protected labels and arrays
        >>> print(sorted(list(pn.keys(element='pore', mode='all'))))
        ['pore.all', 'pore.coords', 'pore.geo1']

        """
        protected = ['pore.all', 'throat.all', 'pore.coords', 'throat.conns']
        allowed = ['props', 'labels', 'model_data', 'all']
        mode = self._parse_mode(mode=mode, allowed=allowed)
        if 'model_data' in mode:
            for item in self.models.keys():
                print('deleting ' + item)
                del self[item]
            mode.remove('model_data')
        for item in self.keys(mode=mode, element=element):
            if item not in protected:
                if item.split('.')[1] not in self.project.names:
                    del self[item]

    def keys(self, element=None, mode=None):
        r"""
        This subclass works exactly like ``keys`` when no arguments are passed,
        but optionally accepts an ``element`` and/or a ``mode``, which filters
        the output to only the requested keys.

        The default behavior is exactly equivalent to the normal ``keys``
        method.

        Parameters
        ----------
        element : string
            Can be either 'pore' or 'throat', which limits the returned list of
            keys to only 'pore' or 'throat' keys.  If neither is given, then
            both are assumed.

        mode : string (optional, default is 'skip')
            Controls which keys are returned.  Options are:

            **``None``** : This mode (default) bypasses this subclassed method
            and just returns the normal KeysView object.

            **'labels'** : Limits the returned list of keys to only 'labels'
            (boolean arrays)

            **'props'** : Limits he return list of keys to only 'props'
            (numerical arrays).

            **'all'** : Returns both 'labels' and 'props'.  This is equivalent
            to sending a list of both 'labels' and 'props'.

        See Also
        --------
        props
        labels

        Notes
        -----
        This subclass can be used to get dictionary keys of specific kinds of
        data.  It's use augments ``props`` and ``labels`` by returning a list
        containing both types, but possibly limited by element type ('pores'
        or 'throats'.)

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic([5, 5, 5])
        >>> pn.keys(mode='props')  # Get all props
        ['pore.coords', 'throat.conns']
        >>> pn.keys(mode='props', element='pore')  # Get only pore props
        ['pore.coords']

        """
        if mode is None:
            return super().keys()
        element = self._parse_element(element=element)
        allowed = ['props', 'labels']
        if 'all' in mode:
            mode = allowed
        mode = self._parse_mode(mode=mode, allowed=allowed)
        keys = super().keys()
        temp = []
        if 'props' in mode:
            temp.extend([i for i in keys if self.get(i).dtype != bool])
        if 'labels' in mode:
            temp.extend([i for i in keys if self.get(i).dtype == bool])
        if element:
            temp = [i for i in temp if i.split('.')[0] in element]
        return temp

    # -------------------------------------------------------------------------
    """Data Query Methods"""
    # -------------------------------------------------------------------------
    def props(self, element=None, mode='all', deep=False):
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

            **'constants'** : returns data values that were *not* generated by
            a model, but manaully created.

        Returns
        -------
        A an alphabetically sorted list containing the string name of all
        pore or throat properties currently defined.  This list is an iterable,
        so is useful for scanning through properties.

        See Also
        --------
        labels
        keys

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[3, 3, 3])
        >>> pn.props('pore')
        ['pore.coords']
        >>> pn.props('throat')
        ['throat.conns']
        >>> pn.props()
        ['pore.coords', 'throat.conns']
        """
        # Parse Inputs
        element = self._parse_element(element=element)
        allowed_modes = ['all', 'constants', 'models']
        mode = self._parse_mode(mode=mode, allowed=allowed_modes, single=True)
        if mode == 'all':
            vals = set(self.keys(mode='props'))
        if mode == 'constants':
            if hasattr(self, 'models'):
                temp = set(self.keys(mode='props'))
                vals = temp.difference(self.models.keys())
            else:
                vals = set(self.keys(mode='props'))
        if mode == 'models':
            if hasattr(self, 'models'):
                temp = set(self.keys(mode='props'))
                vals = temp.intersection(self.models.keys())
            else:
                logger.warning('Object does not have a models attribute')
                vals = set()
        # Deal with hidden props
        hide = set([i for i in self.keys() if i.split('.')[1].startswith('_')])
        vals = vals.difference(hide)
        # Remove values of the wrong element
        temp = set([i for i in vals if i.split('.')[0] not in element])
        vals = set(vals).difference(temp)
        # Convert to nice list for printing
        vals = PrintableList(list(vals))
        return vals

    def _get_labels(self, element, locations, mode):
        r"""
        This is the actual label getter method, but it should not be called
        directly.  Wrapper methods have been created, use ``labels``.
        """
        # Parse inputs
        locations = self._parse_indices(locations)
        allowed = ['none', 'union', 'intersection', 'complement', 'not',
                   'count', 'mask', 'difference']
        mode = self._parse_mode(mode=mode, allowed=allowed, single=True)
        element = self._parse_element(element=element)
        # Collect list of all pore OR throat labels
        labels = self.keys(mode='labels', element=element)
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
            return PrintableList(temp)
        if mode in ['intersection']:
            temp = labels[sp.sum(arr, axis=0) == sp.shape(locations, )[0]]
            temp.tolist()
            return PrintableList(temp)
        if mode in ['not', 'complement', 'difference']:
            temp = labels[sp.sum(arr, axis=0) != sp.shape(locations, )[0]]
            temp.tolist()
            return PrintableList(temp)
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

            **'complement'** : Labels NOT applied to ALL pores (or throats)

            **'count'** : The number of labels on each pores (or throats)

            **'mask'**: returns an N x Lt array, where each row corresponds to
            a pore (or throat) location, and each column contains the truth
            value for the existance of labels as returned from
            ``labels(element='pores')``.

        Returns
        -------
        A list containing the dictionary keys on the object, limited by the
        specified ``mode``.

        See Also
        --------
        props
        keys

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> pn.labels(pores=[0, 1, 5, 6])
        ['pore.all', 'pore.bottom', 'pore.front', 'pore.internal', 'pore.left']
        >>> pn.labels(pores=[0, 1, 5, 6], mode='intersection')
        ['pore.all', 'pore.front', 'pore.internal']
        """
        labels = PrintableList()
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

    def _get_indices(self, element, labels='all', mode='union'):
        r"""
        This is the actual method for getting indices, but should not be called
        directly.  Use pores or throats instead.
        """
        # Parse and validate all input values.
        allowed = ['union', 'intersection', 'not_intersection', 'not',
                   'difference', 'complement']
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
            none = sp.zeros_like(self[element+'.all'], dtype=int)
            for item in labels:  # Iterate over labels and collect all indices
                info = self[element+'.'+item.split('.')[-1]]
                none = none + sp.int8(info)
            ind = (none == 1)
        elif mode in ['not', 'difference', 'complement']:
            none = sp.zeros_like(self[element+'.all'], dtype=int)
            for item in labels:  # Iterate over labels and collect all indices
                info = self[element+'.'+item.split('.')[-1]]
                none = none + sp.int8(info)
            ind = (none == 0)
        # Extract indices from boolean mask
        ind = sp.where(ind)[0]
        ind = ind.astype(dtype=int)
        return ind

    def pores(self, labels='all', mode='union'):
        r"""
        Returns pore indicies where given labels exist, according to the logic
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
            returned. This is equivalent to ``intersection`` which is
            deprecated.

            **'exclusive_or'** : Only pores with exactly ONE of the given
            labels are returned.

            **'complement'** : Only pores with NONE of the given labels are
            returned. This is equivalent to ``not_intersection`` which is
            deprecated.

        Returns
        -------
        A Numpy array containing pore indices where the specified label(s)
        exist.

        See Also
        --------
        throats

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> Ps = pn.pores(labels=['top', 'front'], mode='union')
        >>> Ps[:5]  # Look at first 5 pore indices
        array([0, 1, 2, 3, 4])
        >>> pn.pores(labels=['top', 'front'], mode='intersection')
        array([ 4,  9, 14, 19, 24])
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

            **'all'** : Only throats with ALL the given labels are
            counted.

            **'one'** : Only throats with exactly ONE of the
            given labels are counted.

            **'none'** : Only throats with NONE of the given labels are
            returned.

        Returns
        -------
        A Numpy array containing the throat indices where the specified
        label(s) exist.

        See Also
        --------
        pores

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[3, 3, 3])
        >>> Ts = pn.throats()
        >>> Ts[0:5]  # Look at first 5 throat indices
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

    def _map(self, ids, element, filtered):
        locations = self._get_indices(element=element)
        hash_map = dict(zip(self[element+'._id'], locations))
        ind = sp.array([hash_map.get(i, -1) for i in ids])
        mask = sp.zeros(shape=ids.shape, dtype=bool)
        mask[sp.where(ind >= 0)[0]] = True
        if filtered:
            return ind[mask]
        else:
            t = namedtuple('index_map', ('indices', 'mask'))
            return t(ind, mask)

    def map_pores(self, ids, filtered=True):
        r"""
        Translates pore ids into indices

        Parameters
        ----------
        ids : array_like
            The ids of the pores whose indices are sought

        filtered : boolean (default is ``True``)
            If ``True`` then a ND-array of indices is returned, otherwise
            a named-tuple containing the ``indices`` and the ???
        """
        return self._map(element='pore', ids=ids, filtered=filtered)

    def map_throats(self, ids, filtered=True):
        r"""
        Translates ids to indices

        Parameters
        ----------
        throats : array_like
            The throat indices for which full network indices are sought

        filtered : boolean (default is ``True``)
            If ``True`` then a ND-array of indices is returned, otherwise
            a named-tuple containing the ``indices`` and the ???
        """
        return self._map(element='throat', ids=ids, filtered=filtered)

    def _tomask(self, indices, element):
        r"""
        This is a generalized version of tomask that accepts a string of
        'pore' or 'throat' for programmatic access.
        """
        element = self._parse_element(element, single=True)
        indices = self._parse_indices(indices)
        N = sp.shape(self[element + '.all'])[0]
        ind = sp.array(indices, ndmin=1)
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

        See Also
        --------
        toindices

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> mask = pn.tomask(pores=[0, 10, 20])
        >>> sum(mask)  # 3 non-zero elements exist in the mask (0, 10 and 20)
        3
        >>> len(mask)  # Mask size is equal to the number of pores in network
        125
        >>> mask = pn.tomask(throats=[0, 10, 20])
        >>> len(mask)  # Mask is now equal to number of throats in network
        300

        """
        if (pores is not None) and (throats is None):
            mask = self._tomask(element='pore', indices=pores)
        elif (throats is not None) and (pores is None):
            mask = self._tomask(element='throat', indices=throats)
        else:
            raise Exception('Cannot specify both pores and throats')
        return mask

    def toindices(self, mask):
        r"""
        Convert a boolean mask to a list of pore or throat indices

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

        See Also
        --------
        tomask

        Notes
        -----
        This behavior could just as easily be accomplished by using the mask
        in ``pn.pores()[mask]`` or ``pn.throats()[mask]``.  This method is
        just a convenience function and is a compliment to ``tomask``.

        """
        if sp.amax(mask) > 1:
            raise Exception('Recieved mask is invalid, with values above 1')
        mask = sp.array(mask, dtype=bool)
        indices = self._parse_indices(mask)
        return indices

    def _interleave_data(self, prop, sources):
        r"""
        Retrieves requested property from associated objects, to produce a full
        Np or Nt length array.

        Parameters
        ----------
        prop : string
            The property name to be retrieved

        Returns
        -------
        A full length (Np or Nt) array of requested property values.

        Notes
        -----
        This makes an effort to maintain the data 'type' when possible; however
        when data are missing this can be tricky.  Data can be missing in two
        different ways: A set of pores is not assisgned to a geometry or the
        network contains multiple geometries and data does not exist on all.
        Float and boolean data is fine, but missing ints are converted to float
        when nans are inserted.

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[2, 2, 2])
        >>> Ps = pn['pore.top']
        >>> Ts = pn.find_neighbor_throats(pores=Ps)
        >>> g1 = op.geometry.GenericGeometry(network=pn, pores=Ps, throats=Ts)
        >>> Ts = ~pn.tomask(throats=Ts)
        >>> g2 = op.geometry.GenericGeometry(network=pn, pores=~Ps, throats=Ts)
        >>> g1['pore.value'] = 1
        >>> print(g1['pore.value'])
        [1 1 1 1]
        >>> print(g2['pore.value'])  # 'pore.value' is defined on g1, not g2
        [nan nan nan nan]
        >>> print(pn['pore.value'])
        [nan  1. nan  1. nan  1. nan  1.]
        >>> g2['pore.value'] = 20
        >>> print(pn['pore.value'])
        [20  1 20  1 20  1 20  1]
        >>> pn['pore.label'] = False
        >>> print(g1['pore.label'])  # 'pore.label' is defined on pn, not g1
        [False False False False]
        """
        element = self._parse_element(prop.split('.')[0], single=True)
        N = self.project.network._count(element)

        # Fetch sources list depending on object type?
        # proj = self.project
        # if self._isa('network'):
        #     sources = list(proj.geometries().values())
        # elif self._isa('phase'):
        #     sources = list(proj.phases().values())
        # elif self._isa('physics')
        #     sources = list(proj.physics().values())
        # elif self._isa('physics'):
        #     sources = list(proj.geometries().values())
        # else:
        #     pass

        # Attempt to 'get' the requested array from each object
        # Use 'get' so that missing keys return None, instead of KeyError
        arrs = [item.get(prop, None) for item in sources]
        locs = [self._get_indices(element, item.name) for item in sources]
        sizes = [sp.size(a) for a in arrs]
        if sp.all([item is None for item in arrs]):  # prop not found anywhere
            raise KeyError(prop)
#        if sp.any([i is None for i in arrs]):  # prop not found everywhere
#            logger.warning('\''+prop+'\' not found on at least one object')
#        if sp.sum(sizes) < self._count(element):
#            logger.warning('Not all '+element+'s are assigned to an object')
#            N_missing = self._count(element) - sp.sum(sizes)
#            arrs.append(sp.zeros(shape=(N_missing,), dtype=float)*sp.nan)

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

        # Convert int arrays to float IF NaNs are expected
        if temp_arr.dtype.name.startswith('int') and \
           (sp.any([i is None for i in arrs]) or sp.sum(sizes) != N):
            temp_arr = temp_arr.astype(float)
            temp_arr.fill(sp.nan)

        # Fill new array with values in the corresponding locations
        for vals, inds in zip(arrs, locs):
            if vals is not None:
                temp_arr[inds] = vals
            else:
                temp_arr[inds] = dummy_val[atype[0]]
        return temp_arr

    def interpolate_data(self, propname):
        r"""
        Determines a pore (or throat) property as the average of it's
        neighboring throats (or pores)

        Parameters
        ----------
        propname: string
            The dictionary key to the values to be interpolated.

        Returns
        -------
        An array containing interpolated pore (or throat) data

        Notes
        -----
        - This uses an unweighted average, without attempting to account for
        distances or sizes of pores and throats.

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[3, 1, 1])
        >>> pn['pore.value'] = [1, 2, 3]
        >>> pn.interpolate_data('pore.value')
        array([1.5, 2.5])
        """
        mro = self._mro()
        if 'GenericNetwork' in mro:
            net = self
            Ts = net.throats()
            Ps = net.pores()
            label = 'all'
        elif ('GenericPhase' in mro) or ('GenericAlgorithm' in mro):
            net = self.project.network
            Ts = net.throats()
            Ps = net.pores()
            label = 'all'
        elif ('GenericGeometry' in mro) or ('GenericPhysics' in mro):
            net = self.project.network
            Ts = net.throats(self.name)
            Ps = net.pores(self.name)
            label = self.name
        if propname.startswith('throat'):
            # Upcast data to full network size
            temp = sp.ones((net.Nt,))*sp.nan
            temp[Ts] = self[propname]
            data = temp
            temp = sp.ones((net.Np,))*sp.nan
            for pore in Ps:
                neighborTs = net.find_neighbor_throats(pore)
                neighborTs = net.filter_by_label(throats=neighborTs,
                                                 labels=label)
                temp[pore] = sp.mean(data[neighborTs])
            values = temp[Ps]
        elif propname.startswith('pore'):
            # Upcast data to full network size
            data = sp.ones((net.Np, ))*sp.nan
            data[Ps] = self[propname]
            Ps12 = net['throat.conns'][Ts]
            values = sp.mean(data[Ps12], axis=1)
        return values

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

            **'exlusive_or'** : Only locations with exactly one of the
            given labels are kept.

            **'complement'** : Only locations with none of the given labels
            are kept.

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
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> pn.filter_by_label(pores=[0, 1, 5, 6], labels='left')
        array([0, 1])
        >>> Ps = pn.pores(['top', 'bottom', 'front'], mode='union')
        >>> pn.filter_by_label(pores=Ps, labels=['top', 'front'],
        ...                    mode='intersection')
        array([ 4,  9, 14, 19, 24])
        """
        allowed = ['union', 'intersection', 'exclusive_or', 'complement']
        mode = self._parse_mode(mode=mode, allowed=allowed, single=True)
        # Convert inputs to locations and element
        if (sp.size(throats) > 0) and (sp.size(pores) > 0):
            raise Exception('Can only filter either pores OR labels per call')
        if sp.size(pores) > 0:
            element = 'pore'
            locations = self._parse_indices(pores)
        elif sp.size(throats) > 0:
            element = 'throat'
            locations = self._parse_indices(throats)
        else:
            return(sp.array([], dtype=int))
        # Do it
        labels = self._parse_labels(labels=labels, element=element)
        labels = [element+'.'+item.split('.')[-1] for item in labels]
        all_locs = self._get_indices(element=element, labels=labels, mode=mode)
        mask = self._tomask(indices=all_locs, element=element)
        ind = mask[locations]
        return locations[ind]

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

            **intersection** : Only pores with ALL the given labels are
            counted.

            **'exclusive_or'** : Only pores with exactly one of the given
            labels are counted.

            **'complement'** : Only pores with none of the given labels are
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
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> pn.num_pores()
        125
        >>> pn.num_pores(labels=['top'])
        25
        >>> pn.num_pores(labels=['top', 'front'], mode='union')
        45
        >>> pn.num_pores(labels=['top', 'front'], mode='intersection')
        5

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

            **'exclusive_or'** : Only throats with exactly one of the given
            labels are counted.

            **'complement'** : Only throats with none of the given labels are
            counted.

        Returns
        -------
        Nt : int
            Number of throats with the specified labels

        See Also
        --------
        num_pores
        count

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
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> pn._count('pore')
        125
        >>> pn._count('throat')
        300
        """
        element = self._parse_element(element=element, single=True)
        temp = sp.size(super(Base, self).__getitem__(element+'.all'))
        return temp

    def show_hist(self, props=[], bins=20, **kwargs):
        r"""
        Show a quick plot of key property distributions.

        Parameters
        ----------
        props : string or list of strings
            The pore and/or throat properties to be plotted as histograms

        bins : int or array_like
            The number of bins to use when generating the histogram.  If an
            array is given they are used as the bin spacing instead.

        Notes
        -----
        Other keyword arguments are passed to the ``matplotlib.pyplot.hist``
        function.
        """
        if type(props) is str:
            props = [props]
        N = len(props)
        if N == 1:
            r = 1
            c = 1
        elif N < 4:
            r = 1
            c = N
        else:
            r = int(sp.ceil(N**0.5))
            c = int(sp.floor(N**0.5))

        for i in range(len(props)):
            plt.subplot(r, c, i+1)
            plt.hist(self[props[i]], bins=bins, **kwargs)

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
        >>> import openpnm
        >>> pn = openpnm.network.Cubic(shape=[5, 5, 5])
        >>> h = pn.check_data_health()
        >>> h.health
        True
        """
        health = HealthDict()
        if props == []:
            props = self.props(element)
        else:
            if type(props) == str:
                props = [props]
        for item in props:
            health[item] = []
            if sp.sum(sp.isnan(self[item])) > 0:
                health[item] = 'Has NaNs'
            elif sp.shape(self[item])[0] != self._count(item.split('.')[0]):
                health[item] = 'Wrong Length'
        return health

    def _parse_indices(self, indices):
        r"""
        This private method accepts a list of pores or throats and returns a
        properly structured Numpy array of indices.

        Parameters
        ----------
        indices : multiple options
            This argument can accept numerous different data types including
            boolean masks, integers and arrays.

        Returns
        -------
        A Numpy array of indices.

        Notes
        -----
        This method should only be called by the method that is actually using
        the locations, to avoid calling it multiple times.
        """
        if indices is None:
            indices = sp.array([], ndmin=1, dtype=int)
        locs = sp.array(indices, ndmin=1)
        # If boolean array, convert to indices
        if locs.dtype == bool:
            if sp.size(locs) == self.Np:
                locs = self.Ps[locs]
            elif sp.size(locs) == self.Nt:
                locs = self.Ts[locs]
            else:
                raise Exception('Mask of locations must be either ' +
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
                raise Exception('Invalid element received: '+item)
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

    def _parse_prop(self, propname, element):
        r"""
        """
        element = self._parse_element(element, single=True)
        return element + '.' + propname.split('.')[-1]

    def __str__(self):
        horizonal_rule = '―' * 78
        lines = [horizonal_rule]
        lines.append(self.__module__.replace('__', '') + ': \t' + self.name)
        lines.append(horizonal_rule)
        lines.append("{0:<5s} {1:<45s} {2:<10s}".format('#',
                                                        'Properties',
                                                        'Valid Values'))
        lines.append(horizonal_rule)
        props = list(set(self.keys()).difference(set(self.labels())))
        props.sort()
        for i, item in enumerate(props):
            prop = item
            required = self._count(item.split('.')[0])
            if len(prop) > 35:  # Trim overly long prop names
                prop = prop[0:32] + '...'
            if self[item].dtype == object:  # Print objects differently
                invalid = [i for i in self[item] if i is None]
                defined = sp.size(self[item]) - len(invalid)
                fmt = "{0:<5d} {1:<45s} {2:>5d} / {3:<5d}"
                lines.append(fmt.format(i + 1, prop, defined, required))
            elif '._' not in prop:
                a = sp.isnan(self[item])
                defined = sp.shape(self[item])[0] \
                    - a.sum(axis=0, keepdims=(a.ndim-1) == 0)[0]
                fmt = "{0:<5d} {1:<45s} {2:>5d} / {3:<5d}"
                lines.append(fmt.format(i + 1, prop, defined, required))
        lines.append(horizonal_rule)
        lines.append("{0:<5s} {1:<45s} {2:<10s}".format('#',
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
                fmt = "{0:<5d} {1:<45s} {2:<10d}"
                lines.append(fmt.format(i + 1, prop, sp.sum(self[item])))
        lines.append(horizonal_rule)
        return '\n'.join(lines)

    def _mro(self):
        mro = [c.__name__ for c in self.__class__.__mro__]
        return mro

    def _isa(self, obj_type=None):
        if obj_type is None:
            prefix = 'base'
            if 'GenericNetwork' in self._mro():
                prefix = 'network'
            elif 'GenericGeometry' in self._mro():
                prefix = 'geometry'
            elif 'GenericPhase' in self._mro():
                prefix = 'phase'
            elif 'GenericPhysics' in self._mro():
                prefix = 'physics'
            elif 'GenericAlgorithm' in self._mro():
                prefix = 'algorithm'
            return prefix
        else:
            mro = [s.lower() for s in self._mro()]
            temp = [s.replace('generic', '') for s in mro
                    if s.startswith('generic')]
            mro.extend(temp)
            flag = False
            if obj_type.lower() in mro:
                flag = True
            return flag
