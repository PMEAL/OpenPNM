import warnings
import uuid
import numpy as np
from collections import namedtuple
from openpnm.utils import Workspace, logging
from openpnm.utils.misc import PrintableList, SettingsDict, Docorator
docstr = Docorator()
logger = logging.getLogger(__name__)
ws = Workspace()


@docstr.get_sections(base='Base', sections=['Parameters'])
class Base(dict):
    r"""
    Contains methods for working with the data in the OpenPNM dict objects

    Parameters
    ----------
    project : OpenPNM Project object, optional
        The Project with which the object should be assigned.  If not supplied
        then a new Project is created
    name : string, optional
        The unique name of the object.  If not given one will be generated.

    Np : int, default is 0
        The total number of pores to be assigned to the object
    Nt : int, default is 0
        The total number of throats to be assigned to the object

    Notes
    -----
    This Base class is used as the template for all other OpenPNM objects,
    including Networks, Geometries, Phases, Physics, and Algorithms.  This
    class is a subclass of the standard ``dict`` so has the usual methods such
    as ``pop`` and ``keys``, and has extra methods for working specifically
    with OpenPNM data.  These are outlined briefly in the following table:

    +----------------------+--------------------------------------------------+
    | Method or Attribute  | Functionality                                    |
    +======================+==================================================+
    | ``props``            | List of keys containing numerical arrays         |
    +----------------------+--------------------------------------------------+
    | ``labels``           | List of key containing boolean arrays            |
    +----------------------+--------------------------------------------------+
    | ``pores``            | Returns pore / throat indices that have given    |
    |                      | labels                                           |
    | ``throats``          |                                                  |
    +----------------------+--------------------------------------------------+
    | ``Ps``, ``Ts``       | Indices for ALL pores and throats on object      |
    +----------------------+--------------------------------------------------+
    | ``num_pores`` ,      | Counts the number of pores or throats with a     |
    |                      | given label                                      |
    | ``num_throats``      |                                                  |
    +----------------------+--------------------------------------------------+
    | ``Np``, ``Nt``       | Total number of pores and throats on the object  |
    +----------------------+--------------------------------------------------+
    | ``tomask``           | Converts a list of pore or throat indices to a   |
    |                      | boolean mask                                     |
    +----------------------+--------------------------------------------------+
    | ``toindices``        | Converts a boolean mask to pore or throat indices|
    +----------------------+--------------------------------------------------+
    | ``map_pores`` ,      | Given indices on object B returns corresponding  |
    |                      | indices on object A                              |
    | ``map_throats``      |                                                  |
    +----------------------+--------------------------------------------------+
    | ``interleave_data``  | Fetches data from associated objects into a      |
    |                      | single array                                     |
    +----------------------+--------------------------------------------------+
    | ``interpolate_data`` | Given pore or throat data, interpolate the other |
    +----------------------+--------------------------------------------------+
    | ``filter_by_label``  | Given indices find those with specific labels    |
    +----------------------+--------------------------------------------------+
    | ``show_hist``        | Method for quickly plotting histograms of data   |
    +----------------------+--------------------------------------------------+
    | ``check_data_health``| Ensures all data arrays are valid and complete   |
    +----------------------+--------------------------------------------------+


    In addition to the above methods, there are a few attributes which provide
    access to useful items:

    +----------------+--------------------------------------------------------+
    | Attribute      | Functionality                                          |
    +================+========================================================+
    | ``name``       | The string name of the object, unique to each Project  |
    +----------------+--------------------------------------------------------+
    | ``settings``   | A dictionary containing various setting values         |
    +----------------+--------------------------------------------------------+
    | ``project``    | A handle to the Project containing the object          |
    +----------------+--------------------------------------------------------+

    Examples
    --------
    It is possible to create an instance of Base, although it is not very
    useful except for demonstration purposes as done here.

    >>> import openpnm as op
    >>> obj = op.core.Base(Np=4, Nt=5)

    Now query the object for its basic properties:

    >>> obj.Np, obj.Nt  # Number of pores and throats
    (4, 5)

    Add a label to the object, as a boolean with True where the label applies:

    >>> obj['pore.new_label'] = [ True, False, False, True]

    See list of available labels and confirm new_label was added:

    >>> print(obj.labels())
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
    1     : pore.all
    2     : pore.new_label
    3     : throat.all
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――

    Use the label to fetch pores where it was applied:

    >>> Ps = obj.pores('new_label')
    >>> print(Ps)
    [0 3]

    Find the number of pores with label:

    >>> print(obj.num_pores('new_label'))
    2

    Convert between indices and boolean mask

    >>> mask = obj.tomask(throats=[0, 2, 4])
    >>> print(mask)
    [ True False  True False  True]
    >>> inds = obj.toindices(mask)
    >>> print(inds)
    [0 2 4]

    """

    def __new__(cls, *args, **kwargs):
        instance = super(Base, cls).__new__(cls, *args, **kwargs)
        # It is necessary to set the SettingsDict here since some classes
        # use it before calling super.__init__()
        instance.settings = SettingsDict()
        instance.settings['name'] = None
        instance.settings['_uuid'] = str(uuid.uuid4())
        return instance

    def __init__(self, Np=0, Nt=0, name=None, project=None, settings={}):
        self.settings.setdefault('prefix', 'base')
        self.settings.update(settings)
        super().__init__()
        if project is None:
            project = ws.new_project()
        if name is None:
            name = project._generate_name(self)
        project._validate_name(name)
        project.extend(self)
        self.settings['name'] = name
        self.update({'pore.all': np.ones(shape=(Np, ), dtype=bool)})
        self.update({'throat.all': np.ones(shape=(Nt, ), dtype=bool)})

    def __repr__(self):
        return '<%s object at %s>' % (self.__class__.__module__, hex(id(self)))

    def __eq__(self, other):
        return hex(id(self)) == hex(id(other))

    def __setitem__(self, key, value):
        r"""
        This is a subclass of the default __setitem__ behavior.  The main aim
        is to limit what type and shape of data can be written to protect
        the integrity of the network.  Specifically, this means only Np or Nt
        long arrays can be written, and they must be called 'pore.***' or
        'throat.***'.  Also, any scalars are cast into full length vectors.

        """
        if value is None:
            return
        # Check 1: If value is a dictionary, break it into constituent arrays
        # and recursively call __setitem__ on each
        if hasattr(value, 'keys'):
            for item in value.keys():
                prop = item.replace('pore.', '').replace('throat.', '')
                self.__setitem__(key+'.'+prop, value[item])
            return

        # Check 3: Enforce correct dict naming
        element = key.split('.')[0]
        if element not in ['pore', 'throat']:
            raise Exception('All keys must start with either pore or throat')

        # Check 2: If adding a new key, make sure it has no conflicts
        if self.project:
            proj = self.project
            boss = proj.find_full_domain(self)
            keys = boss.keys(mode='all', deep=True)
        else:
            boss = None
            keys = self.keys()
        # Prevent 'pore.foo.bar' when 'pore.foo' present
        long_keys = [i for i in keys if i.count('.') > 1]
        key_root = '.'.join(key.split('.')[:2])
        if (key.count('.') > 1) and (key_root in keys):
            raise Exception('Cannot create ' + key + ' when '
                            + key_root + ' is already defined')
        # Prevent 'pore.foo' when 'pore.foo.bar' is present
        if (key.count('.') == 1) and any([i.startswith(key) for i in long_keys]):
            hit = [i for i in keys if i.startswith(key)][0]
            raise Exception('Cannot create ' + key + ' when '
                            + hit + ' is already defined')
        # Prevent writing pore.foo on boss when present on subdomain
        if boss:
            if boss is self and (key not in ['pore.all', 'throat.all']):
                if (key in keys) and (key not in self.keys()):
                    raise Exception('Cannot create ' + key + ' when it is'
                                    + ' already defined on a subdomain')

        # This check allows subclassed numpy arrays through, eg. with units
        if not isinstance(value, np.ndarray):
            value = np.array(value, ndmin=1)  # Convert value to an ndarray

        # Skip checks for 'coords', 'conns'
        if key in ['pore.coords', 'throat.conns']:
            super(Base, self).__setitem__(key, value)
            return

        # Skip checks for protected props, and prevent changes if defined
        protected_keys = ['all']
        if key.split('.')[1] in protected_keys:
            if key in self.keys():
                if np.shape(self[key]) == (0, ):
                    super(Base, self).__setitem__(key, value)
                else:
                    warnings.warn(key+' is already defined.')
            else:
                super(Base, self).__setitem__(key, value)
            return

        # Write value to dictionary
        if np.shape(value)[0] == 1:  # If value is scalar
            value = np.ones((self._count(element), ), dtype=value.dtype)*value
            super(Base, self).__setitem__(key, value)
        elif np.shape(value)[0] == self._count(element):
            super(Base, self).__setitem__(key, value)
        else:
            if self._count(element) == 0:
                self.update({key: value})
            else:
                raise Exception('Provided array is wrong length for ' + key)

    def __getitem__(self, key):
        element, prop = key.split('.', 1)
        if key in self.keys():
            # Get values if present on self
            vals = super().__getitem__(key)
        elif key in self.keys(mode='all', deep=True):
            # Interleave values from geom if found there
            vals = self.interleave_data(key)
        elif any([k.startswith(key + '.') for k in self.keys()]):
            # Create a subdict of values present on self
            vals = {}
            keys = self.keys()
            vals.update({k: self.get(k) for k in keys if k.startswith(key + '.')})
        elif any([k.startswith(key + '.') for k in self.keys(mode='all',
                                                             deep=True)]):
            # Create a subdict of values in subdomains by interleaving
            vals = {}
            keys = self.keys(mode='all', deep=True)
            vals.update({k: self.interleave_data(k) for k in keys
                         if k.startswith(key + '.')})
        # Attempt to run model when missing data.
        elif hasattr(self, 'models') and key in self.models:
            self.regenerate_models(key)
            vals = super().__getitem__(key)
        # The following is probably a bad idea, but just trying it for fun
        # elif self.settings['interpolation_mode'] is not None:
        #     temp = list(set(['pore', 'throat']).difference(set([element])))[0]
        #     vals = self.interpolate_data(temp + '.' + prop,
        #                                  mode=self.settings['interpolation_mode'])
        #     self[element + '.' + prop] = vals
        else:
            raise KeyError(key)
        return vals

    def _set_name(self, name, validate=True):
        old_name = self.settings['name']
        if name == old_name:
            return
        if name is None:
            name = self.project._generate_name(self)
        if validate:
            self.project._validate_name(name)
        self.settings['name'] = name
        # Rename any label arrays in other objects
        for item in self.project:
            if 'pore.' + old_name in item.keys():
                item['pore.'+name] = item.pop('pore.' + old_name)
            if 'throat.' + old_name in item.keys():
                item['throat.' + name] = item.pop('throat.' + old_name)

    def _get_name(self):
        return self.settings['name']

    name = property(_get_name, _set_name)

    def _get_project(self):
        for proj in ws.values():
            if self in proj:
                return proj

    project = property(fget=_get_project)

    @property
    def network(self):
        r"""
        A shortcut to get a handle to the associated network
        There can only be one so this works
        """
        return self.project.network

    def clear(self, element=None, mode='all'):
        r"""
        A subclassed version of the standard dict's clear method.  This can be
        used to selectively clear certain data from the object, including
        properties and/or labels.  Importantly, it does NOT clear items that
        are required to maintain the integrity of the simulation.  These are
        arrays that define the topology (ie. 'pore.all', 'pore.coords',
        'throat.all', 'throat.conns'), as well as arrays that indicate
        associations bewteen objects (ie. 'pore.geo_01').

        Parameters
        ----------
        element : string or list of strings
            Can be either 'pore' or 'throat', which specifies whether 'pore'
            and/or 'throat' arrays should be cleared.  The default is both.

        mode : string or list of strings
            This controls what is cleared from the object.  Options are:

            **'props'** : Removes all numerical property values from the object
            dictionary

            **'model_data'** : Removes only numerical data that were produced
            by an associated model

            **'labels'** : Removes all labels from the object dictionary,
            except those relating to the pore and throat locations of
            associated objects

            **'all'** : Removes both 'props' and 'labels'

        Notes
        -----
        If you wish to selectively remove some properties but not others, use
        something like ``del object['pore.blah']`` at the Python prompt. This
        can also be done in a for-loop to remove a list of items.

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> len(pn.labels())  # There are 10 total labels on the network
        12
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
            for item in list(self.keys()):
                temp = '.'.join(item.split('.')[0:2])
                if temp in self.models.keys():
                    logger.info('deleting ' + item)
                    del self[item]
            mode.remove('model_data')
        for item in self.keys(mode=mode, element=element):
            if item not in protected:
                if item.split('.')[1] not in self.project.names:
                    del self[item]

    def keys(self, element=None, mode=None, deep=False):
        r"""
        This subclass works exactly like ``keys`` when no arguments are passed,
        but optionally accepts an ``element`` and a ``mode``, which filters
        the output to only the requested keys.

        The default behavior is exactly equivalent to the normal ``keys``
        method.

        Parameters
        ----------
        element : string
            Can be either 'pore' or 'throat', which limits the returned list of
            keys to only 'pore' or 'throat' keys.  If neither is given, then
            both are assumed.

        mode : string (optional)
            Controls which keys are returned.  Options are:

            **'labels'** : Limits the returned list of keys to only 'labels'
            (boolean arrays)

            **'props'** : Limits he return list of keys to only 'props'
            (numerical arrays).

            **'all'** : Returns both 'labels' and 'props'.  This is equivalent
            to sending a list of both 'labels' and 'props'.

            If no mode is specified then the normal KeysView object is
            returned.

        deep : Boolean
            If set to ``True`` then the keys on all associated subdomain
            objects are returned as well.

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

        if deep:
            if self._isa('phase'):
                for item in self.project.find_physics(phase=self):
                    temp += item.keys(element=element, mode=mode, deep=False)
            if self._isa('network'):
                for item in self.project.geometries().values():
                    temp += item.keys(element=element, mode=mode, deep=False)

        return temp

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

            **'all'** : Returns all properties on the object (default)

            **'models'** : Returns only properties that are associated with a
            model

            **'constants'** : returns data values that were *not* generated by
            a model, but manaully created.

        deep : Boolean
            If set to ``True`` then the props on all associated subdomain
            objects are returned as well.

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
        # Repeat for associated objects if deep is True
        if deep:
            if self._isa('phase'):
                for item in self.project.find_physics(phase=self):
                    vals += item.props(element=element, mode=mode, deep=False)
            if self._isa('network'):
                for item in self.project.geometries().values():
                    vals += item.props(element=element, mode=mode, deep=False)
        return vals

    def _get_labels(self, element, locations, mode):
        r"""
        This is the actual label getter method, but it should not be called
        directly.  Use ``labels`` instead.
        """
        # Parse inputs
        locations = self._parse_indices(locations)
        element = self._parse_element(element=element)
        # Collect list of all pore OR throat labels
        labels = self.keys(mode='labels', element=element)
        labels.sort()
        labels = np.array(labels)  # Convert to ND-array for following checks
        # Make an 2D array with locations in rows and labels in cols
        arr = np.vstack([self[item][locations] for item in labels]).T
        num_hits = np.sum(arr, axis=0)  # Number of locations with each label
        if mode in ['or', 'union', 'any']:
            temp = labels[num_hits > 0]
        elif mode in ['and', 'intersection']:
            temp = labels[num_hits == locations.size]
        elif mode in ['xor', 'exclusive_or']:
            temp = labels[num_hits == 1]
        elif mode in ['nor', 'not', 'none']:
            temp = labels[num_hits == 0]
        elif mode in ['nand']:
            temp = labels[num_hits == (locations.size - 1)]
        elif mode in ['xnor', 'nxor']:
            temp = labels[num_hits > 1]
        else:
            raise Exception('Unrecognized mode:'+str(mode))
        return PrintableList(temp)

    def labels(self, pores=[], throats=[], element=None, mode='union'):
        r"""
        Returns a list of labels present on the object

        Additionally, this function can return labels applied to a specified
        set of pores or throats

        Parameters
        ----------
        element : string
            Controls whether pore or throat labels are returned.  If empty then
            both are returned (default).

        pores (or throats) : array_like
            The pores (or throats) whose labels are sought.  If left empty a
            list containing all pore and throat labels is returned.

        mode : string, optional
            Controls how the query should be performed.  Only applicable
            when ``pores`` or ``throats`` are specified:

            **'or', 'union', 'any'**: (default) Returns the labels that are
            assigned to *any* of the given locations.

            **'and', 'intersection', 'all'**: Labels that are present on *all*
            the given locations.

            **'xor', 'exclusive_or'** : Labels that are present on *only one*
            of the given locations.

            **'nor', 'none', 'not'**: Labels that are *not* present on any of
            the given locations.

            **'nand'**: Labels that are present on *all but one* of the given
            locations

            **'xnor'**: Labels that are present on *more than one* of the given
            locations.  'nxor' is also accepted.

        Returns
        -------
        A list containing the labels on the object.  If ``pores`` or
        ``throats`` are given, the results are filtered according to the
        specified ``mode``.

        See Also
        --------
        props
        keys

        Notes
        -----
        Technically, *'nand'* and *'xnor'* should also return pores with *none*
        of the labels but these are not included.  This makes the returned list
        more useful.

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> pn.labels(pores=[11, 12])
        ['pore.all', 'pore.internal', 'pore.left', 'pore.surface']
        """
        # Short-circuit query when no pores or throats are given
        if (np.size(pores) == 0) and (np.size(throats) == 0):
            labels = PrintableList(self.keys(element=element, mode='labels'))
        elif (np.size(pores) > 0) and (np.size(throats) > 0):
            raise Exception('Cannot perform label query on pores and '
                            + 'throats simultaneously')
        elif np.size(pores) > 0:
            labels = self._get_labels(element='pore', locations=pores,
                                      mode=mode)
        elif np.size(throats) > 0:
            labels = self._get_labels(element='throat', locations=throats,
                                      mode=mode)
        return labels

    def set_label(self, label, pores=None, throats=None, mode='add'):
        r"""
        Creates or updates a label array

        Parameters
        ----------
        label : string
                The label to apply to the specified locations
        pores : array_like
            A list of pore indices or a boolean mask of where given label
            should be added or removed (see ``mode``)
        throats : array_like
            A list of throat indices or a boolean mask of where given label
            should be added or removed (see ``mode``)
        mode : string
            Controls how the labels are handled.  Options are:

            *'add'* - Adds the given label to the specified locations while
            keeping existing labels (default)

            *'overwrite'* - Removes existing label from all locations before
            adding the label in the specified locations

            *'remove'* - Removes the  given label from the specified locations
            leaving the remainder intact.

            *'purge'* - Removes the specified label from the object

        """
        if mode == 'purge':
            if label.split('.')[0] in ['pore', 'throat']:
                if label in self.labels():
                    del self[label]
                else:
                    logger.warning(label + ' is not a label, skpping')
            else:
                self.set_label(label='pore.'+label, mode='purge')
                self.set_label(label='throat.'+label, mode='purge')
        else:
            if label.split('.')[0] in ['pore', 'throat']:
                label = label.split('.', 1)[1]
            if pores is not None:
                pores = self._parse_indices(pores)
                if (mode == 'overwrite') or ('pore.'+label not in self.labels()):
                    self['pore.' + label] = False
                if mode in ['remove']:
                    self['pore.' + label][pores] = False
                else:
                    self['pore.' + label][pores] = True
            if throats is not None:
                throats = self._parse_indices(throats)
                if (mode == 'overwrite') or ('throat.'+label not in self.labels()):
                    self['throat.' + label] = False
                if mode in ['remove']:
                    self['throat.' + label][throats] = False
                else:
                    self['throat.' + label][throats] = True
            if pores is None and throats is None:
                del self

    def _get_indices(self, element, labels='all', mode='or'):
        r"""
        This is the actual method for getting indices, but should not be called
        directly.  Use ``pores`` or ``throats`` instead.
        """
        # Parse and validate all input values.
        element = self._parse_element(element, single=True)
        labels = self._parse_labels(labels=labels, element=element)
        if element+'.all' not in self.keys():
            raise Exception('Cannot proceed without {}.all'.format(element))

        # Begin computing label array
        if mode in ['or', 'any', 'union']:
            union = np.zeros_like(self[element+'.all'], dtype=bool)
            for item in labels:  # Iterate over labels and collect all indices
                union = union + self[element+'.'+item.split('.')[-1]]
            ind = union
        elif mode in ['and', 'all', 'intersection']:
            intersect = np.ones_like(self[element+'.all'], dtype=bool)
            for item in labels:  # Iterate over labels and collect all indices
                intersect = intersect*self[element+'.'+item.split('.')[-1]]
            ind = intersect
        elif mode in ['xor', 'exclusive_or']:
            xor = np.zeros_like(self[element+'.all'], dtype=int)
            for item in labels:  # Iterate over labels and collect all indices
                info = self[element+'.'+item.split('.')[-1]]
                xor = xor + np.int8(info)
            ind = (xor == 1)
        elif mode in ['nor', 'not', 'none']:
            nor = np.zeros_like(self[element+'.all'], dtype=int)
            for item in labels:  # Iterate over labels and collect all indices
                info = self[element+'.'+item.split('.')[-1]]
                nor = nor + np.int8(info)
            ind = (nor == 0)
        elif mode in ['nand']:
            nand = np.zeros_like(self[element+'.all'], dtype=int)
            for item in labels:  # Iterate over labels and collect all indices
                info = self[element+'.'+item.split('.')[-1]]
                nand = nand + np.int8(info)
            ind = (nand < len(labels)) * (nand > 0)
        elif mode in ['xnor', 'nxor']:
            xnor = np.zeros_like(self[element+'.all'], dtype=int)
            for item in labels:  # Iterate over labels and collect all indices
                info = self[element+'.'+item.split('.')[-1]]
                xnor = xnor + np.int8(info)
            ind = (xnor > 1)
        else:
            raise Exception('Unsupported mode: '+mode)
        # Extract indices from boolean mask
        ind = np.where(ind)[0]
        ind = ind.astype(dtype=int)
        return ind

    def pores(self, labels='all', mode='or', asmask=False, target=None):
        r"""
        Returns pore indicies where given labels exist, according to the logic
        specified by the ``mode`` argument.

        Parameters
        ----------
        labels : string or list of strings
            The label(s) whose pores locations are requested.  This argument
            also accepts '*' for wildcard searches.

        mode : string
            Specifies how the query should be performed.  The options are:

            **'or', 'union', 'any'** : (default) Pores with *one or more* of
            the given labels are returned.

            **'and', 'intersection', 'all'** : Pores with *all* of the given
            labels are returned.

            **'xor', 'exclusive_or'** : Pores with *only one* of the given
            labels are returned.

            **'nor', 'none', 'not'** : Pores with *none* of the given labels
            are returned.

            **'nand'** : Pores with *not all* of the given labels are
            returned.

            **'xnor'** : Pores with *more than one* of the given labels are
            returned.

        asmask : boolean
            If ``True`` then a boolean array of length Np is returned with
            ``True`` values indicating the pores that satisfy the query.

        target : OpenPNM Base object
            If given, the returned indices will be indexed relative to the
            ``target`` object.  This can be used to determine how indices on
            one object map onto another object.

        Returns
        -------
        A Numpy array containing pore indices filtered by the logic specified
        in ``mode``.

        See Also
        --------
        throats
        map_pores

        Notes
        -----
        Technically, *nand* and *xnor* should also return pores with *none* of
        the labels but these are not included.  This makes the returned list
        more useful.

        To perform more complex or compound queries, you can opt to receive
        the result a a boolean mask (``asmask=True``), then manipulate the
        arrays manually.

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> Ps = pn.pores(labels=['top', 'back'], mode='union')
        >>> Ps[:5]  # Look at first 5 pore indices
        array([ 4,  9, 14, 19, 20])
        >>> pn.pores(labels=['top', 'back'], mode='xnor')
        array([ 24,  49,  74,  99, 124])
        """
        ind = self._get_indices(element='pore', labels=labels, mode=mode)
        if target is not None:
            ind = target.map_pores(pores=ind, origin=self, filtered=True)
        if asmask:
            if target is not None:
                ind = target.tomask(pores=ind)
            else:
                ind = self.tomask(pores=ind)
        return ind

    @property
    def Ps(self):
        r"""
        A shortcut to get a list of all pores on the object
        """
        return np.arange(0, self.Np)

    def throats(self, labels='all', mode='or', asmask=False, target=None):
        r"""
        Returns throat locations where given labels exist, according to the
        logic specified by the ``mode`` argument.

        Parameters
        ----------
        labels : string or list of strings
            The throat label(s) whose locations are requested.  If omitted,
            'all' throat inidices are returned.  This argument also accepts
            '*' for wildcard searches.

        mode : string
            Specifies how the query should be performed.  The options are:

            **'or', 'union', 'any'** : (default) Throats with *one or more* of
            the given labels are returned.

            **'and', 'intersection', 'all'** : Throats with *all* of the given
            labels are returned.

            **'xor', 'exclusive_or'** : Throats with *only one* of the given
            labels are returned.

            **'nor', 'none', 'not'** : Throats with *none* of the given labels
            are returned.

            **'nand'** : Throats with *not all* of the given labels are
            returned.

            **'xnor'** : Throats with *more than one* of the given labels are
            returned.

        asmask : boolean
            If ``True`` then a boolean array of length Nt is returned with
            ``True`` values indicating the throats that satisfy the query.

        target : OpenPNM Base object
            If given, the returned indices will be indexed relative to the
            ``target`` object.  This can be used to determine how indices on
            one object map onto another object.

        Returns
        -------
        A Numpy array containing throat indices filtered by the logic specified
        in ``mode``.

        See Also
        --------
        pores
        map_throats

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[3, 3, 3])
        >>> Ts = pn.throats()
        >>> Ts[0:5]  # Look at first 5 throat indices
        array([0, 1, 2, 3, 4])

        """
        ind = self._get_indices(element='throat', labels=labels, mode=mode)
        if target is not None:
            ind = target.map_throats(throats=ind, origin=self, filtered=True)
        if asmask:
            if target is not None:
                ind = target.tomask(throats=ind)
            else:
                ind = self.tomask(throats=ind)
        return ind

    @property
    def Ts(self):
        r"""
        A shortcut to get a list of all throats on the object
        """
        return np.arange(0, self.Nt)

    def _map(self, ids, element, filtered):
        ids = np.array(ids, dtype=np.int64)
        locations = self._get_indices(element=element)
        self_in_ids = np.isin(ids, self[element+'._id'], assume_unique=True)
        ids_in_self = np.isin(self[element+'._id'], ids, assume_unique=True)
        mask = np.zeros(shape=ids.shape, dtype=bool)
        mask[self_in_ids] = True
        ind = np.ones_like(mask, dtype=np.int64) * -1
        ind[self_in_ids] = locations[ids_in_self]
        if filtered:
            return ind[mask]
        t = namedtuple('index_map', ('indices', 'mask'))
        return t(ind, mask)

    def map_pores(self, pores, origin, filtered=True):
        r"""
        Given a list of pore on a target object, finds indices of those pores
        on the calling object

        Parameters
        ----------
        pores : array_like
            The indices of the pores on the object specifiedin ``origin``

        origin : OpenPNM Base object
            The object corresponding to the indices given in ``pores``

        filtered : boolean (default is ``True``)
            If ``True`` then a ND-array of indices is returned with missing
            indices removed, otherwise a named-tuple containing both the
            ``indices`` and a boolean ``mask`` with ``False`` indicating
            which locations were not found.

        Returns
        -------
        Pore indices on the calling object corresponding to the same pores
        on the ``origin`` object.  Can be an array or a tuple containing an
        array and a mask, depending on the value of ``filtered``.

        See Also
        --------
        pores
        map_throats

        """
        ids = origin['pore._id'][pores]
        return self._map(element='pore', ids=ids, filtered=filtered)

    def map_throats(self, throats, origin, filtered=True):
        r"""
        Given a list of throats on a target object, finds indices of
        those throats on the calling object

        Parameters
        ----------
        throats : array_like
            The indices of the throats on the object specified in ``origin``

        origin : OpenPNM Base object
            The object corresponding to the indices given in ``throats``

        filtered : boolean (default is ``True``)
            If ``True`` then a ND-array of indices is returned with missing
            indices removed, otherwise a named-tuple containing both the
            ``indices`` and a boolean ``mask`` with ``False`` indicating
            which locations were not found.

        Returns
        -------
        Throat indices on the calling object corresponding to the same throats
        on the ``origin`` object.  Can be an array or a tuple containing an
        array and a mask, depending on the value of ``filtered``.

        See Also
        --------
        throats
        map_pores

        """
        ids = origin['throat._id'][throats]
        return self._map(element='throat', ids=ids, filtered=filtered)

    def _tomask(self, indices, element):
        r"""
        This is a generalized version of tomask that accepts a string of
        'pore' or 'throat' for programmatic access.
        """
        element = self._parse_element(element, single=True)
        indices = self._parse_indices(indices)
        N = np.shape(self[element + '.all'])[0]
        ind = np.array(indices, ndmin=1)
        mask = np.zeros((N, ), dtype=bool)
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
        just a convenience function and is a complement to ``tomask``.

        """
        if np.amax(mask) > 1:
            raise Exception('Received mask does not appear to be boolean')
        mask = np.array(mask, dtype=bool)
        indices = self._parse_indices(mask)
        return indices

    def interleave_data(self, prop):
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
        # Fetch sources list depending on type of self
        proj = self.project
        if self._isa() in ['network', 'geometry']:
            sources = list(proj.geometries().values())
        elif self._isa() in ['phase', 'physics']:
            sources = list(proj.find_physics(phase=self))
        elif self._isa() in ['algorithm', 'base']:
            sources = [self]
        else:
            raise Exception('Unrecognized object type, cannot find dependents')

        # Get generalized element and array length
        element = self._parse_element(prop.split('.')[0], single=True)
        N = self.project.network._count(element)

        # Attempt to fetch the requested array from each object
        arrs = [obj.get(prop, None) for obj in sources]

        # Check for missing sources, and add None to arrs if necessary
        if N > sum([obj._count(element) for obj in sources]):
            arrs.append(None)

        # Obtain list of locations for inserting values
        locs = [self._get_indices(element, item.name) for item in sources]

        if np.all([item is None for item in arrs]):  # prop not found anywhere
            raise KeyError(prop)

        # Let's start by handling the easy cases first
        if not any([a is None for a in arrs]):
            # All objs present and array found on all objs
            shape = list(arrs[0].shape)
            shape[0] = N
            types = [a.dtype for a in arrs]
            if len(set(types)) == 1:
                # All types are the same
                temp_arr = np.ones(shape, dtype=types[0])
                for vals, inds in zip(arrs, locs):
                    temp_arr[inds] = vals
                return temp_arr  # Return early because it's just easier
            if all([a.dtype in [float, int, bool] for a in arrs]):
                # All types are numeric, make float
                temp_arr = np.ones(shape, dtype=float)
                for vals, inds in zip(arrs, locs):
                    temp_arr[inds] = vals
                return temp_arr  # Return early because it's just easier

        # Now handle the complicated cases
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
        dummy_val = {'numeric': np.nan, 'boolean': False, 'other': None}

        # Create an empty array of the right type and shape
        for item in arrs:
            if item is not None:
                if len(item.shape) == 1:
                    temp_arr = np.zeros((N, ), dtype=item.dtype)
                else:
                    temp_arr = np.zeros((N, item.shape[1]), dtype=item.dtype)
                temp_arr.fill(dummy_val[atype[0]])

        sizes = [np.size(a) for a in arrs]
        # Convert int arrays to float IF NaNs are expected
        if temp_arr.dtype.name.startswith('int') and \
           (np.any([i is None for i in arrs]) or np.sum(sizes) != N):
            temp_arr = temp_arr.astype(float)
            temp_arr.fill(np.nan)

        # Fill new array with values in the corresponding locations
        for vals, inds in zip(arrs, locs):
            if vals is not None:
                temp_arr[inds] = vals
            else:
                temp_arr[inds] = dummy_val[atype[0]]

        # Check if any arrays have units, if so then apply them to result
        # Importing unyt significantly adds to our import time, we also
        # currently don't use this package extensively, so we're not going
        # to support it for now.

        # if any([hasattr(a, 'units') for a in arrs]):
        #     [a.convert_to_mks() for a in arrs if hasattr(a, 'units')]
        #     units = [a.units.__str__() for a in arrs if hasattr(a, 'units')]
        #     if len(units) > 0:
        #         if len(set(units)) == 1:
        #             temp_arr *= np.array([1]) * getattr(unyt, units[0])
        #         else:
        #             raise Exception('Units on the interleaved array are not equal')

        return temp_arr

    def interpolate_data(self, propname, mode='mean'):
        r"""
        Determines a pore (or throat) property as the average of it's
        neighboring throats (or pores)

        Parameters
        ----------
        propname: string
            The dictionary key to the values to be interpolated.
        mode : string
            The method used for interpolation.  Options are 'mean' (default),
            'min', and 'max'.

        Returns
        -------
        vals : ND-array
            An array containing interpolated pore (or throat) data

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[3, 1, 1])
        >>> pn['pore.value'] = [1, 2, 3]
        >>> pn.interpolate_data('pore.value')
        array([1.5, 2.5])

        """
        from openpnm.models.misc import from_neighbor_throats, from_neighbor_pores
        if propname.startswith('throat'):
            values = from_neighbor_throats(target=self, prop=propname, mode=mode)
        elif propname.startswith('pore'):
            values = from_neighbor_pores(target=self, prop=propname, mode=mode)
        if hasattr(self[propname], 'units'):
            values *= self[propname].units
        return values

    def filter_by_label(self, pores=[], throats=[], labels=None, mode='or'):
        r"""
        Returns which of the supplied pores (or throats) has the specified
        label(s)

        Parameters
        ----------
        pores, or throats : array_like
            List of pores or throats to be filtered

        labels : list of strings
            The labels to apply as a filter

        mode : string

            Controls how the filter is applied.  Options include:

            **'or', 'union', 'any'**: (default) Returns a list of the given
            locations where *any* of the given labels exist.

            **'and', 'intersection', 'all'**: Only locations where *all* the
            given labels are found.

            **'xor', 'exclusive_or'**: Only locations where exactly *one* of
            the given labels are found.

            **'nor', 'none', 'not'**: Only locations where *none* of the given
            labels are found.

            **'nand'** : Only locations with *some but not all* of the given
            labels are returned.

            **'xnor'** : Only locations with *more than one* of the given
            labels are returned.

        Returns
        -------
        A list of pores (or throats) that have been filtered according the
        given criteria.  The returned list is a subset of the received list of
        pores (or throats).

        See Also
        --------
        pores
        throats

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> pn.filter_by_label(pores=[0, 1, 25, 32], labels='left')
        array([0, 1])
        >>> Ps = pn.pores(['top', 'bottom', 'back'], mode='or')
        >>> pn.filter_by_label(pores=Ps, labels=['top', 'back'],
        ...                    mode='and')
        array([ 24,  49,  74,  99, 124])
        """
        # Convert inputs to locations and element
        if (np.size(throats) > 0) and (np.size(pores) > 0):
            raise Exception('Can only filter either pores OR labels')
        if np.size(pores) > 0:
            element = 'pore'
            locations = self._parse_indices(pores)
        elif np.size(throats) > 0:
            element = 'throat'
            locations = self._parse_indices(throats)
        else:
            return(np.array([], dtype=int))
        labels = self._parse_labels(labels=labels, element=element)
        labels = [element+'.'+item.split('.')[-1] for item in labels]
        all_locs = self._get_indices(element=element, labels=labels, mode=mode)
        mask = self._tomask(indices=all_locs, element=element)
        ind = mask[locations]
        return locations[ind]

    def num_pores(self, labels='all', mode='or'):
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

            **'or', 'union', 'any'** : (default) Pores with *one or more* of
            the given labels are counted.

            **'and', 'intersection', 'all'** : Pores with *all* of the given
            labels are counted.

            **'xor', 'exclusive_or'** : Pores with *only one* of the given
            labels are counted.

            **'nor', 'none', 'not'** : Pores with *none* of the given labels
            are counted.

            **'nand'** : Pores with *some but not all* of the given labels are
            counted.

            **'xnor'** : Pores with *more than one* of the given labels are
            counted.

        Returns
        -------
        Np : int
            Number of pores with the specified labels

        See Also
        --------
        num_throats
        count

        Notes
        -----
        Technically, *'nand'* and *'xnor'* should also count pores with *none*
        of the labels, however, to make the count more useful these are not
        included.

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> pn.num_pores()
        125
        >>> pn.num_pores(labels=['top'])
        25
        >>> pn.num_pores(labels=['top', 'front'], mode='or')
        45
        >>> pn.num_pores(labels=['top', 'front'], mode='xnor')
        5

        """
        # Count number of pores of specified type
        Ps = self._get_indices(labels=labels, mode=mode, element='pore')
        Np = np.shape(Ps)[0]
        return Np

    @property
    def Np(self):
        r"""
        A shortcut to query the total number of pores on the object'
        """
        return np.shape(self.get('pore.all'))[0]

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

            **'or', 'union', 'any'** : (default) Throats with *one or more* of
            the given labels are counted.

            **'and', 'intersection', 'all'** : Throats with *all* of the given
            labels are counted.

            **'xor', 'exclusive_or'** : Throats with *only one* of the given
            labels are counted.

            **'nor', 'none', 'not'** : Throats with *none* of the given labels
            are counted.

            **'nand'** : Throats with *some but not all* of the given labels
            are counted.

            **'xnor'** : Throats with *more than one* of the given labels are
            counted.

        Returns
        -------
        Nt : int
            Number of throats with the specified labels

        See Also
        --------
        num_pores
        count

        Notes
        -----
        Technically, *'nand'* and *'xnor'* should also count throats with
        *none* of the labels, however, to make the count more useful these are
        not included.

        """
        # Count number of pores of specified type
        Ts = self._get_indices(labels=labels, mode=mode, element='throat')
        Nt = np.shape(Ts)[0]
        return Nt

    @property
    def Nt(self):
        r"""
        A shortcut to query the total number of throats on the object'
        """
        return np.shape(self.get('throat.all'))[0]

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
        temp = np.size(super(Base, self).__getitem__(element+'.all'))
        return temp

    def show_hist(self,
                  props=['pore.diameter', 'throat.diameter', 'throat.length'],
                  bins=20, fontsize=14, **kwargs):
        r"""
        Show a quick plot of key property distributions.

        Parameters
        ----------
        props : string or list of strings
            The pore and/or throat properties to be plotted as histograms.  By
            default this function will show 'pore.diameter', 'throat.diameter',
            and 'throat.length'.

        bins : int or array_like
            The number of bins to use when generating the histogram.  If an
            array is given they are used as the bin spacing instead.

        fontsize : int
            Sets the font size temporarily.  The default size of matplotlib is
            10, which is too small for many screens.  This function has a
            default of 22, which does not overwrite the matplotlib setting.
            Note that you can override matplotlib setting globally with
            ``matplotlib.rcParams['font.size'] = 22``.

        Notes
        -----
        Other keyword arguments are passed to the ``matplotlib.pyplot.hist``
        function.
        """
        import matplotlib.pyplot as plt

        temp = plt.rcParams['font.size']
        plt.rcParams['font.size'] = fontsize
        if isinstance(props, str):
            props = [props]
        N = len(props)
        color = plt.cm.tab10(range(10))
        if N <= 3:
            r, c = 1, N
        elif N == 4:
            r, c = 2, 2
        else:
            r, c = N // 3 + 1, 3
        fig, ax = plt.subplots(r, c, figsize=(3*c, 3*r))
        axs = np.array(ax).flatten()
        i = None
        for i, _ in enumerate(props):
            try:
                # Update kwargs with some default values
                if 'edgecolor' not in kwargs.keys():
                    kwargs.update({'edgecolor': 'k'})
                if 'facecolor' not in kwargs:
                    kwargs.update({'facecolor': color[np.mod(i, 10)]})
                axs[i].hist(self[props[i]], bins=bins, **kwargs)
                axs[i].set_xlabel(props[i])
            except KeyError:
                pass
        # Hide unpopulated subplots from the grid
        for j in range(i + 1, len(axs)):
            axs[j].set_axis_off()
        plt.rcParams['font.size'] = temp
        plt.tight_layout(h_pad=0.9, w_pad=0.9)

    def check_data_health(self):
        r"""
        Check the health of pore and throat data arrays.

        Returns
        -------
        health: HealthDict object
            A  basic dictionary with an added ``health`` attribute that is
            ``True`` if all entries in the dict are deemed healthy
            (empty lists), or ``False`` otherwise.

        Examples
        --------
        >>> import openpnm
        >>> pn = openpnm.network.Cubic(shape=[5, 5, 5])
        >>> h = pn.check_data_health()
        >>> h.health
        True
        """
        health = self.project.check_data_health(obj=self)
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
            indices = np.array([], ndmin=1, dtype=int)
        locs = np.array(indices, ndmin=1)
        # If boolean array, convert to indices
        if locs.dtype == bool:
            if np.size(locs) == self.Np:
                locs = self.Ps[locs]
            elif np.size(locs) == self.Nt:
                locs = self.Ts[locs]
            else:
                raise Exception('Mask of locations must be either '
                                + 'Np nor Nt long')
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
        When ``single`` is ``False`` (default) a list containing the element(s)
        is returned.  When ``single`` is ``True`` a bare string containing the
        element is returned.
        """
        if element is None:
            element = ['pore', 'throat']
        # Convert element to a list for subsequent processing
        if isinstance(element, str):
            element = [element]
        # Convert 'pore.prop' and 'throat.prop' into just 'pore' and 'throat'
        element = [item.split('.')[0] for item in element]
        # Make sure all are lowercase
        element = [item.lower() for item in element]
        # Deal with an plurals
        element = [item.rsplit('s', maxsplit=1)[0] for item in element]
        for item in element:
            if item not in ['pore', 'throat']:
                raise Exception('All keys must start with either pore or throat')
        # Remove duplicates if any
        _ = [element.remove(L) for L in element if element.count(L) > 1]
        if single:
            if len(element) > 1:
                raise Exception('Both elements recieved when single element '
                                + 'allowed')
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
        if isinstance(labels, str):
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
            _ = [parsed_labels.remove(L) for L in parsed_labels
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
        if isinstance(mode, str):
            mode = [mode]
        for item in mode:
            if (allowed is not None) and (item not in allowed):
                raise Exception('\'mode\' must be one of the following: '
                                + allowed.__str__())
        # Remove duplicates, if any
        _ = [mode.remove(L) for L in mode if mode.count(L) > 1]
        if single:
            if len(mode) > 1:
                raise Exception('Multiple modes received when only one mode '
                                + 'is allowed by this method')
            mode = mode[0]
        return mode

    def _parse_prop(self, propname, element):
        r"""
        """
        element = self._parse_element(element, single=True)
        return element + '.' + propname.split('.')[-1]

    def __str__(self):
        horizontal_rule = '―' * 78
        lines = [horizontal_rule]
        lines.append(self.__module__.replace('__', '') + ' : ' + self.name)
        lines.append(horizontal_rule)
        lines.append("{0:<5s} {1:<45s} {2:<10s}".format('#',
                                                        'Properties',
                                                        'Valid Values'))
        fmt = "{0:<5d} {1:<45s} {2:>5d} / {3:<5d}"
        lines.append(horizontal_rule)
        props = self.props()
        props.sort()
        for i, item in enumerate(props):
            prop = item
            required = self._count(item.split('.')[0])
            if len(prop) > 35:  # Trim overly long prop names
                prop = prop[0:32] + '...'
            if self[item].dtype == object:  # Print objects differently
                invalid = [i for i in self[item] if i is None]
                defined = np.size(self[item]) - len(invalid)
                lines.append(fmt.format(i + 1, prop, defined, required))
            elif '._' not in prop:
                a = np.isnan(self[item])
                defined = np.shape(self[item])[0] \
                    - a.sum(axis=0, keepdims=(a.ndim-1) == 0)[0]
                lines.append(fmt.format(i + 1, prop, defined, required))
        lines.append(horizontal_rule)
        lines.append("{0:<5s} {1:<45s} {2:<10s}".format('#',
                                                        'Labels',
                                                        'Assigned Locations'))
        lines.append(horizontal_rule)
        labels = self.labels()
        labels.sort()
        fmt = "{0:<5d} {1:<45s} {2:<10d}"
        for i, item in enumerate(labels):
            prop = item
            if len(prop) > 35:
                prop = prop[0:32] + '...'
            if '._' not in prop:
                lines.append(fmt.format(i + 1, prop, np.sum(self[item])))
        lines.append(horizontal_rule)
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
        mro = [s.lower() for s in self._mro()]
        temp = [s.replace('generic', '') for s in mro
                if s.startswith('generic')]
        mro.extend(temp)
        flag = False
        if obj_type.lower() in mro:
            flag = True
        return flag
