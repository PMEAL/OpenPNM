import warnings
import uuid
from copy import deepcopy
import numpy as np
from openpnm.utils import Workspace, logging
from openpnm.utils import SettingsAttr
from openpnm.utils.misc import PrintableList, Docorator
docstr = Docorator()
logger = logging.getLogger(__name__)
ws = Workspace()

__all__ = ['Base']


@docstr.get_sections(base='BaseSettings', sections=['Parameters'])
class BaseSettings:
    r"""
    The default settings to use on instance of Base

    Parameters
    ----------
    prefix : str
        The default prefix to use when generating a name
    name : str
        The name of the object, which will be generated if not given
    uuid : str
        A universally unique identifier for the object to keep things straight

    """
    prefix = 'base'
    name = ''
    uuid = ''


@docstr.get_sections(base='Base', sections=['Parameters'])
class Base(dict):
    r"""
    Contains methods for working with the data in the OpenPNM dict objects

    Parameters
    ----------
    network : OpenPNM network object
        The network to which this object is associated
    settings : dataclass-like or dict, optional
        User defined settings for the object to override defaults. Can be a
        dataclass-type object with settings stored as attributes or a python
        dicionary of key-value pairs. Settings are stored in the ``settings``
        attribute of the object.
    name : string, optional
        A unique name to assign to the object for easier identification.  If
        not given one will be generated.

    Np : int, default is 0
        The total number of pores to be assigned to the object
    Nt : int, default is 0
        The total number of throats to be assigned to the object

    Notes
    -----
    This ``Base`` class is used as the template for all other OpenPNM objects,
    including Networks, Geometries, Phase, Physics, and Algorithms. This
    class is a subclass of the standard ``dict`` so has the usual methods such
    as ``pop`` and ``keys``, and has extra methods for working specifically
    with OpenPNM data.

    """

    def __new__(cls, *args, **kwargs):
        instance = super(Base, cls).__new__(cls, *args, **kwargs)
        instance._settings = None
        instance._settings_docs = None
        return instance

    def __init__(self, Np=0, Nt=0, network=None, name=None, project=None,
                 settings=None):
        super().__init__()
        self.settings = SettingsAttr(BaseSettings, settings)

        if project is None:
            if network is None:
                project = ws.new_project()
            else:
                project = network.project
        if name is None:
            name = project._generate_name(self)
        project._validate_name(name)
        project.extend(self)
        self.settings['name'] = name
        self.settings.uuid = str(uuid.uuid4())
        self.update({'pore.all': np.ones(shape=(Np, ), dtype=bool)})
        self.update({'throat.all': np.ones(shape=(Nt, ), dtype=bool)})

    def __repr__(self):
        module = self.__module__
        module = ".".join([x for x in module.split(".") if not x.startswith("_")])
        cname = self.__class__.__name__
        return f'<{module}.{cname} at {hex(id(self))}>'

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
            raise Exception('All keys must start with either pore, or throat')

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
        else:
            raise KeyError(key)
        return vals

    def __delitem__(self, key):
        try:
            super().__delitem__(key)
        except KeyError as e:
            d = self[key]  # if key is a nested dict, get all values
            for item in d.keys():
                super().__delitem__(item)

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
        """String representing the name of the object"""
        try:
            return self.settings['name']
        except AttributeError:
            return None

    name = property(_get_name, _set_name)

    def _get_project(self):
        """A shortcut to get a handle to the associated project."""
        for proj in ws.values():
            if self in proj:
                return proj

    project = property(fget=_get_project)

    def _set_settings(self, settings):
        self._settings = deepcopy(settings)
        if (self._settings_docs is None) and (settings.__doc__ is not None):
            self._settings_docs = settings.__doc__

    def _get_settings(self):
        """Dictionary containing object settings."""
        if self._settings is None:
            self._settings = SettingsAttr()
        if self._settings_docs is not None:
            self._settings.__dict__['__doc__'] = self._settings_docs
        return self._settings

    def _del_settings(self):
        self._settings = None

    settings = property(fget=_get_settings, fset=_set_settings, fdel=_del_settings)

    @property
    def _domain(self):
        try:
            return self.phase
        except AttributeError:
            return self.network

    @property
    def network(self):
        r"""
        A shortcut to get a handle to the associated network.
        There can only be one so this works.
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
        element : str or List[str]
            Can be either 'pore' or 'throat', which specifies whether 'pore'
            and/or 'throat' arrays should be cleared.  The default is both.
        mode : str or List[str]
            This controls what is cleared from the object.  Options are:
                **'props'** : Removes all numerical property values from
                the object dictionary
                **'model_data'** : Removes only numerical data that were
                produced by an associated model
                **'labels'** : Removes all labels from the object
                dictionary, except those relating to the pore and throat
                locations of associated objects
                **'all'** : Removes both 'props' and 'labels'

        Notes
        -----
        If you wish to selectively remove some properties but not others,
        use something like ``del object['pore.blah']`` at the Python
        prompt. This can also be done in a for-loop to remove a list of
        items.

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
        element : str
            Can be either 'pore' or 'throat', which limits the returned list of
            keys to only 'pore' or 'throat' keys.  If neither is given, then
            both are assumed.
        mode : str, optional
            Controls which keys are returned. Options are:

                **'labels'**
                    Limits the returned list of keys to only 'labels'
                    (boolean arrays)
                **'props'**
                    Limits he return list of keys to only 'props'
                    (numerical arrays).
                **'all'**
                    Returns both 'labels' and 'props'. This is equivalent
                    to sending a list of both 'labels' and 'props'.

            If no mode is specified then the normal KeysView object is
            returned.
        deep : bool
            If set to ``True`` then the keys on all associated subdomain
            objects are returned as well.

        See Also
        --------
        props
        labels

        Notes
        -----
        This subclass can be used to get dictionary keys of specific
        kinds of data.  It's use augments ``props`` and ``labels`` by
        returning a list containing both types, but possibly limited by
        element type ('pores' or 'throats'.)

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
        element : str, optional
            Can be either 'pore' or 'throat' to specify what properties are
            returned.  If no element is given, both are returned
        mode : str, optional
            Controls what type of properties are returned.  Options are:
                **'all'** : Returns all properties on the object (default)

                **'models'** : Returns only properties that are associated with a
                model

                **'constants'** : returns data values that were *not* generated by
                a model, but manaully created.
        deep : bool
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

    @property
    def Np(self):
        r"""A shortcut to query the total number of pores on the object"""
        return np.shape(self.get('pore.all'))[0]

    @property
    def Nt(self):
        r"""A shortcut to query the total number of throats on the object"""
        return np.shape(self.get('throat.all'))[0]

    @property
    def Ps(self):
        r"""A shortcut to get a list of all pores on the object"""
        return np.arange(0, self.Np)

    @property
    def Ts(self):
        r"""A shortcut to get a list of all throats on the object"""
        return np.arange(0, self.Nt)

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

    def to_mask(self, pores=None, throats=None):
        r"""
        Convert a list of pore or throat indices into a boolean mask of the
        correct length.

        Parameters
        ----------
        pores : array_like
            List of pore indices. Only one of these can be specified at a
            time, and the returned result will be of the corresponding
            length.
        throats : array_like
            List of throat indices. Only one of these can be specified at
            a time, and the returned result will be of the corresponding
            length.

        Returns
        -------
        ndarray
            A boolean mask of length Np or Nt with ``True`` in the
            specified pore or throat locations.

        See Also
        --------
        to_indices

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> mask = pn.to_mask(pores=[0, 10, 20])
        >>> sum(mask)  # 3 non-zero elements exist in the mask (0, 10 and 20)
        3
        >>> len(mask)  # Mask size is equal to the number of pores in network
        125
        >>> mask = pn.to_mask(throats=[0, 10, 20])
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

    def to_indices(self, mask):
        r"""
        Converts a boolean mask to a list of pore or throat indices.

        Parameters
        ----------
        mask : array_like of booleans
            A boolean array with ``True`` at locations where indices are
            desired. The appropriate indices are returned based an the length
            of mask, which must be either Np or Nt long.

        Returns
        -------
        ndarray
            A list of pore or throat indices corresponding the locations
            where the received mask was ``True``.

        See Also
        --------
        to_mask

        Notes
        -----
        This behavior could just as easily be accomplished by using the mask
        in ``pn.pores()[mask]`` or ``pn.throats()[mask]``.  This method is
        just a convenience function and is a complement to ``to_mask``.

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
        prop : str
            The property name to be retrieved

        Returns
        -------
        ndarray
            A full length (Np or Nt) array of requested property values.

        Notes
        -----
        This makes an effort to maintain the data 'type' when possible;
        however when data are missing this can be tricky. Data can be
        missing in two different ways: A set of pores is not assisgned to
        a geometry or the network contains multiple geometries and data
        does not exist on all. Float and boolean data is fine, but missing
        ints are converted to float when nans are inserted.

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[2, 2, 2])
        >>> Ps = pn['pore.top']
        >>> Ts = pn.find_neighbor_throats(pores=Ps)
        >>> g1 = op.geometry.GenericGeometry(network=pn, pores=Ps, throats=Ts)
        >>> Ts = ~pn.to_mask(throats=Ts)
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
        # Fetch subdomains list depending on type of self
        proj = self.project
        if self._isa() in ['network', 'geometry']:
            subdomains = list(proj.geometries().values())
        elif self._isa() in ['phase', 'physics']:
            subdomains= list(proj.find_physics(phase=self))
        elif self._isa() in ['algorithm', 'base']:
            subdomains= [self]
        else:
            raise Exception('Unrecognized object type, cannot find dependents')

        # Get generalized element and array length
        element = self._parse_element(prop.split('.')[0], single=True)
        N = self.project.network._count(element)

        # Attempt to fetch the requested array from each object
        arrs = [obj.get(prop, None) for obj in subdomains]

        # Check for missing sources, and add None to arrs if necessary
        if N > sum([obj._count(element) for obj in subdomains]):
            arrs.append(None)

        # Obtain list of locations for inserting values
        locs = [self._get_indices(element, item.name) for item in subdomains]

        if np.all([item is None for item in arrs]):  # prop not found anywhere
            raise KeyError(prop)

        # Let's start by handling the easy cases first
        if not any([a is None for a in arrs]):
            # If any arrays are more than 1 column, then make they all are
            if np.any(np.array([a.ndim for a in arrs]) > 1):
                arrs = [np.atleast_2d(a).T if a.ndim == 1 else a for a in arrs]
            # All objs are present and array found on all objs
            try:
                W = max([a.shape[1] for a in arrs])  # Width of array
                shape = [N, W]
            except IndexError:
                shape = [N, ]
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

        return temp_arr

    def interpolate_data(self, propname, mode='mean'):
        r"""
        Determines a pore (or throat) property as the average of it's
        neighboring throats (or pores)

        Parameters
        ----------
        propname: str
            The dictionary key to the values to be interpolated.
        mode : str
            The method used for interpolation.  Options are 'mean' (default),
            'min', and 'max'.

        Returns
        -------
        vals : ndarray
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

    def get_conduit_data(self, poreprop, throatprop=None, mode='mean'):
        r"""
        Combines requested data into a single 3-column array.

        Parameters
        ----------
        poreprop : str
            The dictionary key to the pore property of interest
        throatprop : str, optional
            The dictionary key to the throat property of interest. If not
            given then the same property as ``poreprop`` is assumed.  So
            if poreprop = 'pore.foo' (or just 'foo'), then throatprop is
            set to 'throat.foo').
        mode : string
            How interpolation should be peformed for missing values. If values
            are present for both pores and throats, then this argument is
            ignored.  The ``interpolate`` data method is used.  Options are:

                * 'mean' (default)

                **'mean'** (default):
                    Finds the mean value of the neighboring pores (or throats)
                * 'min'
                    Finds the minimum of the neighboring pores (or throats)
                * 'max'
                    Finds the maximum of the neighboring pores (or throats)

        Returns
        -------
        conduit_data : ndarray
            An Nt-by-3 array with each column containing the requested
            property for each pore-throat-pore conduit.

        """
        # Deal with various args
        if not poreprop.startswith('pore'):
            poreprop = 'pore.' + poreprop
        if throatprop is None:
            throatprop = 'throat.' + poreprop.split('.', 1)[1]
        if not throatprop.startswith('throat'):
            throatprop = 'throat.' + throatprop
        # Generate array
        conns = self.network.conns
        domain = self._domain
        try:
            T = domain[throatprop]
            try:
                P1, P2 = domain[poreprop][conns.T]
            except KeyError:
                P = domain.interpolate_data(propname=throatprop, mode=mode)
                P1, P2 = P[conns.T]
        except KeyError:
            P1, P2 = domain[poreprop][conns.T]
            T = domain.interpolate_data(propname=poreprop, mode=mode)
        mask = self.throats(to_global=True)
        return np.vstack((P1[mask], T[mask], P2[mask])).T

    def _count(self, element=None):
        r"""
        Returns a dictionary containing the number of pores and throats in
        the network, stored under the keys 'pore' or 'throat'

        Parameters
        ----------
        element : str, optional
            Can be either 'pore' , 'pores', 'throat' or 'throats', which
            specifies which count to return.

        Returns
        -------
        dict
            A dictionary containing the number of pores and throats under
            the 'pore' and 'throat' key respectively.

        See Also
        --------
        num_pores
        num_throats

        Notes
        -----
        The ability to send plurals is useful for some types of 'programmatic'
        access. For instance, the standard argument for locations is pores
        or throats. If these are bundled up in a **kwargs dict then you can
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
        temp = np.size(self.__getitem__(element+'.all'))
        return temp

    def show_hist(self,
                  props=['pore.diameter', 'throat.diameter', 'throat.length'],
                  bins=20, fontsize=14, **kwargs):
        r"""
        Shows a quick plot of key property distributions.

        Parameters
        ----------
        props : str or List[str]
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

    def _parse_indices(self, indices):
        r"""
        This private method accepts a list of pores or throats and returns a
        properly structured Numpy array of indices.

        Parameters
        ----------
        indices : int or array_like
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
        element : str or List[str]
            The element argument to check.  If is None is recieved, then a list
            containing both \'pore\' and \'throat\' is returned.
        single : bool (default is False)
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
        labels : str or List[str]
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
        mode : str or List[str]
            The mode(s) to be parsed
        allowed : List[str]
            A list containing the allowed modes.  This list is defined by the
            calling method.  If any of the received modes are not in the
            allowed list an exception is raised.
        single : bool (default is False)
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
        element = self._parse_element(element, single=True)
        return element + '.' + propname.split('.')[-1]

    def __str__(self):
        module = self.__module__
        module = ".".join([x for x in module.split(".") if not x.startswith("_")])
        cname = self.__class__.__name__
        horizontal_rule = '―' * 78
        lines = [horizontal_rule]
        lines.append(f"{module}.{cname} : {self.name}")
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
