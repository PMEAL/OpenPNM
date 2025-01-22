import fnmatch
import numpy as np
from copy import deepcopy, copy
from collections.abc import Iterable



__all__ = [
    # Parsing functions
    "_parse_indices",
    "_parse_element",
    "_parse_labels",
    "_parse_mode",
    "_parse_prop",
    # Getters
    "get_data",
    "get_param_data",
    "get_prop_data",
    "get_label_data",
    # Setters
    "set_data",
    "set_label",
    # Indices
    "get_pores",
    "pores",
    "get_throats",
    "throats",
    "get_indices",
    "filter_by_label",
    "num_throats",
    "num_pores",
    "count",
    # Helpers
    "build_conduit_data",
    "interpolate_data",
    "flatten_dict",
    "flatten_list",
    "fold_dict",
]


def get_param_data(target, key=None):
    r"""
    Returns a dictionary containing only scalar numerical values
    """
    if key is None:
        key = '*.*'
    d = get_data(target, key)  # Get ALL data that matches key, then trim it down
    params = {}
    for k, v in d.items():
        if np.isscalar(v):
            params[k] = v
    return params


def get_prop_data(target, element=['pore', 'throat', 'conduit']):
    r"""
    Returns a dictionary containing only numerical arrays of the specified element(s)
    """
    if isinstance(element, str):
        element = [element]
    d = get_data(target, key)
    props = {}
    for k, v in d.items():
        el, prop = k.split('.', 1)
        if (el in element) and (v.dtype != bool):
            props[k] = v
    return props


def get_label_data(target, element=['pore', 'throat', 'conduit']):
    r"""
    Returns a dictionary containing only booolean arrays
    """
    if isinstance(element, str):
        element = [element]
    d = get_data(target, '*')
    labels = {}
    for k, v in d.items():
        el, prop = k.split('.', 1)
        if (el in element) and (v.dtype == bool):
            labels[k] = v
    return labels


def flatten_dict(group, delim='/'):
    r"""
    Takes a nested `dict` and returns a flattened version with the specified
    delimiter between keys

    Parameters
    ----------
    group : dict
        A nested dictionary
    delim : str
        The character to use when joining nested keys
    """
    # Taken from this SO answer https://stackoverflow.com/a/64717285, which has
    # LOTS of other suggestions.
    stack = list(group.items())
    ans = {}
    while stack:
        key, val = stack.pop()
        if isinstance(val, dict):
            for sub_key, sub_val in val.items():
                stack.append((f"{key + delim + sub_key}", sub_val))
        else:
            ans[key] = val
    return ans


def _merge_dicts(d1, d2):
    for k, v in d1.items():
        if isinstance(v, dict):
            _merge_dicts(v, d2.setdefault(k, {}))
        else:
            d2[k] = v
    return d2


def _expand_dicts(d, k, v, delim='/'):
    if '/' in k:
        group, k = k.split(delim, 1)
        d[group] = _expand_dicts({}, k, v)
    else:
        d[k] = v
    return d


def fold_dict(group, delim='/'):
    d = {}
    for k, v in group.items():
        d = _merge_dicts(d, _expand_dicts({}, k, v, delim=delim))
    return d


def flatten_list(input_list):
    r"""
    Given a list of nested lists of arbitrary depth, returns a single
    level or 'flat' list.
    """
    def _flatten(L):
        for el in L:
            if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
                yield from _flatten(el)
            else:
                yield el

    return list(_flatten(input_list))


def get_data(target, key):
    r"""
    Retrieves data from the target dictionary utilizing `pnmlib`'s syntax

    Parameters
    ----------
    target : dict
        The dictionary containing the simulation data
    key : str
        The name of the arrays to be retrieved. There are a number of syntactic
        shortcuts available:
        * If `key` is a numerical value it is returned directly, unchanged.
        * If `key` is blank, then target is returned directly.
        * If `key` contains `@<label>` then only values where `<label> == True` are
          returned
        * If `key` contains `*` then it is treated as a wildcard and a dictionary
          containing all arrays which match the `key` are returned.
        * If `key` contains `|`, the keys in the returned dictionary will have the
          text preceeding the '|' trimmed.

    Returns
    -------
    data : depends
        The returned data will depend on what was requested. Could be a numpy
        array or a dictionary of arrays.  Could also be a scalar or other object
        if parameters or attributes are request.
    """

    # Return non-string key directly
    if not isinstance(key, str):
        return key

    # Return specific key directly
    if key in target.keys():
        return target.get(key)

    # Return entire target if key is wildcardy
    if key in ['', '*', '/']:
        return target

    # Strip leading and trailing '/'
    key = key.strip('/').rstrip('/')

    # If key contains @ then return a subset of values at the requested locations
    if '@' in key:
        key, domain = key.split('@')
        try:
            group, prop = key.rsplit('/', 1)
        except ValueError:
            group, prop = '', key
        element, prop = prop.split('.', 1)
        vals = get_data(target, group + '/' + element + '.' + prop)
        try:  # Try domain on group
            locs = get_data(target, group + '/' + element + '.' + domain)
        except KeyError:
            try:  # Look for domain on project
                locs = get_data(target, element + '.' + domain)
            except KeyError:
                msg = f"{element + '.' + domain} not found in group or project"
                raise KeyError(msg)
        return vals[locs]

    # Slice the group name or prefix off of the keys
    if '|' in key:
        pre, post = key.split('|', 1)
        new_key = key.replace('|', '/')
        d = get_data(target, new_key)
        d = {k.split('/', pre.strip('/').count('/')+1)[-1]: v for k, v in d.items()}
        return d

    # Deal with wildcard in key
    if '*' in key:
        data = {}  # Dict to collect all hits
        hits = fnmatch.filter(target.keys(), key)
        for k in hits:
            data[k] = get_data(target, k)
        return data

    raise KeyError(key)


def set_data(target, key, value, locs=...):
    r"""
    Writes data to the target dictionary according to `pnmlib`'s rules

    Parameters
    ----------
    target : dict
        The dictionary containing the simulation data
    key : str
        The dictionary key to use when storing the data.
    value : scalar or array-like
        The numerical value(s) to write. If a scalar then it is applied everywhere,
        if a vector then it must be the correct size (i.e. $N_p-long$ or $N_t-long$),
        or its size must match `locs` if given. If `value` is a `dict` then each
        keywork is joined with `key` using a `/` meaning that `key` is the group.
    locs : ndarray or tuple
        The locations to write `value`. If `value` is a 1D array, then `locs` can
        be an array-like object of indices. If `value` is >1D, then `locs` must be
        a tuple of index arrays similar to that returned by `np.where`, meaning the
        first array in the tuple specifies the the first axis of `value` and so on.

    Returns
    -------
    This function operates "in-place" so nothing is returned.

    """
    from pnmlib import _reserved_prefixes as reserved_prefixes
    from pnmlib import _reserved_delimiters as reserved_delimiters

    # Trim leading and trailing /'s
    key = key.strip('/').rstrip('/')

    # If value is a dict, treat each one separately
    if hasattr(value, 'keys'):
        for k, v in value.items():
            set_data(target, key+'/'+k, v, locs=locs)
        return

    # If value is None then delete the given array
    if value is None:
        try:  # If value not present, silently return
            d = get_data(target, key)
        except KeyError:
            return
        if hasattr(d, 'keys'):
            for k, v in d.items():
                _ = target.pop(k)
        else:
            _ = target.pop(key)
        return

    # Intercept @ symbol
    if '@' in key:
        key, domain = key.split('@')
        if '/' in key:
            group, prop = key.rsplit('/', 1)
        else:
            group, prop = '', key
        element, prop = prop.split('.', 1)
        try:  # If prop is already on target then get it
            vals = get_data(target, group + '/' + element + '.' + prop)
        except KeyError:  # If not, create it
            vals = _initialize_empty_array_like(target, value, element)
            set_data(target, group + '/' + element + '.' + prop, vals, locs)
        try:  # If label exist on current group then use it
            _locs = get_data(target, group+'/'+element+'.'+domain)
        except KeyError:  # If not, check top level group
            _locs = get_data(target, element+'.'+domain)
        # Write values specified locations in array
        vals[_locs] = value
        # Finally write array to dict
        set_data(target, group + '/' + element + '.' + prop, vals, locs)
        return

    # Make sure group name, if present, is valid and retrieve it
    if '/' in key:
        groups = key.split('/')
        hits = set(groups).intersection(set(reserved_prefixes))
        if len(hits) > 0:
            raise Exception(f'{list(hits)[0]} is a reserved prefix, cannot use as a group')
        for item in groups[:-1]:
            if '.' in item:
                raise Exception('Group name cannot contain dots')
        group, name = key.rsplit('/', 1)
    else:
        group, name = '', key

    # Parse property name
    try:
        element, prop = name.split('.', 1)
    except ValueError:
        raise Exception(f'{name} is not a valid property name')

    # Make sure element is valid
    if element not in reserved_prefixes:
        raise Exception(f"Unsupported prefix '{element}'")

    # Intercept param prefix, ensure they're scalar numbers, and write
    if element == 'param':
        value = np.array(value).squeeze()
        if (np.size(value) == 1) and isinstance(value.item(), (int, float)):
            target.update({key: value.item()})
            return
        else:
            raise Exception("The 'param' prefix only accepts scalar values")

    # Intercept attr prefix, ensure they're serializable (not implemented), and write
    if element == 'attr':
        target.update({key: value})
        return

    # Confirm shape of conduit arrays
    if element == 'conduit':
        if (np.size(value) > 1) and (locs == ...):  # if not a scalar
            if (value.shape[1] != 3):
                raise Exception("conduit arrays must be 3 columns wide")

    # Convert value to ndarray for use below
    if hasattr(value, '__iter__'):
        value = np.array(list(value), ndmin=1)
    elif np.size(value) == 1:
        value = np.array(value, ndmin=1)

    # Try determining the lengh of element arrays for use below
    try:
        size = count(target, element=element)
    except Exception:  # If size cant be determined, given array will become default
        if value.shape[0] > 1:  # Ensure array is not a scalar meant to be broadcast
            target.update({key: value})
            return
        else:
            msg = f"Unable to determine size of array to hold {key}"
            raise Exception(msg)

    # Finally write arrays
    if np.size(value) == 1:  # Broadcast scalars to full length
        if key in target.keys():
            temp = target.get(key)
        else:
            temp = _initialize_empty_array_like(target, value, element)
        temp[locs] = value
        target.update({key: temp})
    elif np.shape(value)[0] == size:
        target.update({key: value})  # Array is already correct length
    elif locs != ...:  # Locs were provided
        temp = target.get(key)
        if temp is None:
            temp = _initialize_empty_array_like(target, value, element)
        temp[locs] = value
        target.update({key: temp})
    else:
        raise Exception('Provided array is wrong length for ' + key)


def set_label(target, label, pores=None, throats=None, mode='add'):
    r"""
    Creates or updates a label array

    Parameters
    ----------
    label : str
        The label to apply to the specified locations
    pores : array_like
        A list of pore indices or a boolean mask of where given label
        should be added or removed (see ``mode``)
    throats : array_like
        A list of throat indices or a boolean mask of where given label
        should be added or removed (see ``mode``)
    mode : str
        Controls how the labels are handled.  Options are:

        =========== ======================================================
        mode        description
        =========== ======================================================
        'add'       (default) Adds the given label to the specified
                    locations while keeping existing labels

        'overwrite' Removes existing label from all locations before
                    adding the label in the specified locations

        'remove'    Removes the given label from the specified locations
                    leaving the remainder intact

        'purge'     Removes the specified label from the object completely.
                    This ignores the ``pores`` and ``throats`` arguments.

        'clear'     Sets all the labels to ``False`` but does not remove
                    the label array
        =========== ======================================================

    """
    _parse_mode(
        target=target,
        mode=mode,
        allowed=['add', 'overwrite', 'remove', 'purge', 'clear'],
    )

    if label.split('.', 1)[0] in ['pore', 'throat']:
        label = label.split('.', 1)[1]

    if (pores is not None) and (throats is not None):
        set_label(target=target, label=label, pores=pores, mode=mode)
        set_label(target=target, label=label, throats=throats, mode=mode)
        return
    elif pores is not None:
        locs = _parse_indices(target, pores)
        element = 'pore'
        N = num_pores(target)
    elif throats is not None:
        locs = _parse_indices(target, throats)
        element = 'throat'
        N = num_throats(target)

    blank = np.zeros(N, dtype=bool)
    if mode == 'add':
        if element + '.' + label not in target.keys():
            target[element + '.' + label] = blank
        target[element + '.' + label][locs] = True
    if mode == 'overwrite':
        target[element + '.' + label] = blank
        target[element + '.' + label][locs] = True
    if mode == 'remove':
        target[element + '.' + label][locs] = False
    if mode == 'clear':
        target['pore' + '.' + label] = False
        target['throat' + '.' + label] = False
    if mode == 'purge':
        _ = target.pop('pore.' + label, None)
        _ = target.pop('throat.' + label, None)


def count(target, element):
    if '/' in element:
        element = element.rsplit('/', 1)[1]
    element = 'throat' if element == 'conduit' else element
    counts = set()
    for k, v in target.items():
        if element in k:
            try:
                counts.add(np.shape(v)[0])
            except IndexError:
                counts.add(np.size(v))
    if len(counts) == 0:
        raise Exception(f'No {element} arrays found')
    elif len(counts) == 1:
        return list(counts)[0]
    else:
        raise Exception('Multiple arrays with different length')


def get_indices(target, element, labels='all', mode='or'):
    r"""

    """
    # Parse and validate all input values.
    element = _parse_element(target=target, element=element, single=True)
    labels = _parse_labels(target=target, labels=labels, element=element)

    # Begin computing label array
    if mode in ['or', 'any', 'union']:
        union = np.zeros([count(target, element), ], dtype=bool)
        for item in labels:  # Iterate over labels and collect all indices
            union = union + target[element+'.'+item.split('.', 1)[-1]]
        ind = union
    elif mode in ['and', 'all', 'intersection']:
        intersect = np.ones([count(target, element), ], dtype=bool)
        for item in labels:  # Iterate over labels and collect all indices
            intersect = intersect*target[element+'.'+item.split('.', 1)[-1]]
        ind = intersect
    elif mode in ['xor', 'exclusive_or']:
        xor = np.zeros([count(target, element), ], dtype=int)
        for item in labels:  # Iterate over labels and collect all indices
            info = target[element+'.'+item.split('.', 1)[-1]]
            xor = xor + np.int8(info)
        ind = (xor == 1)
    elif mode in ['nor', 'not', 'none']:
        nor = np.zeros([count(target, element), ], dtype=int)
        for item in labels:  # Iterate over labels and collect all indices
            info = target[element+'.'+item.split('.', 1)[-1]]
            nor = nor + np.int8(info)
        ind = (nor == 0)
    elif mode in ['nand']:
        nand = np.zeros([count(target, element), ], dtype=int)
        for item in labels:  # Iterate over labels and collect all indices
            info = target[element+'.'+item.split('.', 1)[-1]]
            nand = nand + np.int8(info)
        ind = (nand < len(labels)) * (nand > 0)
    elif mode in ['xnor', 'nxor']:
        xnor = np.zeros([count(target, element), ], dtype=int)
        for item in labels:  # Iterate over labels and collect all indices
            info = target[element+'.'+item.split('.', 1)[-1]]
            xnor = xnor + np.int8(info)
        ind = (xnor > 1)
    else:
        raise Exception('Unsupported mode: '+mode)
    # Extract indices from boolean mask
    ind = np.where(ind)[0]
    ind = ind.astype(dtype=int)
    return ind


def get_pores(target, labels='all', mode='or', asmask=False):
    r"""
    Returns pore indices where given labels exist, according to the logic
    specified by the ``mode`` argument.

    Parameters
    ----------
    labels : str or list[str]
        The label(s) whose pores locations are requested.  This argument
        also accepts '*' for wildcard searches.
    mode : str
        Specifies how the query should be performed.  The options are:

        ==============  ===================================================
        mode            meaning
        ==============  ===================================================
        'or'            Returns the labels that are assigned to *any* of
                        the given locations. Also accepts 'union' and 'any'
        'and'           Labels that are present on all the given locations.
                        also accepts 'intersection' and 'all'
        'xor'           Labels that are present on *only one*
                        of the given locations.Also accepts 'exclusive_or'
        'nor'           Labels that are *not* present on any of
                        the given locations. Also accepts 'not' and 'none'
        'nand'          Labels that are present on *all but one* of the
                        given locations
        'xnor'          Labels that are present on *more than one* of the
                        given locations.
        ==============  ===================================================

    asmask : bool
        If ``True`` then a boolean array of length Np is returned with
        ``True`` values indicating the pores that satisfy the query.

    Returns
    -------
    A Numpy array containing pore indices filtered by the logic specified
    in ``mode``.

    See Also
    --------
    throats

    Notes
    -----
    Technically, *nand* and *xnor* should also return pores with *none* of
    the labels but these are not included.  This makes the returned list
    more useful.

    To perform more complex or compound queries, you can opt to receive
    the result a a boolean mask (``asmask=True``), then manipulate the
    arrays manually.

    """
    ind = get_indices(target=target, element='pore', labels=labels, mode=mode)
    if asmask:
        ind = to_mask(target=target, pores=ind)
    return ind


def get_throats(target, labels='all', mode='or', asmask=False):
    r"""
    Returns throat locations where given labels exist, according to the
    logic specified by the ``mode`` argument.

    Parameters
    ----------
    labels : str or list[str]
        The throat label(s) whose locations are requested.  If omitted,
        'all' throat indices are returned.  This argument also accepts
        '*' for wildcard searches.
    mode : str
        Specifies how the query should be performed. The options are:

        ==============  ===================================================
        mode            meaning
        ==============  ===================================================
        'or'            Returns the labels that are assigned to *any* of
                        the given locations. Also accepts 'union' and 'any'
        'and'           Labels that are present on all the given locations.
                        also accepts 'intersection' and 'all'
        'xor'           Labels that are present on *only one*
                        of the given locations.Also accepts 'exclusive_or'
        'nor'           Labels that are *not* present on any of
                        the given locations. Also accepts 'not' and 'none'
        'nand'          Labels that are present on *all but one* of the
                        given locations
        'xnor'          Labels that are present on *more than one* of the
                        given locations.
        ==============  ===================================================

    asmask : bool
        If ``True`` then a boolean array of length Nt is returned with
        ``True`` values indicating the throats that satisfy the query.

    Returns
    -------
    A Numpy array containing throat indices filtered by the logic specified
    in ``mode``.

    See Also
    --------
    pores

    """
    ind = get_indices(target=target, element='throat', labels=labels, mode=mode)
    if asmask:
        ind = to_mask(target=target, throats=ind)
    return ind


def pores(target, labels='all', mode='or', asmask=False):
    Ps = get_pores(target=target, labels=labels, mode=mode, asmask=asmask)
    return Ps


def throats(target, labels='all', mode='or', asmask=False):
    Ts = get_throats(target=target, labels=labels, mode=mode, asmask=asmask)
    return Ts


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
    mode : str
        Controls how the filter is applied. The default value is
        'or'. Options include:

        ==============  ===================================================
        mode            meaning
        ==============  ===================================================
        'or'            Returns the labels that are assigned to *any* of
                        the given locations. Also accepts 'union' and 'any'
        'and'           Labels that are present on all the given locations.
                        also accepts 'intersection' and 'all'
        'xor'           Labels that are present on *only one*
                        of the given locations.Also accepts 'exclusive_or'
        'nor'           Labels that are *not* present on any of
                        the given locations. Also accepts 'not' and 'none'
        'nand'          Labels that are present on *all but one* of the
                        given locations
        'xnor'          Labels that are present on *more than one* of the
                        given locations.
        ==============  ===================================================

    Returns
    -------
    A list of pores (or throats) that have been filtered according the
    given criteria. The returned list is a subset of the received list of
    pores (or throats).

    See Also
    --------
    pores
    throats

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
        return np.array([], dtype=int)
    labels = self._parse_labels(labels=labels, element=element)
    labels = [element+'.'+item.split('.', 1)[-1] for item in labels]
    all_locs = self._get_indices(element=element, labels=labels, mode=mode)
    mask = self._tomask(indices=all_locs, element=element)
    ind = mask[locations]
    return locations[ind]


def num_pores(target, labels='all', mode='or'):
    r"""
    Returns the number of pores of the specified labels

    Parameters
    ----------
    labels : list of strings, optional
        The pore labels that should be included in the count.
        If not supplied, all pores are counted.
    labels : list of strings
        Label of pores to be returned
    mode : str, optional
        Specifies how the count should be performed. The options are:

        ==============  ===================================================
        mode            meaning
        ==============  ===================================================
        'or'            Returns the labels that are assigned to *any* of
                        the given locations. Also accepts 'union' and 'any'
        'and'           Labels that are present on all the given locations.
                        also accepts 'intersection' and 'all'
        'xor'           Labels that are present on *only one*
                        of the given locations.Also accepts 'exclusive_or'
        'nor'           Labels that are *not* present on any of
                        the given locations. Also accepts 'not' and 'none'
        'nand'          Labels that are present on *all but one* of the
                        given locations
        'xnor'          Labels that are present on *more than one* of the
                        given locations.
        ==============  ===================================================

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

    """
    # Count number of pores of specified type
    Ps = get_indices(target=target, labels=labels, mode=mode, element='pore')
    Np = np.shape(Ps)[0]
    return Np


def num_throats(target, labels='all', mode='union'):
    r"""
    Return the number of throats of the specified labels

    Parameters
    ----------
    labels : list of strings, optional
        The throat labels that should be included in the count.
        If not supplied, all throats are counted.
    mode : str, optional
        Specifies how the count should be performed.  The options are:

        ==============  ===================================================
        mode            meaning
        ==============  ===================================================
        'or'            Returns the labels that are assigned to *any* of
                        the given locations. Also accepts 'union' and 'any'
        'and'           Labels that are present on all the given locations.
                        also accepts 'intersection' and 'all'
        'xor'           Labels that are present on *only one*
                        of the given locations.Also accepts 'exclusive_or'
        'nor'           Labels that are *not* present on any of
                        the given locations. Also accepts 'not' and 'none'
        'nand'          Labels that are present on *all but one* of the
                        given locations
        'xnor'          Labels that are present on *more than one* of the
                        given locations.
        ==============  ===================================================

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
    Ts = get_indices(target=target, labels=labels, mode=mode, element='throat')
    Nt = np.shape(Ts)[0]
    return Nt


def to_mask(target, pores=None, throats=None):
    r"""
    Generates a boolean mask with `True` values in the given locations

    Parameters
    ----------
    pores : array_like
        The pore indices where `True` values will be placed. If `pores` is
        given the `throats` is ignored.
    throats : array_like
        The throat indices where `True` values will be placed. If `pores` is
        given the `throats` is ignored.

    Returns
    -------
    mask : ndarray, boolean
        A boolean array of length Np is `pores` was given or Nt if
        `throats` was given.

    """
    if pores is not None:
        indices = np.array(pores, ndmin=1)
        N = num_pores(target)
    elif throats is not None:
        indices = np.array(throats, ndmin=1)
        N = num_throats(target)
    else:
        raise Exception('Must specify either pores or throats')
    mask = np.zeros((N, ), dtype=bool)
    mask[indices] = True
    return mask


def to_indices(target, mask):
    r"""
    Converts a boolean mask to pore or throat indices

    Parameters
    ----------
    mask : ndarray
        A boolean mask with `True` values indicating either pore or
        throat indices. This array must either be Nt or Np long, otherwise
        an Exception is raised.

    Returns
    -------
    indices : ndarray
        An array containing numerical indices of where `mask` was `True`.

    Notes
    -----
    This function is equivalent to just calling `np.where(mask)[0]` but
    does check to ensure that `mask` is a valid length.
    """
    mask = np.array(mask, dtype=bool)
    if mask.shape[0] not in [num_pores(target), num_throats(target)]:
        raise Exception('Mask must be either Nt or Np long')
    return np.where(mask)[0]


def interpolate_data(target, propname, mode='mean'):
    r"""
    Generates an array of the requested pore/throat data by interpolating
    the neighboring throat/pore data.

    Parameters
    ----------
    propname : str
        The data to be generated.
    mode : str
        Dictate how the interpolation is done. Options are 'mean', 'min',
        and 'max'.

    Returns
    -------
    data : ndarray
        An ndarray containing the interpolated data.  E.g. Requesting
        'throat.temperature' will read the values of 'pore.temperature'
        in each of the neighboring pores and compute the average
        (if `mode='mean'`).
    """
    from openpnm.models.misc import from_neighbor_throats, from_neighbor_pores
    element, prop = propname.split('.', 1)
    if element == 'throat':
        if target['pore.'+prop].dtype == bool:
            raise Exception('The requested datatype is boolean, cannot interpolate')
        values = from_neighbor_pores(target, prop='pore.'+prop, mode=mode)
    elif element == 'pore':
        if target['throat.'+prop].dtype == bool:
            raise Exception('The requested datatype is boolean, cannot interpolate')
        values = from_neighbor_throats(target, prop='throat.'+prop, mode=mode)
    return values


def _initialize_empty_array_like(target, value, key):
    value = np.array(value)
    size = count(target, element=key.split('.', 1)[0])
    if value.dtype == bool:
        temp = np.zeros([size, *value.shape[1:]],
                        dtype=bool)
    elif key.split('.', 1)[0] in ['conduit']:
        temp = np.zeros([size, 3],
                        dtype=float)*np.nan
    else:
        temp = np.zeros([size, *value.shape[1:]],
                        dtype=float)*np.nan
    return temp


def build_conduit_data(target, network, propname):
    r"""
    Fetches an Nt-by-3 array of the requested property

    Parameters
    ----------
    propname : str
        The dictionary key of the property to fetch.

    Returns
    -------
    data : ndarray
        An Nt-by-3 array with each column containing the requrested data
        for pore1, throat, and pore2 respectively.

    """
    poreprop = 'pore.' + propname.split('.', 1)[-1]
    throatprop = 'throat.' + propname.split('.', 1)[-1]
    conns = network['throat.conns']
    try:
        T = target[throatprop]
        if T.ndim > 1:
            raise Exception(f'{throatprop} must be a single column wide')
    except KeyError:
        T = np.ones([num_throats(target), ], dtype=float)*np.nan
    try:
        P1, P2 = target[poreprop][conns.T]
    except KeyError:
        P1 = np.ones([num_throats(target), ], dtype=float)*np.nan
        P2 = np.ones([num_throats(target), ], dtype=float)*np.nan
    vals = np.vstack((P1, T, P2)).T
    if np.isnan(vals).sum() == vals.size:
        raise KeyError(f'{propname} not found')
    return vals


def _parse_indices(target, indices):
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
        if np.size(locs) == num_pores(target):
            locs = pores(target)[locs]
        elif np.size(locs) == num_throats(target):
            locs = throats(target)[locs]
        else:
            raise Exception('Mask of locations must be either '
                            + 'Np nor Nt long')
    locs = locs.astype(dtype=int)
    return locs


def _parse_element(target, element, single=False):
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
    element = [item.split('.', 1)[0] for item in element]
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


def _parse_labels(target, labels, element):
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
            label = label.split('.', 1)[-1]
        # Deal with wildcards
        if '*' in label:
            Ls = [L.split('.', 1)[-1] for L in labels(target, element=element)]
            if label.startswith('*'):
                temp = [L for L in Ls if L.endswith(label.strip('*'))]
            if label.endswith('*'):
                temp = [L for L in Ls if L.startswith(label.strip('*'))]
            temp = [element+'.'+L for L in temp]
        elif element+'.'+label in target.keys():
            temp = [element+'.'+label]
        else:
            temp = [element+'.'+label]
        parsed_labels.extend(temp)
        # Remove duplicates if any
        _ = [parsed_labels.remove(L) for L in parsed_labels
             if parsed_labels.count(L) > 1]
    return parsed_labels


def _parse_mode(target, mode, allowed=None, single=False):
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


def _parse_prop(target, propname, element):
    element = _parse_element(target, element, single=True)
    if propname.split('.', 1)[0] in ['pore', 'throat']:
        propname = propname.split('.', 1)[-1]
    return element + '.' + propname


if __name__ == "__main__":
    from openpnm.pnmlib.inspect import tree
    d1 = {
        'pore.all': np.ones(10, dtype=bool),
        'throat.all': np.ones(10, dtype=bool),
        'param.test': 2.2,
        'phase1': {
            'pore.all': np.ones(10, dtype=bool),
            'pore.test1': np.ones(10, dtype=int),
            'throat.all': np.ones(20, dtype=bool),
            'throat.test1': np.ones(20, dtype=float),
            'param.test1': 2.2,
            'phase3': {
                'pore.all': np.ones(10, dtype=bool),
                'pore.test3': np.ones(10, dtype=int),
                'throat.all': np.ones(20, dtype=bool),
                'throat.test3': np.ones(20, dtype=float),
                'param.test3': 2.2,
            },
        },
        'phase2': {
            'pore.all': np.ones(10, dtype=bool),
            'pore.test2': np.ones(10, dtype=int),
            'throat.all': np.ones(20, dtype=bool),
            'throat.test2': np.ones(20, dtype=float),
            'param.test2': 2.2,
            'phase3': {
                'pore.all': np.ones(10, dtype=bool),
                'pore.test4': np.ones(10, dtype=int),
                'throat.all': np.ones(20, dtype=bool),
                'throat.test4': np.ones(20, dtype=float),
                'param.test4': 2.2,
            },
        },
    }
    d2 = flatten_dict(d1)
    key = '*.all'
    print(key)
    print('─'*10)
    tree(get_data(target=d2, key=key))
    print('─'*10)
    tree(get_data(target=d2, key='phase1/*.*'))

    # def test_