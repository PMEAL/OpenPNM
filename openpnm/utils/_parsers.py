import numpy as np


__all__ = [
    'parse_indices',
    'parse_element',
    'parse_labels',
    'parse_mode',
    'parse_prop',
]


def parse_indices(obj, indices):
    r"""
    Accepts a list of pores or throats and returns a properly structured
    Numpy array of indices.

    Parameters
    ----------
    indices : int or array_like
        This argument can accept numerous different data types including
        boolean masks, integers and arrays.

    Returns
    -------
    A Numpy array of indices

    """
    if indices is None:
        indices = np.array([], ndmin=1, dtype=int)
    locs = np.array(indices, ndmin=1)
    # If boolean array, convert to indices
    if locs.dtype == bool:
        if np.size(locs) == obj.Np:
            locs = obj.Ps[locs]
        elif np.size(locs) == obj.Nt:
            locs = obj.Ts[locs]
        else:
            raise Exception('Mask of locations must be either '
                            + 'Np nor Nt long')
    locs = locs.astype(dtype=int)
    return locs


def parse_element(obj, element, single=False):
    r"""
    Parses the keyword \'element\'

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


def parse_labels(obj, labels, element):
    r"""
    Converts \'labels\' to a proper format, including dealing with wildcards (\*).

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
            Ls = [L.split('.', 1)[-1] for L in obj.labels(element=element)]
            if label.startswith('*'):
                temp = [L for L in Ls if L.endswith(label.strip('*'))]
            if label.endswith('*'):
                temp = [L for L in Ls if L.startswith(label.strip('*'))]
            temp = [element+'.'+L for L in temp]
        elif element+'.'+label in obj.keys():
            temp = [element+'.'+label]
        else:
            temp = [element+'.'+label]
        parsed_labels.extend(temp)
        # Remove duplicates if any
        _ = [parsed_labels.remove(L) for L in parsed_labels
             if parsed_labels.count(L) > 1]
    return parsed_labels


def parse_mode(obj, mode, allowed=None, single=False):
    r"""
    Checks the \'mode\' argument

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


def parse_prop(obj, propname, element):
    element = obj._parse_element(element, single=True)
    if propname.split('.', 1)[0] in ['pore', 'throat']:
        propname = propname.split('.', 1)[-1]
    return element + '.' + propname
