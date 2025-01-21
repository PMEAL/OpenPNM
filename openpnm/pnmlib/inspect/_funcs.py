import numpy as np
from openpnm.pnmlib.core import (
    _parse_indices,
    _parse_element,
    count,
    get_data,
    set_data,
    num_pores,
    num_throats,
    get_prop_data,
    get_label_data,
    fold_dict,
)
# from openpnm.pnmlib.models import get_model_args


__all__ = [
    # Printout
    "info",
    "data",
    "tree",
    "models",
    "model",
    # Basic Queries
    "get_props",
    "get_labels",
    "get_params",
    "get_attrs",
    # Data
    "stats",
    "props_to_df",
    # Utils
    "PrintableList",
    "PrintableDict",
]


def get_prefixes(target):
    r"""
    Get a list of all defined prefixes on the given target (i.e., 'pore', 'throat')
    """
    elements = set()
    for item in target.keys():
        if '.' in item:
            elements.add(item.split('.', 1)[0])
    return list(elements)


def get_attrs(target, return_values=False):
    r"""
    Returns a list of which attributes are defined on the target

    Parameters
    ----------
    target : dict
        The dictionary for which the list of attributes is desired. If a nested
        dictionary is supplied only the top level is scanned.
    return_values : bool
        If `True` then a dictionary of attributes-values is returned. Otherwise,
        just a list of attribute names is returned.

    Returns
    -------
    params : list or dict
        If `return_values` is `False` (default) a list of attribute names is
        returned.  If `return_values` is `True` then a dictionary of attribute
        names and values is returned.

    """
    d = {}
    for k, v in target.items():
        if k.startswith('attr.'):
            d[k.split('.', 1)[1]] = v
    if not return_values:
        return PrintableList(d.keys())
    else:
        return PrintableDict(d, key="Attribute")


def get_params(target, return_values=False):
    r"""
    Returns a list of which parameters are defined on the target

    Parameters
    ----------
    target : dict
        The dictionary for which the list of parameters is desired. If a nested
        dictionary is supplied only the top level is scanned.
    return_values : bool
        If `True` then a dictionary of parameters-values is returned. Otherwise,
        just a list of parameter names is returned.

    Returns
    -------
    params : list or dict
        If `return_values` is `False` (default) a list of parameter names is
        returned.  If `return_values` is `True` then a dictionary of parameter
        names and values is returned.

    """
    d = {}
    for k, v in target.items():
        if k.startswith('param.'):
            d[k.split('.', 1)[1]] = v
    if not return_values:
        return PrintableList(d.keys())
    else:
        return PrintableDict(d, key="Parameter")


class PrintableList(list):
    r"""
    Simple subclass of ``list`` that has nice printing. Only works flat lists.

    """
    def __str__(self):
        horizontal_rule = "―" * 78
        lines = [horizontal_rule]
        self.sort()
        for i, item in enumerate(self):
            lines.append("{0:<5s} : {1}".format(str(i + 1), item))
        lines.append(horizontal_rule)
        return "\n".join(lines)

    # def __repr__(self):  # pragma: no cover
    #     return self.__str__()


class PrintableDict(dict):
    r"""
    Simple subclass of ``dict`` that has nicer printing.

    """
    def __init__(self, *args, key="Key", value="Value", **kwargs):
        self._value = value
        self._key = key
        super().__init__(*args, **kwargs)

    # def __repr__(self):  # pragma: no cover
    #     return self.__str__()

    def __str__(self):
        header = "―" * 78
        lines = [header, "{0:<35s} {1}".format(self._key, self._value), header]
        for item in list(self.keys()):
            if item.startswith('_'):
                continue
            if isinstance(self[item], np.ndarray):
                lines.append("{0:<35s} {1}".format(item, np.shape(self[item])))
            else:
                lines.append("{0:<35s} {1}".format(item, self[item]))
        lines.append(header)
        return "\n".join(lines)


def get_props(target, element=['pore', 'throat', 'conduit']):
    r"""
    Retrieves a list of keys that contain numerical data (i.e. "properties")

    Parameters
    ----------
    element : str, list of strings
        Indicates whether `'pore'` or `'throat'` properties should be returned.
        The default is `['pore', 'throat']`, so both are returned.

    Returns
    -------
    props : list of strings
        The names of all dictionary keys on the object that contain
        numerical data.
    """
    props = []
    for k, v in target.items():
        if k.split('.', 1)[0] in ['pore', 'throat', 'conduit']:
            el, prop = k.split('.', 1)
            if (v.dtype != bool) and not prop.startswith('_'):
                props.append(k)
    if len(props) == 0:
        for el in element:
            if el in target.keys():
                for k, v in target[el].items():
                    if (not hasattr(v, 'keys')) and (v.dtype != bool) \
                            and not k.startswith('_'):
                        props.append(el + '.' + k)
    props = sorted(props)
    props = PrintableList(props)
    return props


def _labels(target, element, locations, mode):
    r"""
    This is the actual label getter method, but it should not be called
    directly.  Use ``labels`` instead.
    """
    # Parse inputs
    locations = _parse_indices(target, locations)
    element = _parse_element(target, element=element)
    # Collect list of all pore OR throat labels
    labels = [i for i in target.keys(mode='labels') if i.split('.', 1)[0] in element]
    labels.sort()
    labels = np.array(labels)  # Convert to ndarray for following checks
    # Make an 2D array with locations in rows and labels in cols
    arr = np.vstack([target[item][locations] for item in labels]).T
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


def get_labels(
    target,
    pores=[],
    throats=[],
    element=['pore', 'throat'],
    mode='union'
):
    r"""
    Returns a list of labels present on the object

    Additionally, this function can return labels applied to a specified
    set of pores or throats

    Parameters
    ----------
    element : str
        Controls whether pore or throat labels are returned.  If empty then
        both are returned (default).
    pores (or throats) : array_like
        The pores (or throats) whose labels are sought.  If left empty a
        list containing all pore and throat labels is returned.
    mode : str, optional
        Controls how the query should be performed.  Only applicable
        when ``pores`` or ``throats`` are specified:

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

    """
    # Short-circuit query when no pores or throats are given
    if (np.size(pores) == 0) and (np.size(throats) == 0):
        labels = PrintableList()
        for k, v in target.items():
            if k.startswith('pore.') or k.startswith('throat.'):
                el, prop = k.split('.', 1)
                if (v.dtype == bool) and not prop.startswith('_'):
                    labels.append(k)
        if len(labels) == 0:
            for el in element:
                if el in target.keys():
                    for k, v in target[el].items():
                        if (not hasattr(v, 'keys')) and (v.dtype == bool) \
                                and not k.startswith('_'):
                            labels.append(el + '.' + k)
    elif (np.size(pores) > 0) and (np.size(throats) > 0):
        raise Exception('Cannot perform label query on pores and '
                        + 'throats simultaneously')
    elif np.size(pores) > 0:
        labels = _labels(target=target, element='pore', locations=pores,
                         mode=mode)
    elif np.size(throats) > 0:
        labels = _labels(target=target, element='throat', locations=throats,
                         mode=mode)
    return labels


def get_printable_props(target, suffix='', hr=78*'―'):
    r"""
    This function is used by the __str__ methods on all classes to get a
    nicely formatted list of properties on the object.

    Parameters
    ----------
    target : dict
        The dictionary object with each dictionary key containing a
        numpy array
    suffix : str, optional
        If provided, this will be attached to the end of every dictionary
        key so that 'pore.viscosity' becomes 'pore.viscosity.phase_01'.  This
        is a workaround to enhance the printing of component information on
        mixtures.
    hr : str, optional
        The horizontal rule to use between the table heading and body

    Returns
    -------
    table : str
        A formatted string that will output a 78 character wide table when
        printed

    Notes
    -----
    The table returned by this function only contains items that are numerical
    arrays.  Any boolean arrays are ignored.

    See Also
    --------
    get_printable_labels

    """
    if suffix and not suffix.startswith('.'):
        suffix = '.' + suffix
    header = [' ']*78
    header[2] = '#'
    header[5:15] = 'Properties'
    header[-12:] = 'Valid Values'
    lines = ''.join(header) + '\n' + hr
    i = 0
    for k in get_props(target):
        v = target[k]
        if (v.dtype != bool) and not ('._' in k):
            i += 1
            s = [' ']*78
            s[:3] = str(i+1).rjust(3)
            prop = k + suffix
            s[5:5+len(prop)] = prop
            element = k.split('.', 1)[0]
            nans = np.any(np.isnan(np.atleast_2d(v.T)), axis=0)
            valid = str(np.sum(~nans)) + ' / ' + str(count(target, element))
            s[-20:] = valid.rjust(20)
            a = ''.join(s)
            lines = '\n'.join((lines, a))
    return lines


def get_printable_labels(target, suffix='', hr=78*'―'):
    r"""
    This function is used by the __str__ methods on all classes to get a
    nicely formatted list of labels on the object.

    Parameters
    ----------
    target : dict
        The dictionary object with each dictionary key containing a
        numpy array
    suffix : str, optional
        If provided, this will be attached to the end of every dictionary
        key so that 'pore.viscosity' becomes 'pore.viscosity.phase_01'.  This
        is a workaround to enhance the printing of component information on
        mixtures.
    hr : str, optional
        The horizontal rule to use between the table heading and body

    Returns
    -------
    table : str
        A formatted string that will output a 78 character wide table when
        printed

    Notes
    -----
    The table returned by this function only contains items that boolean
    arrays.  Any numerical arrays are ignored.

    See Also
    --------
    get_printable_props
    """
    if suffix and not suffix.startswith('.'):
        suffix = '.' + suffix
    header = [' ']*78
    header[2] = '#'
    header[5:11] = 'Labels'
    header[-18:] = 'Assigned Locations'
    lines = ''.join(header) + '\n' + hr
    i = 0
    for k in get_labels(target):
        v = target[k]
        if (v.dtype == bool) and not ('._' in k):
            i += 1
            s = [' ']*78
            s[:3] = str(i+1).rjust(3)
            prop = k + suffix
            s[5:5+len(prop)] = prop
            valid = str(np.sum(v))
            s[-12:] = valid.rjust(12)
            a = ''.join(s)
            lines = '\n'.join((lines, a))
    return lines


def info(target):  # pragma: no cover
    hr = '―' * 78
    lines = ''
    lines += '\n' + "═"*78 + '\n'
    lines += get_printable_props(target)
    lines += '\n' + hr + '\n'
    lines += get_printable_labels(target)
    lines += '\n' + hr
    # return lines
    print(lines)


def data(target):
    info(target)


def models(models):
    horizontal_rule = '―' * 85
    lines = [horizontal_rule]
    strg = '{0:<3s} {1:<35s} {2:<25s} {3}'
    lines.append(strg.format('#', 'Property Name', 'Parameter', 'Value'))
    lines.append(horizontal_rule)
    for i, item in enumerate(models.keys()):
        temp = get_model_args(models[item])
        lines.append(strg.format(str(i+1), item, 'model:', item))
        for param in temp.keys():
            lines.append(strg.format('', '', param+':', temp[param]))
        lines.append(horizontal_rule)
    return '\n'.join(lines)


def model(model):
    horizontal_rule = '―' * 78
    lines = [horizontal_rule]
    strg = '{0:<25s} {1}'
    lines.append(strg.format('Parameter', 'Value'))
    lines.append(horizontal_rule)
    temp = get_model_args(model)
    for param in temp.keys():
        lines.append(strg.format(param+':', temp[param]))
    lines.append(horizontal_rule)
    return '\n'.join(lines)


def stats(network):
    hr = '―' * 78
    errmsg = '--------- ! --------- ! ---------\n'
    pad = 45
    lines = hr + '\n'
    # ---
    lines += 'Number of pores | throats: '.ljust(pad)
    lines += f"{num_pores(network)} | {num_throats(network)} \n"
    # ---
    lines += 'Pore Coordination (min | mean | max): '.ljust(pad)
    try:
        x = network['pore.coordination_number']
        lines += f"{x.min():.3E} | {x.mean():.3E} | {x.max():.3E} \n"
    except KeyError:
        lines += errmsg
    # ---
    lines += 'Pore Diameter (min | mean | max): '.ljust(pad)
    try:
        x = network['pore.diameter']
        lines += f"{x.min():.3E} | {x.mean():.3E} | {x.max():.3E} \n"
    except KeyError:
        lines += errmsg
    # ---
    lines += 'Pore Spacing (min | mean | max): '.ljust(pad)
    try:
        x = network['throat.spacing']
        lines += f"{x.min():.3E} | {x.mean():.3E} | {x.max():.3E} \n"
    except KeyError:
        lines += errmsg
    # ---
    lines += 'Throat Diameter (min | mean | max): '.ljust(pad)
    try:
        x = network['throat.diameter']
        lines += f"{x.min():.3E} | {x.mean():.3E} | {x.max():.3E} \n"
    except KeyError:
        lines += errmsg
    # ---
    lines += 'Throat Length (min | mean | max): '.ljust(pad)
    try:
        x = network['throat.length']
        lines += f"{x.min():.3E} | {x.mean():.3E} | {x.max():.3E} \n"
    except KeyError:
        lines += errmsg
    # ---
    lines += 'Throat Volume (min | mean | max): '.ljust(pad)
    try:
        x = network['throat.volume']
        lines += f"{x.min():.3E} | {x.mean():.3E} | {x.max():.3E} \n"
    except KeyError:
        lines += errmsg
    # ---
    lines += hr
    return lines


def tree(target, level=100, hide=[]):
    # print('*')
    temp = fold_dict(target)
    if hide is None:
        hide = []
    _tree(temp, level=0, maxdepth=level, hide=hide)


def _tree(target, level, maxdepth, hide):
    arrs = [k for k in sorted(target.keys()) if not hasattr(target[k], 'keys')
            and 'data' not in hide]
    keys = [k for k in target.keys() if hasattr(target[k], 'keys') and k not in hide]

    if level < maxdepth:
        _print_arrays(target, arrs, level=level, last=len(keys) == 0)
        _print_dicts(target, keys, level=level, maxdepth=maxdepth, hide=hide)


def _print_dicts(target, keys, level, maxdepth, hide):
    for i, k in enumerate(keys):
        bottom = (i == (len(keys) - 1))
        margin = '│' if (level > 0) else ''
        node = '└──' if bottom else '├──'
        space = '   │'*level
        space = space[:-1]
        print(margin + space + node + ' ' + k)
        _tree(target[k], level=level+1, maxdepth=maxdepth, hide=hide)


def _print_arrays(target, arrs, level, last):
    L = len(arrs)
    for i, arr in enumerate(arrs):
        margin = "│"*(level > 0)
        node = "└──" if (i == (L-1)) and last else "├──"
        space = '   │'*level
        space = space[:-1]
        if isinstance(target[arr], np.ndarray):
            print(margin + space + node, f"'{arr}':",
                  target[arr].shape, target[arr].dtype)
        else:
            print(margin + space + node, f"'{arr}':", target[arr])


def props_to_df(target, element):
    r"""
    Creates a pandas DataFrame from the specified data on the target

    Parameters
    ----------
    target : dict
        The dictionary with the data to be converted
    element : str
        The prefix specifying which data to extract (i.e. 'pore', 'throat', etc)

    Returns
    -------
    df : DataFrame
        A pandas DataFrame with each property as a column and each row containing
        the value of the that property in the corresponding element (i.e. pore index)
    """
    from pandas import DataFrame
    df = {}
    data = get_data(target, f'{element}.*')
    for k, arr in data.items():
        if not hasattr(arr, 'keys'):
            if np.ndim(arr) == 0:
                arr = [arr]
            elif np.ndim(arr) > 1:
                arr = arr.tolist()
            elif arr.dtype != bool:
                arr = np.around(arr, 3)
            df[k] = arr
    df = DataFrame(df)
    return df
