import functools
import inspect
import os
import warnings
from collections.abc import Iterable

import numpy as np
import numpy.lib.recfunctions as rf
import scipy.sparse as sparse
from docrep import DocstringProcessor

__all__ = [
    'Docorator',
    'PrintableList',
    'PrintableDict',
    'HealthDict',
    'NestedDict',
    'flat_list',
    'sanitize_dict',
    'methods_to_table',
    'models_to_table',
    'ignore_warnings',
    'is_symmetric',
    'is_valid_propname',
    'is_transient',
    'get_mixture_model_args',
    'dict_to_struct',
    'struct_to_dict',
    'get_printable_props',
    'get_printable_labels',
]


class Docorator(DocstringProcessor):
    """OpenPNM's customized docstring processor."""

    __instance__ = None

    def __new__(cls, *args, **kwargs):
        if Docorator.__instance__ is None:
            Docorator.__instance__ = DocstringProcessor()

        # Add custom parameter type sections
        a = DocstringProcessor.param_like_sections
        Docorator.__instance__.param_like_sections = a + []  # ["Attributes", "Settings"]
        # Add custom text type sections
        a = Docorator.__instance__.text_sections
        Docorator.__instance__.text_sections = a + []

        # Create a single list of all section types
        a = Docorator.__instance__.param_like_sections
        b = Docorator.__instance__.text_sections
        Docorator.__instance__.all_sections = a + b

        return Docorator.__instance__


class PrintableList(list):
    r"""
    Simple subclass of ``list`` that has nice printing. Only works flat lists.

    Examples
    --------
    >>> from openpnm.utils import PrintableList
    >>> temp = ['item1', 'item2', 'item3']
    >>> print(PrintableList(temp))
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
    1     : item1
    2     : item2
    3     : item3
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――

    Each line contains the result of ``print(item)`` on each item in the list

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

    Examples
    --------
    >>> from openpnm.utils import PrintableDict
    >>> from numpy import array as arr
    >>> d = {'item1': 1, 'item2': '1', 'item3': [1, 1], 'item4': arr([1, 1])}
    >>> print(PrintableDict(d))
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
    Key                                 Value
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
    item1                               1
    item2                               1
    item3                               [1, 1]
    item4                               (2,)
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――

    If the item is a Numpy array the value column will contain the items'
    shape, otherwise it will contain the result of ``print(item)``

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


class NestedDict(dict):
    """Brief explanation of 'NestedDict'"""

    def __init__(self, mapping={}, delimiter="/"):
        super().__init__()
        self.delimiter = delimiter
        self.update(mapping)
        self.unravel()

    def __setitem__(self, key, value):
        path = key.split(self.delimiter, 1)
        if len(path) > 1:
            if path[0] not in self.keys():
                self[path[0]] = NestedDict(delimiter=self.delimiter)
            self[path[0]][path[1]] = value
        else:
            super().__setitem__(key, value)

    def __missing__(self, key):
        self[key] = NestedDict(delimiter=self.delimiter)
        return self[key]

    def unravel(self):
        for item in self.keys():
            self[item] = self.pop(item)

    def to_dict(self, dct=None):
        if dct is None:
            dct = self
        plain_dict = dict()
        for key in dct.keys():
            value = dct[key]
            if hasattr(value, "keys"):
                plain_dict[key] = self.to_dict(value)
            else:
                plain_dict[key] = value
        return plain_dict

    def keys(self, dicts=True, values=True):
        k = list(super().keys())
        new_keys = []
        for item in k:
            if hasattr(self[item], "keys"):
                if dicts:
                    new_keys.append(item)
            else:
                if values:
                    new_keys.append(item)
        return new_keys

    def __str__(self):
        def print_level(self, p="", indent="-"):
            for item in self.keys():
                if hasattr(self[item], "keys"):
                    p = print_level(self[item], p=p, indent=indent + indent[0])
                elif indent[-1] != " ":
                    indent = indent + ""
                p = indent + item + "\n" + p
            return p

        p = print_level(self)
        return p


class HealthDict(PrintableDict):
    r"""
    This class adds a 'health' check to a standard dictionary.

    This check looks into the dict values, and considers empty lists as
    healthy and all else as unhealthy.  If one or more entries is
    'unhealthy' the health method returns False.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def _get_health(self):
        health = True
        for item in list(self.keys()):
            try:
                if len(self[item]) > 0:
                    health = False
            except TypeError:
                if self[item]:
                    health = False
        return health

    health = property(fget=_get_health)

    def __bool__(self):
        return self.health


def flat_list(input_list):
    r"""
    Given a list of nested lists of arbitrary depth, returns a single
    level or 'flat' list.
    """
    def _flatten(l):
        for el in l:
            if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
                yield from _flatten(el)
            else:
                yield el

    return list(_flatten(input_list))


def sanitize_dict(input_dict):
    r"""
    Given a nested dictionary, ensures that all nested dicts are normal
    Python dicts.  This is necessary for pickling, or just converting
    an 'auto-vivifying' dict to something that acts normal.
    """
    plain_dict = dict()
    for key in input_dict.keys():
        value = input_dict[key]
        if hasattr(value, "keys"):
            plain_dict[key] = sanitize_dict(value)
        else:
            plain_dict[key] = value
    return plain_dict


def methods_to_table(obj):
    r"""
    Converts a methods on an object to a ReST compatible table

    Parameters
    ----------
    obj : Base
        Any object that has a methods
    params : bool
        Indicates whether or not to include a list of parameter
        values in the table. Set to False for just a list of models, and
        True for a more verbose table with all parameter values.

    """
    parent = obj.__class__.__mro__[1]
    temp = inspect.getmembers(parent, predicate=inspect.isroutine)
    parent_funcs = [i[0] for i in temp if not i[0].startswith("_")]

    temp = inspect.getmembers(obj.__class__, predicate=inspect.isroutine)
    obj_funcs = [i[0] for i in temp if not i[0].startswith("_")]
    funcs = set(obj_funcs).difference(set(parent_funcs))

    row = "+" + "-" * 22 + "+" + "-" * 49 + "+"
    fmt = "{0:1s} {1:20s} {2:1s} {3:47s} {4:1s}"
    lines = []
    lines.append(row)
    lines.append(fmt.format("|", "Method", "|", "Description", "|"))
    lines.append(row.replace("-", "="))
    for i, item in enumerate(funcs):
        try:
            s = getattr(obj, item).__doc__.strip()
            end = s.find("\n")
            if end > 47:
                s = s[:44] + "..."
            lines.append(fmt.format("|", item, "|", s[:end], "|"))
            lines.append(row)
        except AttributeError:
            pass
    return "\n".join(lines)


def models_to_table(obj, params=True):
    r"""
    Converts a all the models on an object to a ReST compatible table

    Parameters
    ----------
    obj : Base
        Any object that has a ``models`` attribute
    params : bool
        Indicates whether or not to include a list of parameter
        values in the table. Set to False for just a list of models, and
        True for a more verbose table with all parameter values.

    """
    if not hasattr(obj, "models"):
        raise Exception("Received object does not have any models")
    row = "+" + "-" * 4 + "+" + "-" * 22 + "+" + "-" * 18 + "+" + "-" * 26 + "+"
    fmt = "{0:1s} {1:2s} {2:1s} {3:20s} {4:1s} {5:16s} {6:1s} {7:24s} {8:1s}"
    lines = []
    lines.append(row)
    lines.append(
        fmt.format("|", "#", "|", "Property Name", "|", "Parameter", "|", "Value", "|")
    )
    lines.append(row.replace("-", "="))
    for i, item in enumerate(obj.models.keys()):
        prop = item
        if len(prop) > 20:
            prop = item[:17] + "..."
        temp = obj.models[item].copy()
        model = str(temp.pop("model")).split(" ")[1]
        lines.append(
            fmt.format("|", str(i + 1), "|", prop, "|", "model:", "|", model, "|")
        )
        lines.append(row)
        if params:
            for param in temp.keys():
                p1 = param
                if len(p1) > 16:
                    p1 = p1[:13] + "..."
                p2 = str(temp[param])
                if len(p2) > 24:
                    p2 = p2[:21] + "..."
                lines.append(fmt.format("|", "", "|", "", "|", p1, "|", p2, "|"))
                lines.append(row)
    return "\n".join(lines)


def ignore_warnings(warning=RuntimeWarning):
    r"""
    Decorator for catching warnings. Useful in pore-scale models where
    nans are inevitable, and numpy gets annoying by throwing lots of
    RuntimeWarnings.

    Parameters
    ----------
    warning : Warning
        Python warning type that you want to temporarily ignore

    Examples
    --------
    >>> from openpnm.utils import ignore_warnings
    >>> @ignore_warnings()
    ... def myfun(x):
    ...     return 1/x

    >>> import numpy as np
    >>> x = np.arange(5)
    >>> myfun(x)
    array([       inf, 1.        , 0.5       , 0.33333333, 0.25      ])

    """

    def _ignore_warning(function):
        @functools.wraps(function)
        def __ignore_warning(*args, **kwargs):
            with warnings.catch_warnings(record=True):
                # Catch all warnings of this type
                warnings.simplefilter("always", warning)
                # Execute the function
                result = function(*args, **kwargs)
            return result

        return __ignore_warning

    return _ignore_warning


def is_symmetric(a, rtol=1e-10):
    r"""
    Is ``a`` a symmetric matrix?

    Parameters
    ----------
    a : ndarray or sparse matrix
        Object to check for being a symmetric matrix.
    rtol : float
        Relative tolerance with respect to the smallest entry in ``a``
        that is used to determine if ``a`` is symmetric.

    Returns
    -------
    bool
        ``True`` if ``a`` is a symmetric matrix, ``False`` otherwise.

    """
    if not isinstance(a, np.ndarray) and not sparse.issparse(a):
        raise Exception("'a' must be either a sparse matrix or an ndarray.")
    if a.shape[0] != a.shape[1]:
        raise Exception("'a' must be a square matrix.")

    atol = np.amin(np.absolute(a.data)) * rtol
    if sparse.issparse(a):
        issym = False if ((a - a.T) > atol).nnz else True
    elif isinstance(a, np.ndarray):
        issym = False if np.any((a - a.T) > atol) else True

    return issym


def get_mixture_model_args(
    phase,
    composition='xs',
    args={
        'mus': 'pore.viscosity',
        'MWs': 'param.molecular_weight',
    }
):
    r"""
    This is used in tests to run models generically
    """
    from openpnm.models.phase.misc import mole_to_mass_fraction
    vals = {}
    if composition in ['ws']:
        temp = np.vstack(list(mole_to_mass_fraction(phase=phase).values()))[:, 0]
        vals[composition] = temp
    else:
        temp = np.vstack(list(phase['pore.mole_fraction'].values()))[:, 0]
        vals[composition] = temp
    for item in args.keys():
        temp = np.vstack(list(phase.get_comp_vals(args[item]).values()))[:, 0]
        vals[item] = temp
    return vals


def dict_to_struct(d):
    r"""
    Converts a dictionary of numpy arrays to a numpy struct

    Parameters
    ----------
    d : dict
        A dictionary wtih numpy arrays in each key. The arrays must be all
        the same size.

    Returns
    -------
    s : numpy struct
        A numpy struct with the fields or names take from the dictionary keys
    """
    struct = rf.unstructured_to_structured(np.vstack(list(d.values())).T,
                                           names=list(d.keys()))
    return struct


def struct_to_dict(s):
    r"""
    Converts a numpy struct array into a dictionary using the struct labels as
    keys

    Parameters
    ----------
    s : numpy struct
        The struct array

    Returns
    -------
    d : dict
        A dictionary with the struct labels or fields as the keys
    """
    d = {}
    for key in s.dtype.names:
        d[key] = s[key]
    return d


def get_printable_props(item, suffix='', hr=78*'―'):
    r"""
    This function is used by the __str__ methods on all classes to get a
    nicely formatted list of properties on the object.

    Parameters
    ----------
    item : dict
        The OpenPNM dictionary object with each dictionary key containing a
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
    for k, v in item.items():
        if (v.dtype != bool) and not ('._' in k):
            i += 1
            s = [' ']*78
            s[:3] = str(i+1).rjust(3)
            prop = k + suffix
            s[5:5+len(prop)] = prop
            element = k.split('.', 1)[0]
            nans = np.any(np.isnan(np.atleast_2d(v.T)), axis=0)
            valid = str(np.sum(~nans)) + ' / ' + str(item._count(element))
            s[-20:] = valid.rjust(20)
            a = ''.join(s)
            lines = '\n'.join((lines, a))
    return lines


def get_printable_labels(item, suffix='', hr=78*'―'):
    r"""
    This function is used by the __str__ methods on all classes to get a
    nicely formatted list of labels on the object.

    Parameters
    ----------
    item : dict
        The OpenPNM dictionary object with each dictionary key containing a
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
    for k, v in item.items():
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


def is_transient(algorithms):
    # check that algorithms is a list
    if type(algorithms) is not list:
        algorithms = [algorithms]
    # return True if any algorithm is transient
    for alg in algorithms:
        if hasattr(alg, 'soln'):
            soln_type = type(alg.soln[alg.settings['quantity']])
            if 'TransientSolution' in str(soln_type):
                return True
    return False


def is_valid_propname(propname):
    r"""
    Checks if ``propname`` is a valid OpenPNM propname, i.e. starts with
    'pore.' or 'throat.'

    Parameters
    ----------
    propname : str
        Property name to check whether it's a valid OpenPNM propname.

    Returns
    -------
    bool
        Whether or not ``propname`` is a valid name

    """
    if not isinstance(propname, str):
        return False
    temp = propname.split(".")
    if temp[0] not in ["pore", "throat"]:
        return False
    if len(temp) == 1:
        return False
    for field in temp:
        if len(field) == 0:
            return False
    return True


def nbr_to_str(nbr, t_precision=12):
    r"""
    Converts a scalar into a string in scientific (exponential) notation
    without the decimal point.

    Parameters
    ----------
    nbr : scalar
        The number to be converted into a scalar.
    t_precision : integer
        The time precision (number of decimal places). Default value is 12.

    Returns
    -------
    num : str
        The string represenation of the given number in scientific notation
    """
    from decimal import Decimal as dc
    n = int(-dc(str(round(nbr, t_precision))).as_tuple().exponent
            * (round(nbr, t_precision) != int(nbr)))
    nbr_str = (str(int(round(nbr, t_precision)*10**n)) + ('e-'+str(n))*(n != 0))
    return nbr_str
