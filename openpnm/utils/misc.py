import copy
import inspect
import warnings
import functools
import numpy as _np
import scipy as _sp
import time as _time
from collections import OrderedDict
from docrep import DocstringProcessor


class Docorator(DocstringProcessor):

    __instance__ = None

    def __new__(cls, *args, **kwargs):
        if Docorator.__instance__ is None:
            Docorator.__instance__ = DocstringProcessor()
        return Docorator.__instance__


class PrintableList(list):
    r"""
    Simple subclass of ``list`` that has nice printing.  Only works flat lists.

    Example
    -------
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

    def __repr__(self):
        self.sort()
        return super().__repr__()


class PrintableDict(OrderedDict):
    r"""
    Simple subclass of ``dict`` that has nicer printing.

    Example
    -------
    >>> from openpnm.utils import PrintableDict
    >>> from numpy import array as arr
    >>> d = {'item1': 1, 'item2': '1', 'item3': [1, 1], 'item4': arr([1, 1])}
    >>> print(PrintableDict(d))
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
    key                                 value
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
    item1                               1
    item2                               1
    item3                               [1, 1]
    item4                               (2,)
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――

    If the item is a Numpy array the value column will contain the items'
    shape, otherwise it will contain the result of ``print(item)``

    """

    def __init__(self, *args, **kwargs):
        self._value = "value"
        self._key = "key"
        super().__init__(*args, **kwargs)

    def __repr__(self):
        text = dict(self).__str__()
        return text

    def __str__(self):
        header = "―" * 78
        lines = [header, "{0:<35s} {1}".format(self._key, self._value), header]
        for item in list(self.keys()):
            if type(self[item]) == _sp.ndarray:
                lines.append("{0:<35s} {1}".format(item, _np.shape(self[item])))
            else:
                lines.append("{0:<35s} {1}".format(item, self[item]))
        lines.append(header)
        return "\n".join(lines)


class SettingsDict(PrintableDict):
    r"""
    The SettingsDict implements the __missing__ magic method, which returns
    None instead of KeyError.  This is useful for checking the value of a
    settings without first ensuring it exists.

    Examples
    --------
    >>> from openpnm.utils import SettingsDict
    >>> sd = SettingsDict()
    >>> sd['test'] = True
    >>> print(sd['test'])
    True
    >>> print(sd['not_a_valid_key'])
    None

    """
    __doc__ = ''

    def __setitem__(self, key, value):
        if hasattr(value, "Np"):
            raise Exception(
                "Cannot store OpenPNM objects in settings, "
                + "store object's name instead"
            )
        super().__setitem__(key, value)

    def __missing__(self, key):
        self[key] = None
        return self[key]

    def _update_settings_and_docs(self, dc):
        if type(dc) is type:  # If dc is class then instantiate it
            dc = dc()
        self.__doc__ = dc.__doc__
        # if dc is a dataclass object.  This step is only necessary to support
        # Python 3.6 which doesn't have the dataclasses module
        if hasattr(dc, '__dict__'):
            dc = copy.deepcopy(dc.__dict__)
        else:
            dc = copy.deepcopy(dc)
        for item in dc.keys():
            self[item] = dc[item]


class GenericSettings:
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for item in dir(self):
            if not item.startswith('__'):
                self.__dict__[item] = getattr(self, item)


class SubDict(dict):
    def __getitem__(self, key):
        for item in self.keys():
            if item.endswith('.' + key):
                key = item
        return super().__getitem__(key)


class NestedDict(dict):
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
    This class adds a 'health' check to a standard dictionary.  This check
    looks into the dict values, and considers empty lists as healthy and all
    else as unhealthy.  If one or more entries is 'unhealthy' the health method
    returns False.

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


def tic():
    r"""
    Homemade version of matlab tic and toc function, tic starts or resets
    the clock, toc reports the time since the last call of tic.

    See Also
    --------
    toc

    """
    global _startTime_for_tictoc
    _startTime_for_tictoc = _time.time()


def toc(quiet=False):
    r"""
    Homemade version of matlab tic and toc function, tic starts or resets
    the clock, toc reports the time since the last call of tic.

    Parameters
    ----------
    quiet : Boolean
        If False (default) then a message is output to the console.  If True
        the message is not displayed and the elapsed time is returned.

    See Also
    --------
    tic

    """
    if "_startTime_for_tictoc" in globals():
        t = _time.time() - _startTime_for_tictoc
        if quiet is False:
            print(f"Elapsed time in seconds: {t:0.2f}")
        else:
            return t
    else:
        raise Exception("Start time not set, call tic first")


def unique_list(input_list):
    r"""
    For a given list (of points) remove any duplicates

    """
    output_list = []
    if len(input_list) > 0:
        dim = _np.shape(input_list)[1]
        for i in input_list:
            match = False
            for j in output_list:
                if dim == 3:
                    if i[0] == j[0] and i[1] == j[1] and i[2] == j[2]:
                        match = True
                elif dim == 2:
                    if i[0] == j[0] and i[1] == j[1]:
                        match = True
                elif dim == 1:
                    if i[0] == j[0]:
                        match = True
            if match is False:
                output_list.append(i)
    return output_list


def flat_list(input_list):
    r"""
    Given a list of nested lists of arbitrary depth, returns a single level or
    'flat' list.

    """
    x = input_list
    if isinstance(x, list):
        return [a for i in x for a in flat_list(i)]
    else:
        return [x]


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
    Converts a ModelsDict object to a ReST compatible table

    Parameters
    ----------
    obj : OpenPNM object
        Any object that has a ``models`` attribute

    params : boolean
        Indicates whether or not to include a list of parameter
        values in the table.  Set to False for just a list of models, and
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
                    p1 = p1[:14] + "..."
                p2 = str(temp[param])
                if len(p2) > 24:
                    p2 = p2[:21] + "..."
                lines.append(fmt.format("|", "", "|", "", "|", p1, "|", p2, "|"))
                lines.append(row)
    return "\n".join(lines)


def catch_module_not_found(function):
    """
    A decorator that wraps the passed in function and catches
    ModuleNotFound exception.
    """
    @functools.wraps(function)
    def wrapper(*args, **kwargs):
        try:
            return function(*args, **kwargs)
        except ModuleNotFoundError:
            pass
    return wrapper


def ignore_warnings(warning=RuntimeWarning):
    r"""
    Decorator for catching warnings. Useful in pore-scale models where nans
    are inevitable, and numpy gets annoying by throwing lots of RuntimeWarnings.

    Parameters
    ----------
    warning : Python Warning object
        Python warning type that you want to temporarily ignore

    Examples
    --------
    >>> from openpnm.utils.misc import ignore_warnings
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


def conduit_lengths(network, throats=None, mode="pore"):
    r"""
    Return the respective lengths of the conduit components defined by the throat
    conns P1 - T - P2

    Notes
    -----
    mode = 'pore' - uses pore coordinates
    mode = 'centroid' uses pore and throat centroids

    """
    if throats is None:
        throats = network.throats()
    Ps = network["throat.conns"]
    pdia = network["pore.diameter"]
    Lt = network["throat.length"]

    if mode == "centroid":
        try:
            pcentroids = network["pore.centroid"]
            tcentroids = network["throat.centroid"]
            if _np.sum(_np.isnan(pcentroids)) + _np.sum(_np.isnan(tcentroids)) > 0:
                mode = "pore"
            else:
                plen1 = (
                    _np.sqrt(_np.sum(_sp.square(pcentroids[Ps[:, 0]] - tcentroids), 1))
                    - Lt / 2
                )
                plen2 = (
                    _np.sqrt(_np.sum(_sp.square(pcentroids[Ps[:, 1]] - tcentroids), 1))
                    - Lt / 2
                )
        except KeyError:
            mode = "pore"
    if mode == "pore":
        # Find half-lengths of each pore
        pcoords = network["pore.coords"]
        # Find the pore-to-pore distance, minus the throat length
        lengths = (
            _np.sqrt(_np.sum(_sp.square(pcoords[Ps[:, 0]] - pcoords[Ps[:, 1]]), 1)) - Lt
        )
        lengths[lengths < 0.0] = 2e-9
        # Calculate the fraction of that distance from the first pore
        try:
            fractions = pdia[Ps[:, 0]] / (pdia[Ps[:, 0]] + pdia[Ps[:, 1]])
            # Don't allow zero lengths
            # fractions[fractions == 0.0] = 0.5
            # fractions[fractions == 1.0] = 0.5
        except Exception:
            fractions = 0.5
        plen1 = lengths * fractions
        plen2 = lengths * (1 - fractions)

    return _np.vstack((plen1, Lt, plen2)).T[throats]


def is_symmetric(a, rtol=1e-10):
    r"""
    Is ``a`` a symmetric matrix?

    Parameters
    ----------
    a : ndarray, sparse matrix
        Object to check for being a symmetric matrix.

    rtol : float
        Relative tolerance with respect to the smallest entry in ``a`` that
        is used to determine if ``a`` is symmetric.

    Returns
    -------
    bool
        ``True`` if ``a`` is a symmetric matrix, ``False`` otherwise.

    """
    if type(a) != _sp.ndarray and not _sp.sparse.issparse(a):
        raise Exception("'a' must be either a sparse matrix or an ndarray.")
    if a.shape[0] != a.shape[1]:
        raise Exception("'a' must be a square matrix.")

    atol = _np.amin(_np.absolute(a.data)) * rtol
    if _sp.sparse.issparse(a):
        issym = False if ((a - a.T) > atol).nnz else True
    elif type(a) == _sp.ndarray:
        issym = False if _np.any((a - a.T) > atol) else True

    return issym


def is_valid_propname(propname):
    r"""
    Check if ``propname`` is a valid OpenPNM propname, i.e. starts with
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
    if type(propname) is not str:
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


def nbr_to_str(nbr, t_precision):
    r"""
    Converts a scalar into a string in scientific (exponential) notation
    without the decimal point.

    Parameters
    ----------
    nbr : scalar
        The number to be converted into a scalar.

    t_precision : integer
        The time precision (number of decimal places). Default value is 12.

    """
    from decimal import Decimal as dc
    n = int(-dc(str(round(nbr, t_precision))).as_tuple().exponent
            * (round(nbr, t_precision) != int(nbr)))
    nbr_str = (str(int(round(nbr, t_precision)*10**n)) + ('e-'+str(n))*(n != 0))
    return nbr_str
