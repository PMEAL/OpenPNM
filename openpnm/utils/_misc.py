import sys
import math
import time
import inspect
import warnings
import functools
import numpy as np
import scipy.sparse as sparse
from collections import OrderedDict
from collections.abc import Iterable
from docrep import DocstringProcessor
from copy import deepcopy
import numpy.lib.recfunctions as rf


__all__ = [
    'Docorator',
    'PrintableList',
    'PrintableDict',
    'SubDict',
    'NestedDict',
    'HealthDict',
    'tic',
    'toc',
    'unique_list',
    'flat_list',
    'flat_list2',
    'sanitize_dict',
    'methods_to_table',
    'models_to_table',
    'catch_module_not_found',
    'ignore_warnings',
    'is_symmetric',
    'is_valid_propname',
    'prettify_logger_message',
    'remove_prop_deep',
    'get_model_collection',
    'dict_to_struct',
    'struct_to_dict',
    'get_mixture_model_args',
]


class Docorator(DocstringProcessor):

    __instance__ = None

    def __new__(cls, *args, **kwargs):
        if Docorator.__instance__ is None:
            Docorator.__instance__ = DocstringProcessor()

        # Add custom parameter type sections
        a = DocstringProcessor.param_like_sections
        Docorator.__instance__.param_like_sections = a + [] # ["Attributes", "Settings"]
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

    # def __repr__(self):
    #     return self.__str__()


class PrintableDict(OrderedDict):
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

    # def __repr__(self):
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


class SubDict(dict):
    """Brief explanation of 'SubDict'"""
    def __getitem__(self, key):
        for item in self.keys():
            if item.endswith('.' + key):
                key = item
        return super().__getitem__(key)


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


"""
BSD 3-Clause License

- Copyright (c) 2008-Present, IPython Development Team
- Copyright (c) 2001-2007, Fernando Perez <fernando.perez@colorado.edu>
- Copyright (c) 2001, Janko Hauser <jhauser@zscout.de>
- Copyright (c) 2001, Nathaniel Gray <n8gray@caltech.edu>

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
def _format_time(timespan, precision=3):
    """Formats the timespan in a human readable form"""

    if timespan >= 60.0:
        # we have more than a minute, format that in a human readable form
        # Idea from http://snipplr.com/view/5713/
        parts = [("d", 60*60*24),("h", 60*60),("min", 60), ("s", 1)]
        time = []
        leftover = timespan
        for suffix, length in parts:
            value = int(leftover / length)
            if value > 0:
                leftover = leftover % length
                time.append(u'%s%s' % (str(value), suffix))
            if leftover < 1:
                break
        return " ".join(time)


    # Unfortunately the unicode 'micro' symbol can cause problems in
    # certain terminals.
    # See bug: https://bugs.launchpad.net/ipython/+bug/348466
    # Try to prevent crashes by being more secure than it needs to
    # E.g. eclipse is able to print a µ, but has no sys.stdout.encoding set.
    units = [u"s", u"ms",u'us',"ns"] # the save value
    if hasattr(sys.stdout, 'encoding') and sys.stdout.encoding:
        try:
            u'\xb5'.encode(sys.stdout.encoding)
            units = [u"s", u"ms",u'\xb5s',"ns"]
        except:
            pass
    scaling = [1, 1e3, 1e6, 1e9]

    if timespan > 0.0:
        order = min(-int(math.floor(math.log10(timespan)) // 3), 3)
    else:
        order = 3
    return u"%.*g %s" % (precision, timespan * scaling[order], units[order])


def tic():
    r"""
    Homemade version of matlab tic and toc function, tic starts or resets
    the clock, toc reports the time since the last call of tic.

    See Also
    --------
    toc

    """
    global _startTime_for_tictoc
    _startTime_for_tictoc = time.time()


def toc(quiet=False):
    r"""
    Homemade version of matlab tic and toc function, tic starts or resets
    the clock, toc reports the time since the last call of tic.

    Parameters
    ----------
    quiet : bool, default is False
        If False then a message is output to the console. If
        True the message is not displayed and the elapsed time is returned.

    See Also
    --------
    tic

    """
    if "_startTime_for_tictoc" not in globals():
        raise Exception("Start time not set, call tic first")
    t = time.time() - _startTime_for_tictoc
    if quiet is False:
        print(f"Elapsed time: {_format_time(t)}")
    return t


def unique_list(input_list):
    r"""
    For a given list (of points) remove any duplicates
    """
    output_list = []
    if len(input_list) > 0:
        dim = np.shape(input_list)[1]
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
    Given a list of nested lists of arbitrary depth, returns a single
    level or 'flat' list.
    """
    x = input_list
    if isinstance(x, list):
        return [a for i in x for a in flat_list(i)]
    return [x]


def flat_list2(input_list):
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
    """Brief explanation of 'methods_to_table'"""
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
                    p1 = p1[:14] + "..."
                p2 = str(temp[param])
                if len(p2) > 24:
                    p2 = p2[:21] + "..."
                lines.append(fmt.format("|", "", "|", "", "|", p1, "|", p2, "|"))
                lines.append(row)
    return "\n".join(lines)


def catch_module_not_found(function):
    r"""
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


def prettify_logger_message(msg):
    """Prettifies logger messages by breaking them up into multi lines"""
    from textwrap import wrap
    linewidth = 75 - 13
    indent = "\n" + " " * 13
    temp = wrap(msg, width=linewidth)
    return indent.join(temp)


def remove_prop_deep(obj, propname):
    """Hierarchically deletes the given propname and its children"""
    for k in list(obj.keys()):
        obj._parse_element(propname)
        if k.startswith(propname):
            del obj[k]


def get_model_collection(collection, regen_mode=None, domain=None):
    d = deepcopy(collection)
    for k, v in d.items():
        if regen_mode:
            v['regen_mode'] = regen_mode
        if domain:
            v['domain'] = domain
    return d


def dict_to_struct(d):
    struct = rf.unstructured_to_structured(np.vstack(list(d.values())),
                                           names=list(d.keys()))
    return struct


def struct_to_dict(s):
    d = {}
    for key in s.dtype.names:
        d[key] = s[key]
    return d


def get_mixture_model_args(
    target,
    composition='xs',
    args={
        'mus': 'pore.viscosity',
        'MWs': 'param.molecular_weight',
    }
):
    from openpnm.models.phase.misc import mole_to_mass_fraction
    vals = {}
    if composition in ['ws']:
        temp = np.vstack(list(mole_to_mass_fraction(target=target).values()))[:, 0]
        vals[composition] = temp
    else:
        temp = np.vstack(list(target['pore.mole_fraction'].values()))[:, 0]
        vals[composition] = temp
    for item in args.keys():
        temp = np.vstack(list(target.get_comp_vals(args[item]).values()))[:, 0]
        vals[item] = temp
    return vals
