import logging
from copy import deepcopy

from openpnm.utils import PrintableDict

logger = logging.getLogger(__name__)


__all__ = [
    'TypedMixin',
    'TypedSet',
    'TypedList',
    'SettingsAttr',
]


class TypedMixin:
    """Based class for enforcing types on lists and sets."""

    def __init__(self, iterable=[], types=[]):
        self._types = types
        if iterable:
            super().__init__(iterable)
            self._set_types()

    def _get_types(self):
        if not hasattr(self, '_types'):
            self._types = []
        if self._types == []:
            self._types = list(set([type(i) for i in self]))
        return self._types

    def _set_types(self):
        if self._types == []:
            self._types = list(set([type(i) for i in self]))
        else:
            raise Exception("Types have already been defined")

    types = property(fget=_get_types, fset=_set_types)

    def _check_type(self, value):
        if (type(value) not in self.types) and (len(self.types) > 0):
            raise TypeError("This list cannot accept values of type "
                            + f"{type(value)}")


class TypedSet(TypedMixin, set):
    """A set that enforces all elements have the same type."""

    def add(self, item):
        self._check_type(item)
        super().add(item)


class TypedList(TypedMixin, list):
    """A list that enforces all elements have the same type."""

    def __setitem__(self, ind, value):
        self._check_type(value)
        super().__setitem__(ind, value)

    def append(self, value):
        self._check_type(value)
        super().append(value)

    def extend(self, iterable):
        for value in iterable:
            self._check_type(value)
        super().extend(iterable)

    def insert(self, index, value):
        self._check_type(value)
        super().insert(index, value)


class SettingsDict(dict):
    def update(self, d):
        for k, v in d:
            self[k] = v
        if "Parameters" in d.__doc__:
            pass

    def __getattr__(self, key):
        return self[key]

    def __setattr__(self, key, value):
        self[key] = value


class SettingsAttr:
    r"""
    A custom data class that holds settings for objects.

    The main function of this custom class is to enforce the datatype of
    values that are assigned to ensure they remain consistent.  For instance
    if ``obj.foo = "bar"``, then ``obj.foo = 456`` will fail.
    """

    def __init__(self, *args):
        for i, item in enumerate(args):
            if i == 0:
                super().__setattr__('__doc__', item.__doc__)
            self._update(item)

    def __setattr__(self, attr, value):
        if hasattr(self, attr):
            # If the the attr is already present, check its type
            if getattr(self, attr) is not None:
                # Ensure the written type is an instance of the existing one
                a = value.__class__.__mro__
                b = getattr(self, attr).__class__.__mro__
                c = object().__class__.__mro__
                check = list(set(a).intersection(set(b)).difference(set(c)))
                if len(check) > 0:  # If they share comment parent class
                    super().__setattr__(attr, value)
                else:  # Otherwise raise an error
                    old = type(getattr(self, attr))
                    new = type(value)
                    raise TypeError(f"Attribute \'{attr}\' can only accept "
                                    + f"values of type {old}, but the recieved"
                                    + f" value was of type {new}")
            else:
                # If the current attr is None, let anything be written
                super().__setattr__(attr, value)
        else:
            # If there is no current attr, let anything be written
            super().__setattr__(attr, value)

    def _update(self, settings, docs=False, override=False):
        if settings is None:
            return
        if isinstance(settings, dict):
            docs = False
            for k, v in settings.items():
                v = deepcopy(v)
                if override:
                    super().__setattr__(k, v)
                else:
                    setattr(self, k, v)
        else:  # Dataclass
            attrs = [i for i in dir(settings) if not i.startswith('_')]
            for k in attrs:
                v = deepcopy(getattr(settings, k))
                if override:
                    super().__setattr__(k, v)
                else:
                    setattr(self, k, v)
        if docs:
            self._getdocs(settings)

    @property
    def _attrs(self):
        a = dir(self)
        b = dir(list())
        attrs = list(set(a).difference(set(b)))
        attrs = [i for i in attrs if not i.startswith('_')]
        attrs = sorted(attrs)
        return attrs

    def _deepcopy(self):
        return deepcopy(self)

    def _getdocs(self, settings):
        super().__setattr__('__doc__', settings.__doc__)

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)

    def __str__(self):  # pragma: no cover
        d = PrintableDict(key="Settings", value="Values")
        d.update(self.__dict__)
        for item in self.__dir__():
            if not item.startswith('_'):
                d[item] = getattr(self, item)
        return d.__str__()

    def __repr__(self):  # pragma: no cover
        return self.__str__()
