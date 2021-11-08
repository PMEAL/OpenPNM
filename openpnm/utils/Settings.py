from traits.api import HasTraits, Trait
from openpnm.utils import PrintableDict
from copy import deepcopy


class TypedList(list):

    def __init__(self, iterable=[], types=[]):
        self._types = types
        if iterable:
            super().__init__(iterable)
            self._set_types()

    def _get_types(self):
        if self._types == []:
            self._types = list(set([type(i) for i in self]))
        return self._types

    def _set_types(self):
        if self._types == []:
            self._types = list(set([type(i) for i in self]))
        else:
            raise Exception("Types have already been defined")

    types = property(fget=_get_types, fset=_set_types)

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

    def _check_type(self, value):
        if (type(value) not in self.types) and (len(self.types) > 0):
            raise Exception("This list cannot accept values of type " +
                            f"{type(value)}")


class SettingsAttr:

    def __init__(self, settings={}):
        super().__setattr__('__doc__', settings.__doc__)
        self._update(settings)

    def __setattr__(self, attr, value):
        if hasattr(value, '__contains__'):
            value = deepcopy(value)
        if hasattr(self, attr):
            a = value.__class__.__mro__
            b = getattr(self, attr).__class__.__mro__
            c = object().__class__.__mro__
            check = list(set(a).intersection(set(b)).difference(set(c)))
            if len(check) > 0:
                super().__setattr__(attr, value)
            else:
                old = type(getattr(self, attr))
                new = type(value)
                raise Exception(f"Attribute \'{attr}\' can only accept " +
                                f"values of type {old}, but the recieved " +
                                f"value was of type {new}")
        else:
            super().__setattr__(attr, value)

    def _update(self, settings, docs=False, override=False):
        if hasattr(settings, 'items'): # Dictionary
            for k, v in settings.items():
                if override:
                    super().__setattr__(k, v)
                else:
                    setattr(self, k, v)
        else:  # Dataclass
            attrs = [i for i in dir(settings) if not i.startswith('_')]
            for k in attrs:
                v = getattr(settings, k)
                if override:
                    super().__setattr__(k, v)
                else:
                    setattr(self, k, v)
        if docs:
            super().__setattr__('__doc__', settings.__doc__)

    @property
    def _attrs(self):
        a = dir(self)
        b = dir(list())
        attrs = list(set(a).difference(set(b)))
        attrs = [i for i in attrs if not i.startswith('_')]
        return attrs

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)

    def __str__(self):
        d = PrintableDict()
        d._key = 'Settings'
        d._value = 'Values'
        d.update(self.__dict__)
        return d.__str__()

    def __repr__(self):
        return self.__str__()


class SettingsData(HasTraits):

    def __str__(self):
        d = PrintableDict()
        d._key = 'Settings'
        d._value = 'Values'
        for item in self.visible_traits():
            d[item] = getattr(self, item)
        return d.__str__()

    # def __repr__(self):
    #     return self.__str__()


class SettingsAttr2:
    r"""
    """

    def __init__(self, settings):
        r"""
        Parameters
        ----------
        settings : key-value collection
            Ideally should be an `SettingsData`` object which is a subclass
            of a the ``HasTraits`` class from the ``traits`` package. The
            docstring for this class is adopted from the received argument only
            **if** a ``HasTraits`` class is received. This is meant to work
            with ``docrep`` so the inherited docstring on the ``SettingsData``
            object will be shown by this wrapper. Optionally ``settings`` can
            be a ``dict`` or any ``dataclasss``-type object with data stored
            as attributes.

        Notes
        -----
        This class is a wrapper around the ``HasTraits`` class from the
        ``traits`` package. The main point is the provide a clean namespace
        on the object so that only the specified attributes show up.  It also
        allows for combining multiple objects into one using the ``_update``
        method.
        """
        # Add blank HasTraits object to self._settings
        if isinstance(settings, SettingsData):
            super().__setattr__('_settings', settings)
        else:
            super().__setattr__('_settings', SettingsData())
            self._update(settings)
        self._settings.__doc__ = settings.__doc__

    def __setattr__(self, attr, val):
        if attr in self._settings.visible_traits():
            setattr(self._settings, attr, val)
        else:
            val = Trait(val, val.__class__)
            self._settings.add_trait(attr, val)

    def __getattr__(self, attr):
        return getattr(self._settings, attr)

    def __delattr__(self, attr):
        self._settings.__delattr__(attr)

    def __getitem__(self, key):
        return getattr(self._settings, key)

    def __setitem__(self, key, value):
        setattr(self._settings, key, value)

    def __dir__(self):
        temp = self._settings.visible_traits()
        temp = [i for i in temp if not i.startswith('_')]
        return temp

    def __contains__(self, attr):
        return attr in self._settings.visible_traits()

    def _update(self, settings, docs=False):
        r"""
        Merges new key-value collections onto the existing object.

        This method allows several collections to be combined. It is hidden
        to avoid polluting the name space, which is reserved for actual data
        attributes.

        Parameters
        ----------
        settings : key-value collection
            Can be either a `SettingsData`` object which is a subclass
            of a the ``HasTraits`` class from the ``traits`` package, or a
            ``dict`` or ``dataclasss`` object.

        Notes
        -----
        In the case of ``dict`` or ``dataclass``-like objects, the data-types
        for the attributes are inferred from the values recieved. In all cases,
        the docs of the ``settings`` are ignored; this can only be set during
        init.

        """
        if hasattr(settings, 'items'): # Dictionary
            for k, v in settings.items():
                setattr(self, k, v)
        elif hasattr(settings, 'visible_traits'):
            for k in settings.visible_traits():
                v = getattr(settings, k)
                self._settings.add_trait(k, Trait(v, v.__class__))
        else:
            attrs = [i for i in dir(settings) if not i.startswith('_')]
            for k in attrs:
                v = getattr(settings, k)
                setattr(self, k, v)
        if docs:
            self._settings.__doc__ = settings.__doc__

    def __str__(self):
        return self._settings.__str__()

    def __repr__(self):
        return self.__str__()

    @property
    def __doc__(self):
        return self._settings.__doc__

    @property
    def _attrs(self):
        return self.__dir__()
