from traits.api import HasTraits, Trait
from openpnm.utils import PrintableDict


class SettingsData(HasTraits):

    def __str__(self):
        d = PrintableDict()
        d._value = 'Settings'
        for item in self._attrs():
            d[item] = getattr(self, item)
        return d.__str__()

    def __repr__(self):
        return self.__str__()


class SettingsAttr:
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
            as attributes. In both cases the data-types for the attributes
            are inferred from the values received, but the docs are ignored.

        Notes
        -----
        This class is a wrapper around the ``HasTraits`` class from the
        ``traits`` package. The main point is the provide a clean namespace
        on the object so that only the specified attributes show up.  It also
        allows for combining multiple objects into one using the ``_update``
        method.
        """
        if hasattr(settings, 'traits'):
            # If HasTraits object, assign it to _settings
            super().__setattr__('_settings', settings)
        else:
            # Otherwise add an empty HasTraits object then fill it
            super().__setattr__('_settings', HasTraits())
            self._update(settings)

    def __setattr__(self, attr, val):
        setattr(self._settings, attr, val)

    def __getattr__(self, attr):
        return getattr(self._settings, attr)

    def __getitem__(self, key):
        return getattr(self._settings, key)

    def __setitem__(self, key, value):
        setattr(self._settings, key, value)
    def __dir__(self):
        return self._settings.visible_traits()

    def _update(self, settings):
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
                self._settings.add_trait(k, Trait(v, v.__class__))
        elif hasattr(settings, 'visible_traits'):
            for k in settings.visible_traits():
                self._settings.add_trait(k, getattr(settings, k))
        else:
            attrs = [i for i in dir(settings) if not i.startswith('_')]
            for k in attrs:
                v = getattr(settings, k)
                self._settings.add_trait(k, Trait(v, v.__class__))

    def __str__(self):
        d = PrintableDict()
        for item in self.__dir__():
            d[item] = getattr(self, item)
        return d.__str__()

    def __repr__(self):
        return self.__str__()

    @property
    def __doc__(self):
        return self._settings.__doc__
