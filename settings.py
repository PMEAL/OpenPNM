from traits.api import HasTraits, Trait, Int, Float
from openpnm.utils import PrintableDict


class SettingsTest:
    r"""
    """

    def __init__(self, settings):
        r"""
        Parameters
        ----------
        settings : key-value collection
            Ideally should be an `SettingsData`` object which is a subclass
            of a the ``HasTraits`` class from the ``traits`` package. The
            docstring for this class is adopted from received argument only
            **if** a ``HasTraits`` class is receieved.  Optionally it can
            be a ``dict`` or any ``dataclasss``-type object with data stored
            as attribute. In both cases the data-types for the attributes
            are inferred from the values received.

        Notes
        -----
        This class is a wrapper around the ``HasTraits`` class from the
        ``traits`` package. The main point is the provide a clean namespace
        on the object so that only the specified attributes show up.  It also
        allows for combining multiple objects into one using the ``_update``
        method.
        """
        if hasattr(settings, 'traits'):
            super().__setattr__('_settings', settings)
        else:
            super().__setattr__('_settings', HasTraits())
            self._update(settings)

    def __setattr__(self, attr, val):
        setattr(self._settings, attr, val)

    def __getattr__(self, attr):
        return getattr(self._settings, attr)

    def __dir__(self):
        return self._settings.visible_traits()

    def _update(self, settings):
        r"""
        Merges new key-value collections onto the existing object.

        This method allows several collections to be combined. It is hidden
        to avoid polluting the name space which is reserved for actual data
        attributes.

        Parameters
        ----------
        settings : key-value collection
            Can be either a `SettingsData`` object which is a subclass
            of a the ``HasTraits`` class from the ``traits`` package, or a
            ``dict`` or ``dataclasss`` object.

        Notes
        -----
        In the case of ``dict`` or ``dataclass``-like objects, the
        data-types for the attributes are inferred from the values
        recieved. In the case of a proper ``SettingsData`` object,
        the ``__doc__`` attribute of the received object is adopted
        by the wrapper class.

        """
        if hasattr(settings, 'visible_traits'):
            for k in settings.visible_traits():
                self._settings.add_trait(k, getattr(settings, k))
                self._settings.__doc__ = settings.__doc__
        elif hasattr(settings, 'items'):
            for k, v in settings.items():
                self._settings.add_trait(k, Trait(v, v.__class__))
        else:
            attrs = [i for i in dir(settings) if not i.startswith('_')]
            for k in attrs:
                self._settings.add_trait(k, getattr(settings, k))

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


class S1(HasTraits):
    r"""
    Blah blah blah
    """
    y = Int(3)
    z = Float(4)

s = SettingsTest(S1())

s._update({'x': 1,
            'y': [2, 3, 5]})

from dataclasses import dataclass

@dataclass
class S2:
    a : int = 55

s._update(S2)

class S3:
    b = 66

s._update(S3())

print(s)




















