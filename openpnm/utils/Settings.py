from traits.api import HasTraits, Trait
from openpnm.utils import PrintableDict


class SettingsData(HasTraits):

    def __str__(self):
        d = PrintableDict()
        d._value = 'Settings'
        for item in self.visible_traits():
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

    def __getitem__(self, key):
        return getattr(self._settings, key)

    def __setitem__(self, key, value):
        setattr(self._settings, key, value)

    def __dir__(self):
        temp = self._settings.visible_traits()
        temp = [i for i in temp if not i.startswith('_')]
        return temp

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
                v = getattr(settings, k)
                self._settings.add_trait(k, Trait(v, v.__class__))
        else:
            attrs = [i for i in dir(settings) if not i.startswith('_')]
            for k in attrs:
                v = getattr(settings, k)
                self._settings.add_trait(k, Trait(v, v.__class__))

    def __str__(self):
        return self._settings.__str__()

    def __repr__(self):
        return self.__str__()

    @property
    def __doc__(self):
        return self._settings.__doc__


if __name__ == "__main__":

    from traits.api import Int, List, Str, Float, TraitError, ListStr

    # %% Standard initialization
    class S1(SettingsData):
        r"""
        This is a docstring
        """
        a = Int(1)
        b = Int(2)
        d = List(Str())

    sets1 = SettingsAttr(S1())
    print(sets1)
    assert "This is a docstring" in sets1.__doc__

    # %% Dataclass-style and testing inferred type is enforced
    class S2:
        r"""
        This is a docstring
        """
        a = 2
        b = 3

    sets2 = SettingsAttr(S2())
    print(sets2)
    assert "This is a docstring" in sets2.__doc__

    # %% Inheritance and testing immutability

    class S3(S2):
        r"""
        Different docstring
        """
        b = 3
        c = 4

    sets3 = SettingsAttr(S3())
    sets3.b = 44
    print(sets3)
    assert dir(sets3) == ['a', 'b', 'c']
    assert sets2.b != sets3.b
    assert "Different docstring" in sets3.__doc__

    try:
        sets3.b = 'nope'
    except TraitError as e:
        print(e)

    # %% Adding new attr and testing type is enforced
    class S4:
        r"""
        This is a docstring
        """
        a = 2
        b = 3.5

    sets4 = SettingsAttr(S4())
    sets4.c = Float(1.5)
    print(sets4)
    try:
        sets4.c = "string"
    except TraitError as e:
        print(e)

    # %% From dict
    d = {'a': 3, 'b': 4.5, 'c': []}

    sets5 = SettingsAttr(d)
    try:
        sets5.c = "string"
    except TraitError as e:
        print(e)
    sets5.c.append(1)
    print(sets5)
