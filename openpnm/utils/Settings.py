from traits.api import HasTraits
from openpnm.utils import Docorator, PrintableDict
import attr

docstr = Docorator()


class SettingsData(HasTraits):

    def __init__(self, settings={}):
        for k, v in settings.items():
            setattr(self, k.replace(' ', '_'), v)

    def __dir__(self):
        return self._attrs()

    def _attrs(self):
        temp = list(self.__base_traits__.keys())
        temp.extend(list(self.__dict__.keys()))
        temp = [i for i in temp if not i.startswith('_')]
        temp = [i for i in temp if not i.startswith('trait')]
        return temp

    def __str__(self):
        d = PrintableDict()
        d._value = 'Settings'
        for item in self._attrs():
            d[item] = getattr(self, item)
        return d.__str__()

    def __repr__(self):
        return self.__str__()


class SettingsAttr:

    _settings = None

    def __init__(self, settings=None):
        if settings is not None:
            super().__setattr__('_settings', settings)
        else:
            super().__setattr__('_settings', SettingsData())

    def _update(self, settings):
        for k, v in settings.items():
            setattr(self, k, v)

    def __setattr__(self, attr, val):
        setattr(self._settings, attr, val)

    def __getattr__(self, attr):
        return getattr(self._settings, attr)

    def __dir__(self):
        a = self._settings
        b = HasTraits()
        return sorted(list(set(dir(a)).difference(set(dir(b)))))

    def __str__(self):
        from openpnm.utils import PrintableDict
        d = PrintableDict()
        d._key = "Setting"
        for item in self.__dir__():
            d[item] = getattr(self, item)
        return d.__str__()

    @property
    def __doc__(self):
        return self._settings.__doc__


@attr.s
class SettingsAttr2:
    x = attr.ib(validator=attr.validators.instance_of(int),
                on_setattr=attr.setters.validate)
