from traits.api import HasTraits
from openpnm.utils import Docorator

docstr = Docorator()


class SettingsData(HasTraits):
    r"""

    """


class SettingsAttr:

    def __init__(self, settings={}, defaults=None):
        if defaults is not None:
            super().__setattr__('_settings', defaults)
        else:
            super().__setattr__('_settings', SettingsData())
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


