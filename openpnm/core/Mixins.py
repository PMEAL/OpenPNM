from openpnm.utils import PrintableDict


class ParamMixin:

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._params = PrintableDict()
        self._params._key = "Parameters"
        self._params._value = "Values"

    def __getitem__(self, key):
        if key.startswith('param'):
            try:
                vals = self._params[key]
            except KeyError:
                vals = self.network._params[key]
        else:
            vals = super().__getitem__(key)
        return vals

    def __setitem__(self, key, value):
        if key.startswith('param'):
            self._params[key] = value
        else:
            super().__setitem__(key, value)

    def __str__(self):
        s = super().__str__()
        s = s.rpartition('\n')[0]
        s = s + '\n' + self._params.__str__()
        return s

    def params(self):
        r"""
        Return parameter names and values in a dictionary
        """
        return self._params
