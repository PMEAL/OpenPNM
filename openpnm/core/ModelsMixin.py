import inspect
import copy
import numpy as np
import matplotlib.pyplot as plt
from openpnm.core import Workspace
from openpnm.utils.misc import PrintableDict
ws = Workspace()


class ModelWrapper(PrintableDict):

    def run(self):
        target = self._find_self()
        vals = self['model'](target=target, **self['kwargs'])
        return vals

    def show_hist(self, fig=None):
        if fig is None:
            fig = plt.figure()
        plt.hist(self.run())
        return fig

    def _find_self(self):
        for sim in ws.values():
            for obj in sim:
                if self in obj.models.values():
                    return obj


class ModelsDict(dict):

    def dependency_tree(self):
        tree = []
        for propname in self.keys():
            if propname not in tree:
                tree.append(propname)
            for dependency in self[propname]['kwargs'].values():
                if dependency in list(self.keys()):
                    tree.insert(tree.index(propname), dependency)

        unique = []
        [unique.append(item) for item in tree if item not in unique]
        return unique

    def __str__(self):
        horizontal_rule = '―' * 60
        lines = [horizontal_rule]
        lines.append('{0:<5s} {1:<30s} {2}'.format('#',
                                                   'Property Name',
                                                   'Regeneration Mode'))
        lines.append(horizontal_rule)
        for i, item in enumerate(self.keys()):
            str = '{0:<5d} {1:<30s} {2:<20s}'
            lines.append(str.format(i + 1, item, self[item]['regen_mode']))
        lines.append(horizontal_rule)
        return '\n'.join(lines)


class ModelsMixin():

    def add_model(self, propname, model, regen_mode='deferred', **kwargs):
        self.models[propname] = ModelWrapper()
        self.models[propname]['model'] = model
        if regen_mode in ['normal', 'deferred', 'constant']:
            self.models[propname]['regen_mode'] = regen_mode
        else:
            raise Exception('Unexpected regeneration mode')

        # Insepct model and extract arguments and default values
        if model.__defaults__ is not None:
            vals = list(inspect.getargspec(model).defaults)
            keys = inspect.getargspec(model).args[-len(vals):]
            # Put defaults into the dict
            for k, v in zip(keys, vals):
                # Skip if argument was given in kwargs
                if k not in kwargs:
                    kwargs.update({k: v})

        self.models[propname]['kwargs'] = kwargs
        if regen_mode in ['normal', 'constant']:
            self._regen(propname)

    def regenerate_models(self, propnames=None):
        # If no props give, then regenerate them all
        if propnames is None:
            propnames = list(self.models.dependency_tree())
        # If only one prop give, as string, put into a list
        elif type(propnames) is str:
            propnames = [propnames]
        for item in propnames:
            if self.models[item]['regen_mode'] == 'constant':
                if item not in self.keys():
                    self._regen(item)
            else:
                    self._regen(item)

    def _regen(self, propname):
        f = self.models[propname]['model']
        values = f(target=self, **self.models[propname]['kwargs'])
        self[propname] = values

    # The use of a property attribute here is because I can't just set
    # self.models= {} in the init, since the damn init won't run!
    def _get_models(self):
        if not hasattr(self, '_dict'):
            self._dict = ModelsDict()
        return self._dict

    def _set_models(self, _dict):
        self._dict = copy.deepcopy(_dict)

    models = property(fget=_get_models, fset=_set_models)
