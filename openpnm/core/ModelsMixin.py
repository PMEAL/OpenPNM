import inspect
import copy
import matplotlib.pyplot as plt
from openpnm.core import Workspace
from openpnm.utils.misc import PrintableDict
ws = Workspace()


class ModelWrapper(PrintableDict):

    def __init__(self, dict_={}):
        super().__init__()
        # Insepct model to extract arguments and default values
        model = dict_.pop('model', None)
        if model and model.__defaults__ is not None:
            vals = list(inspect.getargspec(model).defaults)
            keys = inspect.getargspec(model).args[-len(vals):]
            # Put defaults into dict_
            for k, v in zip(keys, vals):
                # Skip if argument was given in kwargs
                if k not in dict_:
                    dict_.update({k: v})
        self['model'] = model
        # Set regeneration mode
        self['regen_mode'] = dict_.pop('regen_mode', 'deferred')
        # Put the remaining keyword arguments in 'kwargs'
        self['kwargs'] = dict_

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

    def copy(self):
        return copy.deepcopy(self)

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
        # Add model and regen_mode to kwargs dictionary
        kwargs.update({'model': model, 'regen_mode': regen_mode})
        # Send whole dictionary to ModelWrapper's init
        self.models[propname] = ModelWrapper(kwargs)
        # Regenerate model values is necessary
        if self.models[propname]['regen_mode'] in ['normal', 'constant']:
            self._regen(propname)

    def regenerate_models(self, propnames=None):
        # If no props given, then regenerate them all
        if propnames is None:
            propnames = list(self.models.dependency_tree())
        # If only one prop given, as string, put into a list
        elif type(propnames) is str:
            propnames = [propnames]
        for item in propnames:
            # Only regenerate constant items if not already present
            if self.models[item]['regen_mode'] == 'constant':
                if item not in self.keys():
                    self._regen(item)
            # Regenerate everything else
            else:
                self._regen(item)

    def _regen(self, propname):
        f = self.models[propname]['model']
        vals = f(target=self, **self.models[propname]['kwargs'])
        self[propname] = vals

    # The use of a property attribute here is because I can't just set
    # self.models= {} in the init, since the damn init won't run!
    def _get_models(self):
        if not hasattr(self, '_dict'):
            self._dict = ModelsDict()
        return self._dict

    def _set_models(self, _dict):
#        if not hasattr(self, '_dict'):
#            self._dict = ModelsDict()
        self._dict = ModelsDict(copy.deepcopy(_dict))

    models = property(fget=_get_models, fset=_set_models)
