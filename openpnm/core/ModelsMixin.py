import inspect
from openpnm.core import Workspace
from openpnm.utils.misc import PrintableDict
ws = Workspace()


class ModelsDict(dict):

    def dependency_tree(self):
        tree = []
        for propname in self.keys():
            if propname not in tree:
                tree.append(propname)
            args = self[propname].copy()
            args.pop('model')
            args.pop('regen_mode', None)
            for dependency in args:
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
        # Add model and regen_mode to kwargs dictionary
        kwargs.update({'model': model, 'regen_mode': regen_mode})
        # Insepct model to extract arguments and default values
        if model.__defaults__:
            vals = list(inspect.getargspec(model).defaults)
            keys = inspect.getargspec(model).args[-len(vals):]
            # Put defaults into dict_
            for k, v in zip(keys, vals):
                # Skip if argument was given in kwargs
                if k not in kwargs:
                    kwargs.update({k: v})
        # Store all keyword argumnents in model
        self.models[propname] = kwargs
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
        # Scan through list of propnames and regenerate each one
        for item in propnames:
            self._regen(item)

    def _regen(self, prop):
        # Create a temporary dict of all model arguments
        kwargs = self.models[prop].copy()
        # Pop model and regen_mode from temporary dict
        model = kwargs.pop('model')
        regen_mode = kwargs.pop('regen_mode', 'normal')
        # Only regenerate model if regen_mode is correct
        if regen_mode == 'constant':
            if prop not in self.keys():
                self[prop] = model(target=self, **kwargs)
        elif regen_mode in ['normal', 'deferred']:
            self[prop] = model(target=self, **kwargs)

    # The use of a property attribute here is because I can't just set
    # self.models= {} in the init, since the damn init won't run!
    def _get_models(self):
        if not hasattr(self, '_dict'):
            self._dict = ModelsDict()
        return self._dict

    def _set_models(self, dict_):
        self._dict = ModelsDict()
        self._dict.update(dict_)

    models = property(fget=_get_models, fset=_set_models)
