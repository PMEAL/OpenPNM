import inspect
from openpnm.core import Workspace
from openpnm.utils.misc import PrintableDict
ws = Workspace()


class ModelsDict(PrintableDict):

    def dependency_tree(self):
        tree = []
        for propname in self.keys():
            if propname not in tree:
                tree.append(propname)
            kwargs = self[propname].copy()
            kwargs.pop('model')
            kwargs.pop('regen_mode', None)
            for dependency in kwargs.values():
                if dependency in list(self.keys()):
                    tree.insert(tree.index(propname), dependency)

        unique = []
        [unique.append(item) for item in tree if item not in unique]
        return unique

    def __str__(self):
        horizontal_rule = '―' * 78
        lines = [horizontal_rule]
        strg = '{0:<3s} {1:<25s} {2:<25s} {3}'
        lines.append(strg.format('#', 'Property Name', 'Parameter', 'Value'))
        lines.append(horizontal_rule)
        for i, item in enumerate(self.keys()):
            temp = self[item].copy()
            regen_mode = temp.pop('regen_mode', None)
            model = str(temp.pop('model')).split(' ')[1]
            lines.append(strg.format(str(i+1), item, 'model:', model))
            for param in temp.keys():
                lines.append(strg.format('', '', param+':', temp[param]))
            lines.append(strg.format('', '', 'regeneration mode:', regen_mode))
            lines.append(horizontal_rule)
        return '\n'.join(lines)


class ModelsMixin():

    def add_model(self, propname, model, regen_mode='normal', **kwargs):
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
        # Regenerate model values if necessary
        if regen_mode != 'deferred':
            self._regen(propname)

    def regenerate_models(self, propnames=None, exclude=[]):
        # If no props given, then regenerate them all
        if propnames is None:
            propnames = list(self.models.dependency_tree())
        # If only one prop given, as string, put into a list
        elif type(propnames) is str:
            propnames = [propnames]
        [propnames.remove(i) for i in exclude if i in propnames]
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
        else:
            self[prop] = model(target=self, **kwargs)

    def _get_models(self):
        if not hasattr(self, '_models_dict'):
            self._models_dict = ModelsDict()
        return self._models_dict

    def _set_models(self, dict_):
        self._models_dict = ModelsDict()
        self._models_dict.update(dict_)
        # Renerate all models in new dict if regen mode says so
        for model in dict_:
            # In case regen mode is not set, do it now
            dict_[model].setdefault('regen_mode', 'normal')
        if self.settings['freeze_models']:
            pass
        else:
            self.regenerate_models()

    models = property(fget=_get_models, fset=_set_models)
