import inspect
from openpnm.core import Workspace
from openpnm.utils.misc import PrintableDict
ws = Workspace()


class ModelsDict(PrintableDict):

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

    def add_model(self, propname, model, regen_mode='deferred', **kwargs):
        r"""
        Adds a new model to the models dictionary (``object.models``)

        Parameters
        ----------
        propname : string
            The name of the property to be calculated by the model.

        model : function
            A reference (handle) to the function to be used.

        regen_mode : string
            Controls how/when the model is run (See Notes for more details).
            Options are:

            *'normal'* or *'eager'* : The model is run directly upon being
            assiged, and also run every time ``object.regenerate_models`` is
            called.

            *'constant'* : The model is run directly upon being assigned, but
            is not called again, thus making it's data act like a constant.

            *'deferred'* or *'lazy'*: (default) Is not run upon being assigned,
            but is run the first time that it's data is requested.

        Notes
        -----
        The difference between 'eager' and 'lazy' execution is useful to
        understand.  In 'eager' mode the model is run as soon as it's assigned
        to the object, which means that users must be careful to assign
        dependent models first, or else a KeyError will be raised since needed
        data is not present.  In 'lazy' mode models are not run until their
        data is asked for, which creates a cascading call to all necessary
        models that have not been run yet.  The latter behavior is more
        convenient, but can be confusing since a lot happens behind the scenes.

        For example: In 'eager' mode, if 'pore.volume' is assigned before
        'pore.diameter' a KeyError will occur when the model attempts to
        access the 'pore.diameter' data.  In 'lazy' mode the request for
        `pore.volume` data will run the 'pore.volume' model, which attempts
        to access `pore.diameter` data, and upon failing to find it will
        attempt to run the `pore.diameter` model.

        """
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
        if regen_mode not in ['deferred', 'lazy']:
            self._regen(propname)

    def remove_model(self, propname=None, mode=['model', 'data']):
        r"""
        Removes model and data from object.

        Parameters
        ----------
        propname : string or list of strings
            The property or list of properties to remove

        mode : list of strings
            Controls what is removed.  Options are:

            *'model'* : Removes the model but not any numerical data that may
            already exist.

            *'data'* : Removes the data but leaves the model.

        The default is both.

        """
        if type(propname) is str:
            propname = [propname]
        for item in propname:
            if 'model' in mode:
                if item in self.models.keys():
                    del self.models[item]
            if 'data' in mode:
                if item in self.keys():
                    del self[item]

    def regenerate_models(self, propnames=None, exclude=[]):
        r"""
        Re-runs the specified models.

        Parameters
        ----------
        propnames : string or list of strings
            The list of property names to be regenerated.  If None are given
            then ALL models are re-run (except for those whose ``regen_mode``
            is 'constant').

        exclude : list of strings
            Since the default behavior is to run ALL models, this can be used
            to exclude specific models.  It may be more convenient to supply
            as list of 2 models to exclude than to specify 8 models include.
        """
        # If no props given, then regenerate them all
        if propnames is None:
            propnames = list(self.models.keys())
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
        regen_mode = kwargs.pop('regen_mode', None)
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
