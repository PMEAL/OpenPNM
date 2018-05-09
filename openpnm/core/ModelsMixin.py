import inspect
import networkx as nx
from openpnm.core import Workspace, logging
from openpnm.utils.misc import PrintableDict
ws = Workspace()
logger = logging.getLogger()


class ModelsDict(PrintableDict):

    def dependency_list(self):
        r'''
        Returns a list of dependencies in the order with which they should be
        called to ensure data is calculated in the correct order.

        Notes
        -----
        This raises an exception if the graph has cycles which means the
        dependencies are unresolvable (i.e. there is no order which the
        models can be called that will work).  In this case it is possible
        to visually inspect the graph using ``dependency_graph``.

        See Also
        --------
        dependency_graph
        dependency_map

        '''
        dtree = self.dependency_graph()
        cycles = list(nx.simple_cycles(dtree))
        if cycles:
            raise Exception('Cyclic dependency found: ' + ' -> '.join(
                            cycles[0] + [cycles[0][0]]))
        return list(nx.algorithms.topological_sort(dtree))

    def dependency_graph(self):
        r"""
        Returns a NetworkX graph object of the dependencies

        See Also
        --------
        dependency_list
        dependency_map

        Notes
        -----
        To visualize the dependencies, the following NetworkX function and
        settings is helpful:

        nx.draw_spectral(d, arrowsize=50, font_size=32, with_labels=True,
                         node_size=2000, width=3.0, edge_color='lightgrey',
                         font_weight='bold')

        """
        dtree = nx.DiGraph()
        for propname in self.keys():
            dtree.add_node(propname)
            for dependency in self[propname].values():
                if dependency in list(self.keys()):
                    dtree.add_edge(dependency, propname)
        return dtree

    def dependency_map(self):
        r"""
        Create a graph of the dependency graph in a decent format

        See Also
        --------
        dependency_graph
        dependency_list

        """
        dtree = self.dependency_graph()
        fig = nx.draw_spectral(dtree,
                               with_labels=True,
                               arrowsize=50,
                               node_size=2000,
                               edge_color='lightgrey',
                               width=3.0,
                               font_size=32,
                               font_weight='bold')
        return fig


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

            *'normal'* : The model is run directly upon being assiged, and
            also run every time ``regenerate_models`` is called.

            *'constant'* : The model is run directly upon being assigned, but
            is not called again, thus making it's data act like a constant.
            If, however, the data is deleted from the object it will be
            regenerated again.

            *'deferred'* Is not run upon being assigned, but is run the first
            time that ``regenerate_models`` is called.

        """
        if propname in kwargs.values():  # Prevent infinite loops of look-ups
            raise Exception(propname+' can\'t be both dependency and propname')
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
        if regen_mode not in ['deferred']:
            self._regen(propname)

    def regenerate_models(self, propnames=None, exclude=[]):
        r"""
        Re-runs the specified model or models.

        Parameters
        ----------
        propnames : string or list of strings
            The list of property names to be regenerated.  If None are given
            then ALL models are re-run (except for those whose ``regen_mode``
            is 'constant').

        exclude : list of strings
            Since the default behavior is to run ALL models, this can be used
            to exclude specific models.  It may be more convenient to supply
            as list of 2 models to exclude than to specify 8 models to include.

        """
        if type(propnames) is str:  # Convert string to list if necessary
            propnames = [propnames]
        if propnames is None:  # If no props given, then regenerate them all
            propnames = self.models.dependency_list()
            # If some props are to be excluded, remove them from list
            if len(exclude) > 0:
                propnames = [i for i in propnames if i not in exclude]
        else:
            # Re-order given propnames according to dependency list
            all_props = self.models.dependency_list()
            propnames = [i for i in all_props if i in propnames]
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
        if self.settings['freeze_models']:
            # Don't run ANY models if freeze_models is set to True
            pass
        elif regen_mode == 'constant':
            # Only regenerate if data not already in dictionary
            if prop not in self.keys():
                self[prop] = model(target=self, **kwargs)
        else:
            try:
                self[prop] = model(target=self, **kwargs)
                logger.debug('Regnerating model: ' + prop)
            except KeyError:
                logger.warning(prop+' was not run due to missing dependencies')
                self.models[prop]['regen_mode'] = 'deferred'

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

    def _get_models(self):
        if not hasattr(self, '_models_dict'):
            self._models_dict = ModelsDict()
        return self._models_dict

    def _set_models(self, dict_):
        self._models_dict = ModelsDict()
        # Renerate all models in new dict if regen mode says so
        for model in dict_.keys():
            self.add_model(propname=model, **dict_[model])

    models = property(fget=_get_models, fset=_set_models)
