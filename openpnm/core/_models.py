import logging
import openpnm as op
from openpnm.utils import PrintableDict, Workspace
from openpnm.utils import is_valid_propname
logger = logging.getLogger(__name__)
ws = Workspace()


__all__ = ['ModelsDict', 'ModelWrapper']


class ModelsDict(PrintableDict):
    """
    This subclassed dictionary is assigned to the ``models`` attribute of
    all objects that inherit from the ``ModelsMixin`` class.  Each
    dictionary entry corresponds to an entry in the target object's
    dictionary, and contains the models and associated parameters for
    generating the model.

    The main features of this subclass are three methods the help resolve
    the order in which models should be called: ``dependency_list``,
    ``dependency_graph``, and ``dependency_map``.

    """

    def _find_target(self):
        """
        Finds and returns the target object to which this ModelsDict is
        associated.
        """
        for proj in ws.values():
            for obj in proj:
                if hasattr(obj, "models"):
                    if obj.models is self:
                        return obj
        raise Exception("No target object found!")

    def dependency_list(self):
        r"""
        Returns a list of dependencies in the order with which they should
        be called to ensure data is calculated by one model before it's
        asked for by another.

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

        """
        import networkx as nx

        dtree = self.dependency_graph()
        cycles = list(nx.simple_cycles(dtree))
        if cycles:
            msg = 'Cyclic dependency: ' + ' -> '.join(cycles[0] + [cycles[0][0]])
            raise Exception(msg)
        d = nx.algorithms.dag.lexicographical_topological_sort(dtree, sorted)
        return list(d)

    def dependency_graph(self, deep=False):
        """
        Returns a NetworkX graph object of the dependencies

        Parameters
        ----------
        deep : bool, optional
            Defines whether intra- or inter-object dependency graph is
            desired. Default is False, i.e. only returns dependencies
            within the object.

        See Also
        --------
        dependency_list
        dependency_map

        """
        import networkx as nx

        dtree = nx.DiGraph()
        models = list(self.keys())

        for model in models:
            propname = model.split("@")[0]
            dtree.add_node(propname)
            # Filter pore/throat props only
            args = op.utils.flat_list2(self[model].values())
            dependencies = [arg for arg in args if is_valid_propname(arg)]
            # Add dependency from model's parameters
            for d in dependencies:
                dtree.add_edge(d, propname)

        return dtree

    def dependency_map(self,
                       ax=None,
                       figsize=None,
                       deep=False,
                       style='shell'):  # pragma: no cover
        """
        Create a graph of the dependency graph in a decent format

        Parameters
        ----------
        ax : matplotlib.axis, optional
            Matplotlib axis object on which dependency map is to be drawn.
        figsize : tuple, optional
            Tuple containing frame size.

        See Also
        --------
        dependency_graph
        dependency_list

        """
        import networkx as nx
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots()
        if figsize is not None:
            fig.set_size_inches(figsize)

        labels = {}
        node_shapes = {}
        dtree = self.dependency_graph(deep=deep)
        for node in dtree.nodes:
            labels[node] = node.split(".")[1]
            node_shapes[node] = "o" if node.startswith("pore") else "s"
        nx.set_node_attributes(dtree, node_shapes, "node_shape")

        layout = getattr(nx, f"{style}_layout")
        pos = layout(dtree)

        Pprops = [prop for prop in dtree.nodes if prop.startswith("pore")]
        Tprops = [prop for prop in dtree.nodes if prop.startswith("throat")]
        colors = ["yellowgreen", "coral"]
        shapes = ["o", "s"]

        for props, color, shape in zip([Pprops, Tprops], colors, shapes):
            nx.draw(
                dtree,
                pos=pos,
                nodelist=props,
                node_shape=shape,
                labels=labels,
                with_labels=True,
                edge_color='lightgrey',
                node_color=color,
                font_size=12,
                width=2.0
            )

        ax = plt.gca()
        ax.margins(x=0.2, y=0.05)

        return ax

    @property
    def _info(self):
        r"""
        Prints a nicely formatted list of model names and the domain to which
        they apply.

        Notes
        -----
        This is a hidden function for now, but could be exposed if useful.
        """
        names = {}
        for item in self:
            name, domain = item.split('@')
            if name not in names.keys():
                names[name] = []
            names[name].append(domain)
        D = PrintableDict(names, key='Model', value='Doamin')
        return D

    def __str__(self):
        horizontal_rule = '―' * 85
        lines = [horizontal_rule]
        strg = '{0:<3s} {1:<35s} {2:<25s} {3}'
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

    def __delitem__(self, key):
        if '@' in key:
            super().__delitem__(key)
        else:  # Delete all models with the same prefix
            for item in list(self.keys()):
                if item.startswith(key):
                    self.__delitem__(item)

    def __getitem__(self, key):
        try:
            return super().__getitem__(key)
        except KeyError:
            d = PrintableDict(key='Model', value='Args')
            for k, v in self.items():
                if k.startswith(key):
                    d[k] = v
            if len(d) > 0:
                return d
            else:
                raise KeyError(key)

    def update(self, d, domain='all'):
        # Catch un-run function
        if hasattr(d, '__call__'):
            raise Exception('Received dict argument is a function, try running it')
        parent = self._find_target()
        for k, v in d.items():
            parent.add_model(propname=k, domain=domain, **v)


class ModelWrapper(dict):
    """
    This class is used to hold individual models and provide some extra
    functionality, such as pretty-printing and the ability to run itself.
    """

    def __call__(self):
        target = self._find_target()
        model = self['model']
        kwargs = {}
        for k, v in self.items():
            if k not in ['model', 'regen_mode']:
                kwargs[k] = v
        return model(target=target, **kwargs)

    def _find_parent(self):
        r"""
        Finds and returns the parent object to self.
        """
        for proj in ws.values():
            for obj in proj:
                if hasattr(obj, "models"):
                    for mod in obj.models.keys():
                        if obj.models[mod] is self:
                            return obj
        raise Exception("No parent object found!")

    @property
    def name(self):
        for proj in ws.values():
            for obj in proj:
                if hasattr(obj, 'models'):
                    for key, mod in obj.models.items():
                        if mod is self:
                            return key

    @property
    def propname(self):
        return self.name.split('@')[0]

    @property
    def domain(self):
        element, prop = self.name.split('.', 1)
        prop, domain = prop.split('@')
        return element + '.' + domain

    def __str__(self):
        horizontal_rule = '―' * 78
        lines = [horizontal_rule]
        strg = '{0:<25s} {1:<25s} {2}'
        lines.append(strg.format('Property Name', 'Parameter', 'Value'))
        lines.append(horizontal_rule)
        temp = self.copy()
        regen_mode = temp.pop('regen_mode', None)
        model = str(temp.pop('model')).split(' ')[1]
        lines.append(strg.format(self.propname, 'model:', model))
        for param in temp.keys():
            lines.append(strg.format('', param+':', temp[param]))
        lines.append(strg.format('', 'regeneration mode:', regen_mode))
        lines.append(horizontal_rule)
        return '\n'.join(lines)

    def _find_target(self):
        """
        Finds and returns the parent object to self.
        """
        for proj in ws.values():
            for obj in proj:
                if hasattr(obj, "models"):
                    for mod in obj.models.values():
                        if mod is self:
                            return obj
        raise Exception("No target object found!")


# class ModelsMixin:
#     """
#     This class is meant to be combined by the Base class in multiple
#     inheritence. This approach is used since Network and Algorithm do not
#     need to have any ``models`` attribute, while Phase, Geometry, and
#     Physics do. By using a mixin class, all objects can inherit from Base
#     while the model functionality can be added only where needed.
#     """

#     def add_model(self, propname, model, regen_mode='', **kwargs):
#         """
#         Adds a new model to the models dictionary.

#         Parameters
#         ----------
#         propname : str
#             The name of the property to be calculated by the model.
#         model : function
#             A reference (handle) to the function to be used.
#         regen_mode : str
#             Controls how/when the model is run (See Notes for more details).
#             Options are:

#             ===========  =====================================================
#             mode         meaning
#             ===========  =====================================================
#             'normal'     The model is run directly upon being
#                          assigned, and also run every time ``regenerate_models``
#                          is called.
#             'constant'   The model is run directly upon being assigned, but
#                          is not called again, thus making its data act like a
#                          constant. If, however, the data is deleted from the
#                          object it will be regenerated again.
#             'deferred'   Is not run upon being assigned, but is run the first
#                          time that ``regenerate_models`` is called.
#             'explicit'   Is only run if the model name is explicitly passed
#                          to the ``regenerate_models`` method.  This allows
#                          full control of when the model is run.
#             ===========  =====================================================


#         """
#         if propname in kwargs.values():  # Prevent infinite loops of look-ups
#             raise Exception(propname+' can\'t be both dependency and propname')
#         # Look for default regen_mode in settings if present, else use 'normal'
#         if regen_mode == '':
#             if 'regen_mode' in self.settings._attrs:
#                 regen_mode = self.settings['regen_mode']
#             else:
#                 regen_mode = 'normal'
#         # Add model and regen_mode to kwargs dictionary
#         kwargs.update({'model': model, 'regen_mode': regen_mode})
#         # Insepct model to extract arguments and default values
#         if model.__defaults__:
#             vals = list(inspect.getfullargspec(model).defaults)
#             keys = inspect.getfullargspec(model).args[-len(vals):]
#             for k, v in zip(keys, vals):  # Put defaults into kwargs
#                 if k not in kwargs:  # Skip if argument was given in kwargs
#                     kwargs.update({k: v})
#         self.models[propname] = ModelWrapper(kwargs)  # Store all kwargs
#         # Regenerate model values if necessary
#         if regen_mode not in ['deferred', 'explicit']:
#             self._regen(propname)

#     def regenerate_models(self, propnames=None, exclude=[], deep=False):
#         """
#         Re-runs the specified model or models.

#         Parameters
#         ----------
#         propnames : str or list of str
#             The list of property names to be regenerated.  If none are given
#             then ALL models are re-run (except for those whose ``regen_mode``
#             is 'constant').
#         exclude : list of str
#             Since the default behavior is to run ALL models, this can be used
#             to exclude specific models.  It may be more convenient to supply
#             as list of 2 models to exclude than to specify 8 models to include.

#         """
#         # If empty list of propnames was given, do nothing and return
#         if isinstance(propnames, list) and len(propnames) == 0:
#             return
#         if isinstance(propnames, str):  # Convert string to list if necessary
#             propnames = [propnames]
#         if propnames is None:  # If no props given, then regenerate them all
#             propnames = self.models.dependency_list()
#             # If some props are to be excluded, remove them from list
#             for k, v in self.models.items():
#                 if 'regen_mode' not in v:
#                     pass
#                 elif v['regen_mode'] == 'explicit':
#                     exclude.extend([k])
#             propnames = [i for i in propnames if i not in exclude]
#         # Re-order given propnames according to dependency tree
#         self_models = self.models.dependency_list()
#         propnames = [i for i in self_models if i in propnames]
#         for item in propnames:
#             self._regen(item)

#     def _regen(self, prop):
#         # Create a temporary dict of all model arguments
#         model = self.models[prop]['model']
#         regen_mode = self.models[prop]['regen_mode']
#         kwargs = {}
#         for k, v in self.models[prop].items():
#             if k not in ['model', 'regen_mode']:
#                 kwargs[k] = v
#         # Only regenerate model if regen_mode is correct
#         if regen_mode == 'constant':
#             # Only regenerate if data not already in dictionary
#             if prop not in self.keys():
#                 self[prop] = model(target=self, **kwargs)
#         else:
#             try:
#                 self[prop] = model(target=self, **kwargs)
#             except KeyError as e:
#                 msg = (f"{prop} was not run since the following property"
#                        f" is missing: {e}")
#                 logger.error(prettify_logger_message(msg))
#                 self.models[prop]['regen_mode'] = 'deferred'

#     def remove_model(self, propname=None, mode=['model', 'data']):
#         """
#         Removes model and data from object.

#         Parameters
#         ----------
#         propname : str or list[str]
#             The property or list of properties to remove
#         mode : list[str]
#             Controls what is removed. Options are:

#             ===========  =====================================================
#             mode         meaning
#             ===========  =====================================================
#             'model'      Removes the model but not any numerical data that may
#                          already exist.
#             data'        Removes the data but leaves the model.
#             ===========  =====================================================
#         The default is both.

#         """
#         if isinstance(propname, str):
#             propname = [propname]
#         for item in propname:
#             if 'model' in mode:
#                 if item in self.models.keys():
#                     del self.models[item]
#             if 'data' in mode:
#                 if item in self.keys():
#                     del self[item]

#     def _get_models(self):
#         """List of available models on the objects"""
#         if not hasattr(self, '_models_dict'):
#             self._models_dict = ModelsDict()
#         return self._models_dict

#     def _set_models(self, dict_):
#         self._models_dict = ModelsDict()
#         # Renerate all models in new dict if regen mode says so
#         for model in dict_.keys():
#             self.add_model(propname=model, **dict_[model])

#     models = property(fget=_get_models, fset=_set_models)
