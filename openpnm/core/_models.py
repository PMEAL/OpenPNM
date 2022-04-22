import logging
import inspect
import numpy as np
from openpnm.utils import PrintableDict, Workspace
from openpnm.utils import is_valid_propname
from openpnm.utils import prettify_logger_message
logger = logging.getLogger(__name__)
ws = Workspace()

__all__ = ['ModelsDict', 'ModelsMixin', 'ModelWrapper']


class ModelsDict(PrintableDict):
    r"""
    This subclassed dictionary is assigned to the ``models`` attribute of
    all objects that inherit from the ``ModelsMixin`` class.  Each
    dictionary entry corresponds to an entry in the target object's
    dictionary, and contains the models and associated parameters for
    generating the model.

    The main features of this subclass are three methods the help resolve
    the order in which models should be called: ``dependency_list``,
    ``dependency_graph``, and ``dependency_map``.

    """

    def update(self, models):
        target = self._find_parent()
        for item in models.keys():
            target.add_model(propname=item, **models[item])

    def _find_parent(self):
        r"""
        Finds and returns the parent object to self.
        """
        for proj in ws.values():
            for obj in proj:
                if hasattr(obj, "models"):
                    if obj.models is self:
                        return obj
        raise Exception("No parent object found!")

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
            raise Exception('Cyclic dependency found: ' + ' -> '.join(
                            cycles[0] + [cycles[0][0]]))
        d = nx.algorithms.dag.lexicographical_topological_sort(dtree, sorted)
        return list(d)

    def dependency_graph(self, deep=False):
        r"""
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

        Notes
        -----
        To visualize the dependencies, the following NetworkX function and
        settings is helpful:

        .. plot::

           import networkx as nx
           import openpnm as op
           import matplotlib.pyplot as plt

           net = op.network.Cubic(shape=[3, 3, 3])
           geo = op.geometry.SpheresAndCylinders(network=net,
                                                 pores=net.Ps,
                                                 throats=net.Ts)

           dtree = geo.models.dependency_graph()
           nx.draw_spectral(dtree,
                            arrowsize=50,
                            font_size=32,
                            with_labels=True,
                            node_size=2000,
                            width=3.0,
                            edge_color='lightgrey',
                            font_weight='bold')

           plt.show()

        """
        import networkx as nx

        dtree = nx.DiGraph()
        models = list(self.keys())
        # Fetch model-less props: those w/o any model, like temperature
        # otherwise, they won't get picked up in the dependency graph.
        all_props = list(self._find_parent().keys())
        exclude_keys = ["pore.all", "throat.all"]
        pure_props = np.setdiff1d(all_props, models + exclude_keys).tolist()

        for model in models:
            dtree.add_node(model)
            # Filter pore/throat props only
            dependencies = set()
            for param in self[model].values():
                if type(param) == list:
                    for element in param:
                        if is_valid_propname(element):
                            dependencies.add(element)
                else:
                    if is_valid_propname(param):
                        dependencies.add(param)
            # Add depenency from model's parameters
            for d in dependencies:
                if not deep:
                    if d in models + pure_props:
                        dtree.add_edge(d, model)
                else:
                    dtree.add_edge(d, model)

        return dtree

    def dependency_map(self,
                       ax=None,
                       figsize=None,
                       deep=False,
                       style='shell'):  # pragma: no cover
        r"""
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
        dtree = self.dependency_graph(deep=deep)
        for node in dtree.nodes:
            labels[node] = node.split(".")[1]

        node_shapes = {}
        for node in dtree.nodes:
            if node.startswith("pore"):
                node_shapes[node] = "o"
            else:
                node_shapes[node] = "s"
        nx.set_node_attributes(dtree, node_shapes, "node_shape")

        style = style.replace('draw_', '')
        draw = getattr(nx, 'draw_' + style)

        pore_props = [prop for prop in dtree.nodes if prop.startswith("pore")]
        throat_props = [prop for prop in dtree.nodes if prop.startswith("throat")]

        for props, color, shape in zip([pore_props, throat_props],
                                       ["yellowgreen", "coral"],
                                       ["o", "s"]):
            draw(dtree,
                 nodelist=props,
                 node_shape=shape,
                 labels=labels,
                 with_labels=True,
                 edge_color='lightgrey',
                 node_color=color,
                 font_size=12,
                 width=2.0)

        ax = plt.gca()
        ax.margins(x=0.2, y=0.02)

        return ax

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


class ModelWrapper(dict):
    r"""
    This class is used to hold individual models and provide some extra
    functionality, such as pretty-printing.
    """
    @property
    def propname(self):
        for proj in ws.values():
            for obj in proj:
                if hasattr(obj, 'models'):
                    for key, mod in obj.models.items():
                        if mod is self:
                            return key

    def __repr__(self):
        return self.__str__()

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


class ModelsMixin:
    """
    This class is meant to be combined by the Base class in multiple
    inheritence. This approach is used since Network and Algorithm do not
    need to have any ``models`` attribute, while Phase, Geometry, and
    Physics do. By using a mixin class, all objects can inherit from Base
    while the model functionality can be added only where needed.
    """

    def add_model(self, propname, model, regen_mode='', **kwargs):
        r"""
        Adds a new model to the models dictionary.

        Parameters
        ----------
        propname : str
            The name of the property to be calculated by the model.
        model : function
            A reference (handle) to the function to be used.
        regen_mode : str
            Controls how/when the model is run (See Notes for more details).
            Options are:

            ===========  =====================================================
            mode         meaning
            ===========  =====================================================
            'normal'     The model is run directly upon being
                         assigned, and also run every time ``regenerate_models``
                         is called.
            'constant'   The model is run directly upon being assigned, but
                         is not called again, thus making its data act like a
                         constant. If, however, the data is deleted from the
                         object it will be regenerated again.
            'deferred'   Is not run upon being assigned, but is run the first
                         time that ``regenerate_models`` is called.
            'explicit'   Is only run if the model name is explicitly passed
                         to the ``regenerate_models`` method.  This allows
                         full control of when the model is run.
            ===========  =====================================================


        """
        if propname in kwargs.values():  # Prevent infinite loops of look-ups
            raise Exception(propname+' can\'t be both dependency and propname')
        # Look for default regen_mode in settings if present, else use 'normal'
        if regen_mode == '':
            if 'regen_mode' in self.settings._attrs:
                regen_mode = self.settings['regen_mode']
            else:
                regen_mode = 'normal'
        # Add model and regen_mode to kwargs dictionary
        kwargs.update({'model': model, 'regen_mode': regen_mode})
        # Insepct model to extract arguments and default values
        if model.__defaults__:
            vals = list(inspect.getfullargspec(model).defaults)
            keys = inspect.getfullargspec(model).args[-len(vals):]
            for k, v in zip(keys, vals):  # Put defaults into kwargs
                if k not in kwargs:  # Skip if argument was given in kwargs
                    kwargs.update({k: v})
        self.models[propname] = ModelWrapper(kwargs)  # Store all kwargs
        # Regenerate model values if necessary
        if regen_mode not in ['deferred', 'explicit']:
            self._regen(propname)

    def regenerate_models(self, propnames=None, exclude=[], deep=False):
        r"""
        Re-runs the specified model or models.

        Parameters
        ----------
        propnames : str or list of str
            The list of property names to be regenerated.  If none are given
            then ALL models are re-run (except for those whose ``regen_mode``
            is 'constant').
        exclude : list of str
            Since the default behavior is to run ALL models, this can be used
            to exclude specific models.  It may be more convenient to supply
            as list of 2 models to exclude than to specify 8 models to include.

        """
        # If empty list of propnames was given, do nothing and return
        if isinstance(propnames, list) and len(propnames) == 0:
            return
        if isinstance(propnames, str):  # Convert string to list if necessary
            propnames = [propnames]
        if propnames is None:  # If no props given, then regenerate them all
            propnames = self.models.dependency_list()
            # If some props are to be excluded, remove them from list
            for k, v in self.models.items():
                if 'regen_mode' not in v:
                    pass
                elif v['regen_mode'] == 'explicit':
                    exclude.extend([k])
            propnames = [i for i in propnames if i not in exclude]
        # Re-order given propnames according to dependency tree
        self_models = self.models.dependency_list()
        propnames = [i for i in self_models if i in propnames]
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
            # Only regenerate if data not already in dictionary
            if prop not in self.keys():
                self[prop] = model(target=self, **kwargs)
        else:
            try:
                self[prop] = model(target=self, **kwargs)
            except KeyError as e:
                msg = (f"{prop} was not run since the following property"
                       f" is missing: {e}")
                logger.error(prettify_logger_message(msg))
                self.models[prop]['regen_mode'] = 'deferred'

    def remove_model(self, propname=None, mode=['model', 'data']):
        r"""
        Removes model and data from object.

        Parameters
        ----------
        propname : str or list[str]
            The property or list of properties to remove
        mode : list[str]
            Controls what is removed. Options are:

            ===========  =====================================================
            mode         meaning
            ===========  =====================================================
            'model'      Removes the model but not any numerical data that may
                         already exist.
            data'        Removes the data but leaves the model.
            ===========  =====================================================
        The default is both.

        """
        if isinstance(propname, str):
            propname = [propname]
        for item in propname:
            if 'model' in mode:
                if item in self.models.keys():
                    del self.models[item]
            if 'data' in mode:
                if item in self.keys():
                    del self[item]

    def _get_models(self):
        """List of available models on the objects"""
        if not hasattr(self, '_models_dict'):
            self._models_dict = ModelsDict()
        return self._models_dict

    def _set_models(self, dict_):
        self._models_dict = ModelsDict()
        # Renerate all models in new dict if regen mode says so
        for model in dict_.keys():
            self.add_model(propname=model, **dict_[model])

    models = property(fget=_get_models, fset=_set_models)
