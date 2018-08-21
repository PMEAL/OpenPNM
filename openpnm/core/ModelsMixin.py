import inspect
import networkx as nx
from openpnm.utils import PrintableDict, logging, Workspace
ws = Workspace()
logger = logging.getLogger(__name__)


class ModelsDict(PrintableDict):
    r"""
    This subclassed dictionary is assigned to the ``models`` attribute of
    all objects that inherit from the ``ModelsMixin`` class.  Each dictionary
    entry corresponds to an entry in the target object's dictionary, and
    contains the models and associated parameters for generating the model.

    The main features of this subclass are three methods the help resolve the
    order in which models should be called: ``dependency_list``,
    ``dependency_graph``, and ``dependency_map``.

    """

    def dependency_list(self):
        r'''
        Returns a list of dependencies in the order with which they should be
        called to ensure data is calculated by one model before it's asked for
        by another.

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
        d = nx.algorithms.dag.lexicographical_topological_sort(dtree, sorted)
        return list(d)

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
    r"""
    This class is meant to be combined by the Base class in multiple
    inheritence.  This approach is used since Network and Algorithm do not
    need to have any ``models`` attribute, while Phase, Geometry, and Physics
    do.  By using a mixin class, all objects can inherit from Base while
    the model functionality can be added only where needed.

    Notes
    -----
    The following table gives a brief overview of the methods that are added
    to the object by this mixin.  In addition to these methods, a ``models``
    attribute is also added, which is a dictionary that contains all of the
    models and their parameters.

    +----------------------+--------------------------------------------------+
    | Method or Attribute  | Functionality                                    |
    +======================+==================================================+
    | ``add_model``        | Add a given model and parameters to the object   |
    +----------------------+--------------------------------------------------+
    | ``regenerate_model`` | Runs the model(s) to recalculate data            |
    +----------------------+--------------------------------------------------+
    | ``remove_model``     | Removes specified model as well as it's data     |
    +----------------------+--------------------------------------------------+

    Examples
    --------
    >>> import openpnm as op

    Create a Demo class using Base and ModelsMixin:

    >>> class Demo(op.core.Base, op.core.ModelsMixin):
    ...     pass
    >>> temp = Demo(Np=3, Nt=2)

    The new class has the normal Base methods:

    >>> print(temp.num_pores())
    3

    But also has those needed for working with models.  For instance, a simple
    model can be added as follows:

    >>> temp.add_model(propname='pore.test',
    ...                model=op.models.misc.constant,
    ...                value=2)
    >>> print(temp['pore.test'])
    [2 2 2]

    All the models and their respective parameters are stored in the ``models``
    attribute:

    >>> print(temp.models)
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
    #   Property Name             Parameter                 Value
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
    1   pore.test                 model:                    constant
                                  value:                    2
                                  regeneration mode:        normal
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――

    """

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
            for k, v in zip(keys, vals):  # Put defaults into kwargs
                if k not in kwargs:  # Skip if argument was given in kwargs
                    kwargs.update({k: v})
        self.models[propname] = kwargs  # Store all keyword argumnents in model
        # Regenerate model values if necessary
        if regen_mode not in ['deferred']:
            self._regen(propname)

    def regenerate_models(self, propnames=None, exclude=[], deep=False):
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

        deep : boolean
            Specifies whether or not to regenerate models on all associated
            objects.  For instance, if ``True``, then all Physics models will
            be regenerated when method is called on the corresponding Phase.
            The default is ``False``.  The method does not work in reverse,
            so regenerating models on a Physics will not update a Phase.

        """
        # If empty list of propnames was given, do nothing and return
        if type(propnames) is list and len(propnames) == 0:
            return
        if type(propnames) is str:  # Convert string to list if necessary
            propnames = [propnames]
        if propnames is None:  # If no props given, then regenerate them all
            propnames = self.models.dependency_list()
            # If some props are to be excluded, remove them from list
            if len(exclude) > 0:
                propnames = [i for i in propnames if i not in exclude]
        # Re-order given propnames according to dependency tree
        self_models = self.models.dependency_list()
        propnames = [i for i in self_models if i in propnames]

        if deep:
            other_models = None  # Will trigger regen of ALL models
        else:
            # Make list of given propnames that are not in self
            other_models = list(set(propnames).difference(set(self_models)))
        # The following has some redundant lines, but is easier to understand
        if self._isa('phase'):
            # Start be regenerating models on self
            for item in propnames:
                self._regen(item)
            # Then regen models on associated objects, if any in other_models
            for phys in self.project.find_physics(phase=self):
                phys.regenerate_models(propnames=other_models, deep=False)
        elif self._isa('network'):  # Repeat for other object types
            for item in propnames:
                self._regen(item)
            for geom in self.project.geometries().values():
                geom.regenerate_models(propnames=other_models, deep=False)
        else:
            for item in propnames:
                self._regen(item)

    def _regen(self, prop):
        # Create a temporary dict of all model arguments
        try:
            kwargs = self.models[prop].copy()
        except KeyError:
            logger.info(prop+' not found, will retry if deep is True')
            return
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
            except KeyError:
                # Find names of missing dependencies and print nice warning
                missing_deps = []
                for key in kwargs.values():
                    if type(key) == str and key.split('.')[0] in ['pore', 'throat']:
                        try:
                            self[key]
                        except KeyError:
                            missing_deps.append(key)

                logger.warning(prop + ' was not run since the following ' +
                               'properties are missing: ' + str(missing_deps))
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
