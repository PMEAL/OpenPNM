import logging
import inspect
import openpnm as op
import numpy as np
from copy import deepcopy
from openpnm.utils import (
    PrintableDict,
    Workspace,
    is_valid_propname,
)


logger = logging.getLogger(__name__)
ws = Workspace()


__all__ = [
    'ModelsDict',
    'ModelWrapper',
    'ModelsMixin2',
]


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
            args = op.utils.flat_list(self[model].values())
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
    def _info(self):  # Pragma: no cover
        r"""
        Prints a nicely formatted list of model names and the domain to which
        they apply.

        Notes
        -----
        This is a hidden function for now, but could be exposed if useful.
        """
        names = {}
        for item in self:
            name, _, domain = item.partition('@')
            if name not in names.keys():
                names[name] = []
            names[name].append(domain)
        D = PrintableDict(names, key='Model', value='Domain')
        print(D)

    def __str__(self):  # pragma: no cover
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
                    super().__delitem__(item)

    def __getitem__(self, key):
        try:
            return super().__getitem__(key)
        except KeyError:
            d = PrintableDict(key='Model', value='Args')
            for k, v in self.items():
                if k.startswith(key+'@'):
                    d[k] = v
            if len(d) > 0:
                return d
            else:
                raise KeyError(key)

    def update(self, d, domain=None):
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
        model = self['model']
        kwargs = {}
        for k, v in self.items():
            if k not in ['model', 'regen_mode']:
                kwargs[k] = v
        return model(self.target, **kwargs)

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

    def __str__(self):  # pragma: no cover
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

    @property
    def target(self):
        """
        Finds and returns the object to which this model is assigned
        """
        for proj in ws.values():
            for obj in proj:
                if hasattr(obj, "models"):
                    for mod in obj.models.values():
                        if mod is self:
                            return obj
        raise Exception("No target object found!")


class ModelsMixin2:
    r"""
    This class is added to ``Network`` and ``Phase`` objects under the
    ``models`` attribute. It provides the functionality for storing and
    running pore-scale models.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.models = ModelsDict()

    def add_model(self, propname, model, domain='all', regen_mode='normal',
                  **kwargs):
        r"""
        Add a pore-scale model to the object, along with the desired arguments

        Parameters
        ----------
        propname : str
            The name of the property being computed. E.g. if
            ``propname='pore.diameter'`` then the computed results will be stored
            in ``obj['pore.diameter']``.
        model : function handle
            The function that produces the values
        domain : str
            The label indicating the locations in which results generated by
            ``model`` should be stored. See `Notes` for more details.
        regen_mode : str
            How the model should be regenerated. Options are:

            ============ =====================================================
            regen_mode   description
            ============ =====================================================
            normal       (default) The model is run immediately upon being
                         added, and is also run each time ``regenerate_models``
                         is called.
            deferred     The model is NOT run when added, but is run each time
                         ``regenerate_models`` is called. This is useful for
                         models that depend on other data that may not exist
                         yet.
            constant     The model is run immediately upon being added, but is
                         is not run when ``regenerate_models`` is called,
                         effectively turning the property into a constant.
            ============ =====================================================

        kwargs : keyword arguments
            All additional keyword arguments are passed on to the model

        Notes
        -----
        The ``domain`` argument dictates where the results of ``model`` should
        be stored. For instance, given ``propname='pore.diameter'`` and
        ``domain='left'`` then when `model` is run, the results are stored in
        in the pores labelled left. Note that if ``model`` returns ``Np``
        values, then values not belonging to ``'pore.left'`` are discarded.
        The following steps outline the process:

        1. Find the pore indices:

        .. code-block:: python

          Ps = obj.pores('left')

        2. Run the model:

        .. code-block:: python

          vals = model(**kwargs)

        3. If the model returns a full Np-length array, then extract the
        correct values and apply them to the corresponding locations:

        .. code-block:: python

          if len(vals) == obj.Np:
              obj['pore.diameter'][Ps] = vals[Ps]

        4. If the model was designed to return only the subset of values then:

        .. code-block:: python

          if len(vals) == obj.num_pores('left'):
              obj['pore.diameter'][Ps] = vals

        """
        if '@' in propname:
            propname, domain = propname.split('@')
        elif domain is None:
            domain = self.settings['default_domain']
        element, prop = propname.split('.', 1)
        domain = domain.split('.', 1)[-1]

        # Add model and regen_mode to kwargs dictionary
        kwargs.update({'model': model, 'regen_mode': regen_mode})
        # Insepct model to extract arguments and default values
        kwargs.update(self._inspect_model(model, kwargs))
        self.models[propname+'@'+domain] = ModelWrapper(**kwargs)
        if regen_mode != 'deferred':
            self.run_model(propname+'@'+domain)

    def _inspect_model(self, model, kwargs={}):
        if model.__defaults__:
            vals = list(inspect.getfullargspec(model).defaults)
            keys = inspect.getfullargspec(model).args[-len(vals):]
            for k, v in zip(keys, vals):  # Put defaults into kwargs
                if k not in kwargs:  # Skip if argument was given in kwargs
                    kwargs.update({k: v})
        return kwargs

    def add_model_collection(self, models, domain='all', regen_mode='deferred'):
        r"""
        Add a ``collection`` of several models at once

        Parameters
        ----------
        models : dict
            The collection of models to add.
        regen_mode : str
            By default the models are not regenerated upon addition. See the
            docstring for ``add_model`` for more information.
        domain : str
            The label indicating which locations the supplied collection
            of models should be applied to.

        Notes
        -----
        Collections are dictionaries that are formatted the same as
        ``obj.models``. Several model collections are available in
        ``openpnm.models.collections``.
        """
        models = deepcopy(models)
        for k, v in models.items():
            if 'domain' not in v.keys():
                v['domain'] = domain
            if 'regen_mode' not in v.keys():
                v['regen_mode'] = regen_mode
            self.add_model(propname=k, **v)

    def regenerate_models(self, propnames=None, exclude=[]):
        r"""
        Runs all the models stored in the object's ``models`` attribute

        Parameters
        ----------
        propnames : list of strings
            If given then only the specified models are run
        exclude : list of strings
            If given then these models will *not* be run

        Notes
        -----
        This function will ensure that models are called in the correct order
        such that 'pore.diameter' will be run before 'pore.volume', since
        the diameter is required to compute the volume.
        """
        all_models = self.models.dependency_list()
        # Regenerate all properties by default
        if propnames is None:
            propnames = all_models
        else:
            propnames = np.atleast_1d(propnames).tolist()
        # Remove any that are specifically excluded
        propnames = np.setdiff1d(propnames, exclude).tolist()
        # Reorder given propnames according to dependency tree
        tmp = [e.split("@")[0] for e in propnames]
        idx_sorted = [all_models.index(e) for e in tmp]
        propnames = [elem for i, elem in sorted(zip(idx_sorted, propnames))]
        # Now run each on in sequence
        for item in propnames:
            try:
                self.run_model(item)
            except KeyError as e:
                msg = (f"{item} was not run since the following property"
                       f" is missing: {e}")
                logger.warning(msg)
                self.models[item]['regen_mode'] = 'deferred'

    def run_model(self, propname, domain=None):
        r"""
        Runs the requested model and places the result into the correct
        locations

        Parameters
        ----------
        propname : str
            The name of the model to run.
        domain : str
            The label of the domain for which the model should be run. Passing
            ``propname='pore.diameter@domain1`` and ``domain=None`` is
            equivalent to passing ``propname='pore.diameter`` and
            ``domain=domain1``. Passing ``domain=None`` will regenerate
            all models starting with ``propname``.
        """
        if domain is None:
            if '@' in propname:  # Get domain from propname if present
                propname, _, domain = propname.partition('@')
                self.run_model(propname=propname, domain=domain)
            else:  # No domain means run model for ALL domains
                for item in self.models.keys():
                    if item.startswith(propname+"@"):
                        _, _, domain = item.partition("@")
                        self.run_model(propname=propname, domain=domain)
        else:  # domain was given explicitly
            domain = domain.split('.', 1)[-1]
            element, prop = propname.split('@')[0].split('.', 1)
            propname = f'{element}.{prop}'
            mod_dict = self.models[propname+'@'+domain]
            # Collect kwargs
            kwargs = {'domain': f'{element}.{domain}'}
            for item in mod_dict.keys():
                if item not in ['model', 'regen_mode']:
                    kwargs[item] = mod_dict[item]
            # Deal with models that don't have domain argument yet
            if 'domain' not in inspect.getfullargspec(mod_dict['model']).args:
                _ = kwargs.pop('domain', None)
                vals = mod_dict['model'](self, **kwargs)
                if isinstance(vals, dict):  # Handle models that return a dict
                    for k, v in vals.items():
                        v = np.atleast_1d(v)
                        if v.shape[0] == 1:  # Returned item was a scalar
                            v = np.tile(v, self._count(element))
                        vals[k] = v[self[f'{element}.{domain}']]
                elif isinstance(vals, (int, float)):  # Handle models that return a float
                    vals = np.atleast_1d(vals)
                else:  # Index into full domain result for use below
                    vals = vals[self[f'{element}.{domain}']]
            else:  # Model that accepts domain arg
                vals = mod_dict['model'](self, **kwargs)
            # Finally add model results to self
            if isinstance(vals, np.ndarray):  # If model returns single array
                if propname not in self.keys():
                    temp = self._initialize_empty_array_like(vals, element)
                    self[f'{element}.{prop}'] = temp
                self[propname][self[f'{element}.{domain}']] = vals
            elif isinstance(vals, dict):  # If model returns a dict of arrays
                for k, v in vals.items():
                    if f'{propname}.{k}' not in self.keys():
                        temp = self._initialize_empty_array_like(v, element)
                        self[f'{propname}.{k}'] = temp
                    self[f'{propname}.{k}'][self[f'{element}.{domain}']] = v
