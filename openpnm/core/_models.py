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
        return model(target, **kwargs)

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
