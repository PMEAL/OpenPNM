import inspect
from collections.abc import Iterable
from pnmlib.core import set_data, get_group
from pnmlib import _reserved_prefixes


__all__ = [
    "apply_models",
    "get_model_args",
    "_inspect_model",
    "dependency_graph",
    "dependency_list",
    "dependency_map",
]


def apply_models(target, models, group=''):
    keys = dependency_list(models)
    for k in keys:
        for prop in models.keys():
            if prop.startswith(k):
                v = models[prop]
                kw = get_model_args(v)
                kw['target'] = get_group(target, name=group)
                val = v['model'](**kw)
                p = group + '/' + prop
                set_data(target, key=p, value=val)


def _inspect_model(model, kwargs={}):
    if model.__defaults__:
        vals = list(inspect.getfullargspec(model).defaults)
        keys = inspect.getfullargspec(model).args[-len(vals):]
        for k, v in zip(keys, vals):  # Put defaults into kwargs
            if k not in kwargs:  # Skip if argument was given in kwargs
                kwargs.update({k: v})
    return kwargs


def get_model_args(model):
    f = model['model']
    kw = {}
    kw.update(_inspect_model(f, {}))
    args = inspect.getfullargspec(f).args
    for arg in args:
        if arg in model.keys():
            kw[arg] = model[arg]
    return kw


def is_valid_propname(propname):
    r"""
    Checks if ``propname`` is valid (i.e. starts with 'pore.' or 'throat.')

    Parameters
    ----------
    propname : str
        Property name to check whether it's a valid OpenPNM propname.

    Returns
    -------
    bool
        Whether or not ``propname`` is a valid name

    """
    if not isinstance(propname, str):
        return False
    for element in _reserved_prefixes:
        if (element + '.') in propname:
            return True
    return True


def flat_list(input_list):
    r"""
    Given a list of nested lists of arbitrary depth, returns a single
    level or 'flat' list.
    """
    def _flatten(l):
        for el in l:
            if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
                yield from _flatten(el)
            else:
                yield el

    return list(_flatten(input_list))


def dependency_graph(models):
    """
    Returns a NetworkX graph object of the dependencies

    See Also
    --------
    dependency_list
    dependency_map

    """
    import networkx as nx
    from collections.abc import Iterable

    dtree = nx.DiGraph()
    # models = list(models.keys())

    for model in models.keys():
        propname = model.split("@")[0]
        dtree.add_node(propname)
        # Filter pore/throat props only
        args = flat_list(models[model].values())
        dependencies = []
        for arg in args:
            if is_valid_propname(arg):
                dependencies.append(arg)
            elif isinstance(arg, Iterable):
                for a in arg:
                    if is_valid_propname(arg):
                        dependencies.append(arg)

        # dependencies = [arg for arg in args if is_valid_propname(arg)]
        # Add dependency from model's parameters
        for d in dependencies:
            dtree.add_edge(d, propname)

    return dtree


def dependency_list(models):
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

    dtree = dependency_graph(models)
    cycles = list(nx.simple_cycles(dtree))
    if cycles:
        msg = 'Cyclic dependency: ' + ' -> '.join(cycles[0] + [cycles[0][0]])
        raise Exception(msg)
    d = nx.algorithms.dag.lexicographical_topological_sort(dtree, sorted)
    return list(d)


def dependency_map(models,
                   ax=None,
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

    labels = {}
    node_shapes = {}
    dtree = dependency_graph(models)
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
