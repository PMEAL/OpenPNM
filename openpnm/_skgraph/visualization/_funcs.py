import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib import colors as mcolors
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from openpnm._skgraph.tools import dimensionality
from openpnm._skgraph.tools import get_node_prefix, get_edge_prefix


__all__ = [
    'plot_edges',
    'plot_nodes',
]


def plot_edges(network,
               edges=None,
               ax=None,
               size_by=None,
               color_by=None,
               cmap='jet',
               color='b',
               alpha=1.0,
               linestyle='solid',
               linewidth=1,
               **kwargs):  # pragma: no cover
    r"""
    Produce a 3D plot of the network topology

    This shows how edges connect for quick visualization without having
    to export data to veiw in Paraview.

    Parameters
    ----------
    network : dict
        The network dictionary
    edges : array_like (optional)
        The list of edges to plot if only a sub-sample is desired.  This is
        useful for inspecting a small region of the network.  If no edges are
        specified then all are shown.
    fig : Matplotlib figure handle and line property arguments (optional)
        If a ``fig`` is supplied, then the topology will be overlaid on this
        plot.  This makes it possible to combine coordinates and connections,
        and to color edges differently for instance.
    size_by : array_like (optional)
        An ndarray of edge values (e.g. alg['throat.rate']).  These
        values are used to scale the ``linewidth``, so if the lines are too
        thin, then increase ``linewidth``.
    color_by : str or array_like (optional)
        An ndarray of edge values (e.g. alg['edge.rate']).
    cmap : str or cmap object (optional)
        The matplotlib colormap to use if specfying an edge property
        for ``color_by``
    color : str, optional (optional)
        A matplotlib named color (e.g. 'r' for red).
    alpha : float (optional)
        The transparency of the lines, with 1 being solid and 0 being invisible
    linestyle : str (optional)
        Can be one of {'solid', 'dashed', 'dashdot', 'dotted'}.  Default is
        'solid'.
    linewidth : float (optional)
        Controls the thickness of drawn lines.  Is used to scale the thickness
        if ``size_by`` is given. Default is 1. If a value is provided for
        ``size_by`` then they are used to scale the ``linewidth``.
    **kwargs : dict
        All other keyword arguments are passed on to the ``Line3DCollection``
        class of matplotlib, so check their documentation for additional
        formatting options.

    Returns
    -------
    lc : LineCollection or Line3DCollection
        Matplotlib object containing the lines representing the throats.

    Notes
    -----
    To create a single plot containing both coordinates and connections,
    consider creating an empty figure and then passing the ``ax`` object as
    an argument to ``plot_connections`` and ``plot_coordinates``.
    Otherwise, each call to either of these methods creates a new figure.

    See Also
    --------
    plot_coordinates

    """
    node_prefix = get_node_prefix(network)
    edge_prefix = get_edge_prefix(network)
    conns = network[edge_prefix+'.conns']
    coords = network[node_prefix+'.coords']
    Ts = np.arange(conns.shape[0]) if edges is None else edges
    dim = dimensionality(network)
    ThreeD = True if dim.sum() == 3 else False
    # Add a dummy axis for 1D networks
    if dim.sum() == 1:
        dim[np.argwhere(~dim)[0]] = True

    if "fig" in kwargs.keys():
        raise Exception("'fig' argument is deprecated, use 'ax' instead.")
    if ax is None:
        fig, ax = plt.subplots()
    else:
        # The next line is necessary if ax was created using plt.subplots()
        fig, ax = ax.get_figure(), ax.get_figure().gca()
    if ThreeD and ax.name != '3d':
        fig.delaxes(ax)
        ax = fig.add_subplot(111, projection='3d')

    # Collect coordinates
    Ps = np.unique(conns[Ts])
    X, Y, Z = coords[Ps].T
    xyz = coords[:, dim]
    P1, P2 = conns[Ts].T
    throat_pos = np.column_stack((xyz[P1], xyz[P2])).reshape((Ts.size, 2, dim.sum()))

    # Deal with optional style related arguments
    if 'c' in kwargs.keys():
        color = kwargs.pop('c')
    color = mcolors.to_rgb(color) + tuple([alpha])
    # Override colors with color_by if given
    if color_by is not None:
        color = cm.get_cmap(name=cmap)(color_by / color_by.max())
        color[:, 3] = alpha
    if size_by is not None:
        linewidth = size_by / size_by.max() * linewidth

    if ThreeD:
        lc = Line3DCollection(throat_pos, colors=color, cmap=cmap,
                              linestyles=linestyle, linewidths=linewidth,
                              antialiaseds=np.ones_like(Ts), **kwargs)
    else:
        lc = LineCollection(throat_pos, colors=color, cmap=cmap,
                            linestyles=linestyle, linewidths=linewidth,
                            antialiaseds=np.ones_like(Ts), **kwargs)
    ax.add_collection(lc)

    _scale_axes(ax=ax, X=X, Y=Y, Z=Z)
    _label_axes(ax=ax, X=X, Y=Y, Z=Z)
    fig.tight_layout()

    return lc


def plot_nodes(network,
               nodes=None,
               ax=None,
               size_by=None,
               color_by=None,
               cmap='jet',
               color='r',
               alpha=1.0,
               marker='o',
               markersize=10,
               **kwargs):  # pragma: no cover
    r"""
    Produce a 3D plot showing specified nodecoordinates as markers.

    Parameters
    ----------
    network : dict
        The network dictionary
    nodes : array_like (optional)
        The list of nodes to plot if only a sub-sample is desired. This is
        useful for inspecting a small region of the network. If no nodes
        are specified then all are shown.
    ax : Matplotlib axis handle
        If ``ax`` is supplied, then the coordinates will be overlaid.
        This enables the plotting of multiple different sets of nodes as
        well as edge connections from ``plot_connections``.
    size_by : array_like
        An ndarray of node values (e.g. alg['node.radius']). These
        values are normalized by scaled by ``markersize``.
    color_by : array_like
        An ndarray of node values (e.g. alg['node.radius']).
    cmap : str or cmap object
        The matplotlib colormap to use if specfying a node property
        for ``color_by``
    color : str
        A matplotlib named color (e.g. 'r' for red).
    alpha : float
        The transparency of the lines, with 1 being solid and 0 being invisible
    marker : 's'
        The marker to use.  The default is a circle.  Options are explained
        `here <https://matplotlib.org/3.2.1/api/markers_api.html>`_
    markersize : scalar
        Controls size of marker, default is 1.0.  This value is used to scale
        the ``size_by`` argument if given.
    **kwargs
        All other keyword arguments are passed on to the ``scatter``
        function of matplotlib, so check their documentation for additional
        formatting options.

    Returns
    -------
    pc : PathCollection
        Matplotlib object containing the markers representing the pores.

    Notes
    -----
    To create a single plot containing both coordinates and connections,
    consider creating an empty figure and then passing the ``ax`` object as
    an argument to ``plot_edges`` and ``plot_nodes``.
    Otherwise, each call to either of these methods creates a new figure.

    See Also
    --------
    plot_edges

    """
    node_prefix = get_node_prefix(network)
    coords = network[node_prefix+'.coords']
    Ps = np.arange(coords.shape[0]) if nodes is None else nodes
    dim = dimensionality(network)
    ThreeD = True if dim.sum() == 3 else False
    # Add a dummy axis for 1D networks
    if dim.sum() == 1:
        dim[np.argwhere(~dim)[0]] = True
    # Add 2 dummy axes for 0D networks (1 pore only)
    if dim.sum() == 0:
        dim[[0, 1]] = True

    if "fig" in kwargs.keys():
        raise Exception("'fig' argument is deprecated, use 'ax' instead.")
    if ax is None:
        fig, ax = plt.subplots()
    else:
        # The next line is necessary if ax was created using plt.subplots()
        fig, ax = ax.get_figure(), ax.get_figure().gca()
    if ThreeD and ax.name != '3d':
        fig.delaxes(ax)
        ax = fig.add_subplot(111, projection='3d')

    # Collect specified coordinates
    X, Y, Z = coords[Ps].T
    # The bounding box for fig is the entire ntwork (to fix the problem with
    # overwriting figures' axes lim)
    Xl, Yl, Zl = coords.T

    # Parse formatting kwargs
    if 'c' in kwargs.keys():
        color = kwargs.pop('c')
    if 's' in kwargs.keys():
        markersize = kwargs.pop('s')
    if color_by is not None:
        color = cm.get_cmap(name=cmap)(color_by / color_by.max())
    if size_by is not None:
        markersize = size_by / size_by.max() * markersize

    if ThreeD:
        sc = ax.scatter(X, Y, Z, c=color, s=markersize, marker=marker,
                        alpha=alpha, **kwargs)
        _scale_axes(ax=ax, X=Xl, Y=Yl, Z=Zl)
    else:
        _X, _Y = np.column_stack((X, Y, Z))[:, dim].T
        sc = ax.scatter(_X, _Y, c=color, s=markersize, marker=marker,
                        alpha=alpha, **kwargs)
        _scale_axes(ax=ax, X=Xl, Y=Yl, Z=np.zeros_like(Yl))

    _label_axes(ax=ax, X=Xl, Y=Yl, Z=Zl)
    fig.tight_layout()

    return sc


def _label_axes(ax, X, Y, Z):
    labels = ["X", "Y", "Z"]
    dim = np.zeros(3, dtype=bool)
    for i, arr in enumerate([X, Y, Z]):
        if np.unique(arr).size > 1:
            dim[i] = True
    # Add a dummy axis for 1D networks
    if dim.sum() == 1:
        dim[np.argwhere(~dim)[0]] = True
    # Add 2 dummy axes for 0D networks (1 pore only)
    if dim.sum() == 0:
        dim[[0, 1]] = True
    dim_idx = np.argwhere(dim).squeeze()
    ax.set_xlabel(labels[dim_idx[0]])
    ax.set_ylabel(labels[dim_idx[1]])
    if hasattr(ax, "set_zlim"):
        ax.set_zlabel("Z")


def _scale_axes(ax, X, Y, Z):
    max_range = np.ptp([X, Y, Z]).max() / 2
    mid_x = (X.max() + X.min()) * 0.5
    mid_y = (Y.max() + Y.min()) * 0.5
    mid_z = (Z.max() + Z.min()) * 0.5
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    if hasattr(ax, "set_zlim"):
        ax.set_zlim(mid_z - max_range, mid_z + max_range)
    else:
        ax.axis("equal")
