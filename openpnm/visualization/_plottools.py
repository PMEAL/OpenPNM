import logging

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import cm
from tqdm.auto import tqdm

import openpnm as op

logger = logging.getLogger(__name__)


__all__ = [
    'plot_connections',
    'plot_coordinates',
    'plot_networkx',
    'plot_tutorial',
    'plot_notebook',
    'plot_vispy',
    'generate_voxel_image',
    'set_mpl_style',
]


def plot_connections(network,
                     throats=None,
                     ax=None,
                     size_by=None,
                     color_by=None,
                     label_by=None,
                     cmap='jet',
                     color='b',
                     alpha=1.0,
                     linestyle='solid',
                     linewidth=1,
                     **kwargs):  # pragma: no cover
    r"""
    Produce a 3D plot of the network topology.

    This shows how throats connect for quick visualization without having
    to export data to veiw in Paraview.

    Parameters
    ----------
    network : Network
        The network whose topological connections to plot
    throats : array_like (optional)
        The list of throats to plot if only a sub-sample is desired.  This is
        useful for inspecting a small region of the network.  If no throats are
        specified then all throats are shown.
    fig : Matplotlib figure handle and line property arguments (optional)
        If a ``fig`` is supplied, then the topology will be overlaid on this
        plot.  This makes it possible to combine coordinates and connections,
        and to color throats differently for instance.
    size_by : array_like (optional)
        An ndarray of throat values (e.g. alg['throat.rate']).  These
        values are used to scale the ``linewidth``, so if the lines are too
        thin, then increase ``linewidth``.
    color_by : array_like (optional)
        An ndarray of throat values (e.g. alg['throat.rate']).
    label_by : array_like (optional)
        An array or list of values to use as labels
    cmap : str or cmap object (optional)
        The matplotlib colormap to use if specfying a throat property
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
    font : dict
        A dictionary of key-value pairs that are used to control the font
        appearance if `label_by` is provided.
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
    To create a single plot containing both pore coordinates and throats,
    consider creating an empty figure and then pass the ``ax`` object as
    an argument to ``plot_connections`` and ``plot_coordinates``.
    Otherwise, each call to either of these methods creates a new figure.

    See Also
    --------
    plot_coordinates

    Examples
    --------
    >>> import openpnm as op
    >>> import matplotlib as mpl
    >>> import matplotlib.pyplot as plt
    >>> mpl.use('Agg')
    >>> pn = op.network.Cubic(shape=[10, 10, 3])
    >>> pn.add_boundary_pores()
    >>> Ts = pn.throats('*boundary', mode='not')  # find internal throats
    >>> fig, ax = plt.subplots()  # create empty figure
    >>> _ = op.visualization.plot_connections(network=pn,
    ...                                       throats=Ts)  # plot internal throats
    >>> Ts = pn.throats('*boundary')  # find boundary throats
    >>> _ = op.visualization.plot_connections(network=pn,
    ...                                       throats=Ts,
    ...                                       ax=ax,
    ...                                       color='r')  # plot boundary throats in red

    """
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib import colors as mcolors
    from matplotlib.collections import LineCollection
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    from openpnm.topotools import dimensionality

    Ts = network.Ts if throats is None else network._parse_indices(throats)
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
    Ps = np.unique(network['throat.conns'][Ts])
    X, Y, Z = network['pore.coords'][Ps].T
    xyz = network["pore.coords"][:, dim]
    P1, P2 = network["throat.conns"][Ts].T
    throat_pos = np.column_stack((xyz[P1], xyz[P2])).reshape((Ts.size, 2, dim.sum()))

    # Deal with optional style related arguments
    if 'c' in kwargs.keys():
        color = kwargs.pop('c')
    color = mcolors.to_rgb(color) + tuple([alpha])
    if isinstance(cmap, str):
        try:
            cmap = plt.colormaps.get_cmap(cmap)
        except AttributeError:
            cmap = plt.cm.get_cmap(cmap)
    # Override colors with color_by if given
    if color_by is not None:
        color_by = np.array(color_by, dtype=np.float16)
        if len(color_by) != len(Ts):
            color_by = color_by[Ts]
        if not np.all(np.isfinite(color_by)):
            color_by[~np.isfinite(color_by)] = 0
            logger.warning('nans or infs found in color_by array, setting to 0')
        vmin = kwargs.pop('vmin', color_by.min())
        vmax = kwargs.pop('vmax', color_by.max())
        cscale = (color_by - vmin) / (vmax - vmin)
        color = cmap(cscale)
        color[:, 3] = alpha
    if size_by is not None:
        if len(size_by) != len(Ts):
            size_by = size_by[Ts]
        if not np.all(np.isfinite(size_by)):
            size_by[~np.isfinite(size_by)] = 0
            logger.warning('nans or infs found in size_by array, setting to 0')
        linewidth = size_by / size_by.max() * linewidth
    if label_by is not None:
        if len(label_by) != len(Ts):
            label_by = label_by[Ts]
    fontkws = kwargs.pop('font', {})

    if ThreeD:
        lc = Line3DCollection(throat_pos, colors=color, cmap=cmap,
                              linestyles=linestyle, linewidths=linewidth,
                              antialiaseds=np.ones_like(network.Ts), **kwargs)
    else:
        lc = LineCollection(throat_pos, colors=color, cmap=cmap,
                            linestyles=linestyle, linewidths=linewidth,
                            antialiaseds=np.ones_like(network.Ts), **kwargs)
        if label_by is not None:
            for count, (P1, P2) in enumerate(network.conns[Ts, :]):
                i, j, k = np.mean(network.coords[[P1, P2], :], axis=0)
                ax.text(i, j, label_by[count],
                        ha='center', va='center',
                        **fontkws)
    ax.add_collection(lc)

    if np.size(Ts) > 0:
        _scale_axes(ax=ax, X=X, Y=Y, Z=Z)
        _label_axes(ax=ax, X=X, Y=Y, Z=Z)
        fig.tight_layout()

    return lc


def plot_coordinates(network,
                     pores=None,
                     ax=None,
                     size_by=None,
                     color_by=None,
                     label_by=None,
                     cmap='jet',
                     color='r',
                     alpha=1.0,
                     marker='o',
                     markersize=10,
                     **kwargs):  # pragma: no cover
    r"""
    Produce a 3D plot showing specified pore coordinates as markers.

    Parameters
    ----------
    network : Network
        The network whose topological connections to plot.
    pores : array_like (optional)
        The list of pores to plot if only a sub-sample is desired. This is
        useful for inspecting a small region of the network. If no pores
        are specified then all are shown.
    ax : Matplotlib axis handle
        If ``ax`` is supplied, then the coordinates will be overlaid.
        This enables the plotting of multiple different sets of pores as
        well as throat connections from ``plot_connections``.
    size_by : str or array_like
        An ndarray of pore values (e.g. alg['pore.concentration']). These
        values are normalized by scaled by ``markersize``.  Note that this controls
        the marker *area*, so if you want the markers to be proportional to diameter
        you should do `size_by=net['pore.diameter']**2`.
    color_by : str or array_like
        An ndarray of pore values (e.g. alg['pore.concentration']).
    label_by : array_like (optional)
        An array or list of values to use as labels
    cmap : str or cmap object
        The matplotlib colormap to use if specfying a pore property
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
    font : dict
        A dictionary of key-value pairs that are used to control the font
        appearance if `label_by` is provided.
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
    To create a single plot containing both pore coordinates and throats,
    consider creating an empty figure and then pass the ``ax`` object as
    an argument to ``plot_connections`` and ``plot_coordinates``.
    Otherwise, each call to either of these methods creates a new figure.

    See Also
    --------
    plot_connections

    Examples
    --------
    >>> import openpnm as op
    >>> import matplotlib as mpl
    >>> import matplotlib.pyplot as plt
    >>> mpl.use('Agg')
    >>> pn = op.network.Cubic(shape=[10, 10, 3])
    >>> pn['pore.internal'] = True
    >>> pn.add_boundary_pores()
    >>> Ps = pn.pores('internal')  # find internal pores
    >>> fig, ax = plt.subplots()  # create empty figure
    >>> _ = op.visualization.plot_coordinates(network=pn,
    ...                                       pores=Ps,
    ...                                       color='b',
    ...                                       ax=ax)  # plot internal pores
    >>> Ps = pn.pores('*boundary')  # find boundary pores
    >>> _ = op.visualization.plot_coordinates(network=pn,
    ...                                       pores=Ps,
    ...                                       color='r',
    ...                                       ax=ax)  # plot boundary pores in red

    """
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from mpl_toolkits.mplot3d import Axes3D

    from openpnm.topotools import dimensionality

    Ps = network.Ps if pores is None else network._parse_indices(pores)

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
    X, Y, Z = network['pore.coords'][Ps].T
    # The bounding box for fig is the entire ntwork (to fix the problem with
    # overwriting figures' axes lim)
    Xl, Yl, Zl = network['pore.coords'].T

    # Parse formatting kwargs
    if 'c' in kwargs.keys():
        color = kwargs.pop('c')
    if 's' in kwargs.keys():
        markersize = kwargs.pop('s')
    if isinstance(cmap, str):
        try:
            cmap = plt.colormaps.get_cmap(cmap)
        except AttributeError:
            cmap = plt.cm.get_cmap(cmap)
    if color_by is not None:
        color_by = np.array(color_by, dtype=np.float16)
        if len(color_by) != len(Ps):
            color_by = color_by[Ps]
        if not np.all(np.isfinite(color_by)):
            color_by[~np.isfinite(color_by)] = 0
            logger.warning('nans or infs found in color_by array, setting to 0')
        vmin = kwargs.pop('vmin', color_by.min())
        vmax = kwargs.pop('vmax', color_by.max())
        cscale = (color_by - vmin) / (vmax - vmin)
        color = cmap(cscale)
    if size_by is not None:
        if len(size_by) != len(Ps):
            size_by = size_by[Ps]
        if not np.all(np.isfinite(size_by)):
            size_by[~np.isfinite(size_by)] = 0
            logger.warning('nans or infs found in size_by array, setting to 0')
        markersize = size_by / size_by.max() * markersize
    if label_by is not None:
        if len(label_by) != len(Ps):
            label_by = label_by[Ps]
    fontkws = kwargs.pop('font', {})
    if ThreeD:
        sc = ax.scatter(X, Y, Z,
                        c=color,
                        s=markersize,
                        marker=marker,
                        alpha=alpha,
                        **kwargs)
        _scale_axes(ax=ax, X=Xl, Y=Yl, Z=Zl)
    else:
        _X, _Y = np.column_stack((X, Y, Z))[:, dim].T
        sc = ax.scatter(_X, _Y,
                        c=color,
                        s=markersize,
                        marker=marker,
                        alpha=alpha,
                        **kwargs)
        if label_by is not None:
            for count, (i, j, k) in enumerate(network.coords[Ps, :]):
                ax.text(i, j, label_by[count],
                        ha='center', va='center',
                        **fontkws)
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


def plot_networkx(network,
                  plot_throats=True,
                  labels=None,
                  colors=None,
                  scale=1,
                  ax=None,
                  alpha=1.0):  # pragma: no cover
    r"""
    Creates a pretty 2D plot for 2D OpenPNM networks.

    Parameters
    ----------
    network : Network
    plot_throats : bool, optional
        Plots throats as well as pores, if True.
    labels : list, optional
        List of OpenPNM labels
    colors : list, optional
        List of corresponding colors to the given `labels`.
    scale : float, optional
        Scale factor for size of pores.
    ax : matplotlib.Axes, optional
        Matplotlib axes object
    alpha: float, optional
        Transparency value, 1 is opaque and 0 is transparent

    """
    import matplotlib.pyplot as plt
    from matplotlib.collections import PathCollection
    from networkx import Graph, draw_networkx_edges, draw_networkx_nodes

    from openpnm.topotools import dimensionality

    dims = dimensionality(network)
    if dims.sum() > 2:
        raise Exception("NetworkX plotting only works for 2D networks.")
    temp = network['pore.coords'].T[dims].squeeze()
    if dims.sum() == 1:
        x = temp
        y = np.zeros_like(x)
    if dims.sum() == 2:
        x, y = temp
    try:
        node_size = scale * network['pore.diameter']
    except KeyError:
        node_size = np.ones_like(x) * scale * 0.5
    G = Graph()
    pos = {network.Ps[i]: [x[i], y[i]] for i in range(network.Np)}
    if not np.isfinite(node_size).all():
        node_size[~np.isfinite(node_size)] = np.nanmin(node_size)
    node_color = np.array(['k'] * len(network.Ps))

    if labels:
        if not isinstance(labels, list):
            labels = [labels]
        if not isinstance(colors, list):
            colors = [colors]
        if len(labels) != len(colors):
            raise Exception('len(colors) must be equal to len(labels)!')
        for label, color in zip(labels, colors):
            node_color[network.pores(label)] = color

    if ax is None:
        fig, ax = plt.subplots()
    ax.set_aspect('equal', adjustable='datalim')
    offset = node_size.max() * 0.5
    ax.set_xlim((x.min() - offset, x.max() + offset))
    ax.set_ylim((y.min() - offset, y.max() + offset))
    ax.axis("off")

    # Keep track of already plotted nodes
    temp = [id(item) for item in ax.collections if isinstance(item, PathCollection)]

    # Plot pores
    gplot = draw_networkx_nodes(G, ax=ax, pos=pos, nodelist=network.Ps.tolist(),
                                alpha=alpha, node_color=node_color, edgecolors=node_color,
                                node_size=node_size)
    # (Optionally) Plot throats
    if plot_throats:
        draw_networkx_edges(G, pos=pos, edge_color='k', alpha=alpha,
                            edgelist=network['throat.conns'].tolist(), ax=ax)

    spi = 2700  # 1250 was obtained by trial and error
    figwidth, figheight = ax.get_figure().get_size_inches()
    figsize_ratio = figheight / figwidth
    data_ratio = ax.get_data_ratio()
    corr = min(figsize_ratio / data_ratio, 1)
    xrange = np.ptp(ax.get_xlim())
    markersize = np.atleast_1d((corr*figwidth)**2 / xrange**2 * node_size**2 * spi)
    for item in ax.collections:
        if isinstance(item, PathCollection) and id(item) not in temp:
            item.set_sizes(markersize)

    return gplot


def plot_tutorial(network,
                  pore_labels=None,
                  throat_labels=None,
                  font_size=12,
                  line_width=2,
                  node_color='b',
                  edge_color='r',
                  node_size=500):  # pragma: no cover
    r"""
    Generate a network plot suitable for tutorials and explanations.

    Parameters
    ----------
    network : Network
        The network to plot, should be 2D, since the z-coordinate will be
        ignored.
    pore_labels : array_like
        A list of values to use for labeling the pores. If not provided then pore
        index is used.
    throat_labels : array_like
        A list of values to use for labeling the throat. If not provided then throat
        index is used.
    font_size : int
        Size of font to use for labels.
    line_width : int
        Thickness of edge lines and node borders.
    node_color : str
        Color of node border.
    edge_color : str
        Color of edge lines.
    node_size : int
        Size of node circle.

    Returns
    -------
    g : NetworkX plot object

    """
    import matplotlib.pyplot as plt
    import networkx as nx

    from openpnm.io import network_to_networkx

    G = network_to_networkx(network=network)
    pos = {i: network['pore.coords'][i, 0:2] for i in network.Ps}
    if pore_labels is None:
        labels = {i: i for i in network.Ps}
    else:
        labels = {i: pore_labels[i] for i in network.Ps}
    if throat_labels is None:
        edge_labels = {tuple(network['throat.conns'][i, :]): i for i in network.Ts}
    else:
        edge_labels = {tuple(network['throat.conns'][i, :]): throat_labels[i]
                       for i in network.Ts}

    gplot = nx.draw_networkx_nodes(G, pos,
                                   node_size=node_size,
                                   node_color='w',
                                   edgecolors=node_color,
                                   linewidths=line_width)
    nx.draw_networkx_edges(
        G, pos, width=line_width, edge_color=edge_color)
    nx.draw_networkx_labels(
        G, pos, labels=labels, font_size=font_size, font_color='k')
    nx.draw_networkx_edge_labels(
        G, pos, edge_labels=edge_labels, font_size=font_size, font_color='k')

    # Prettify the figure (margins, etc.)
    plt.axis('off')
    ax = plt.gca()
    ax.margins(0.1, 0.1)
    ax.set_aspect("equal")
    fig = plt.gcf()
    fig.tight_layout()
    dims = op.topotools.dimensionality(network)
    xy_range = np.ptp(network.coords, axis=0)[dims]
    aspect_ratio = xy_range[0] / xy_range[1]
    fig.set_size_inches(5, 5 / aspect_ratio)

    return gplot


def plot_notebook(network,
                  node_color=0,
                  edge_color=0,
                  node_size=1,
                  node_scale=20,
                  edge_scale=5,
                  colormap='viridis'):
    r"""
    Visualize a network in 3D using Plotly.

    The pores and throats are scaled and colored by their properties.
    The final figure can be rotated and zoomed.

    Parameters
    ----------
    network : Network
        The network to visualize
    node_color : ndarray
        An array of values used for coloring the pores. If not given, the
        lowest value of the employed colormap is assigned to all markers.
    edge_color : ndarray
        An array of values used for coloring the throats. If not given, the
        lowest value of the employed colormap is assigned to all lines.
    node_size : ndarray
        An array of values controlling the size of the markers.  If not given
        all markers will be the same size
    node_scale : scalar
        A scaler to resize the markers
    edge_scale : scalar
        A scaler to the line thickness
    colormap : str
        The colormap to use

    Returns
    -------
    fig : Plotly graph object
        The graph object containing the generated plots. The object has
        several useful methods.

    Notes
    -----
    **Important**

    a) This does not work in Spyder. It should only be called from a
    Jupyter Notebook.

    b) This is only meant for relatively small networks. For proper
    visualization use Paraview.

    """
    try:
        import plotly.graph_objects as go
    except ImportError:
        raise Exception('Plotly is not installed.'
                        'Please install Plotly using "pip install plotly"')

    # Get xyz coords for points
    x_nodes, y_nodes, z_nodes = network.coords.T

    node_size = np.ones(network.Np)*node_size
    node_color = np.ones(network.Np)*node_color
    edge_color = np.ones(network.Nt)*edge_color

    node_labels = [str(i) + ': ' + str(x) for i, x in
                   enumerate(zip(node_size, node_color))]
    edge_labels = [str(i) + ': ' + str(x) for i, x in enumerate(edge_color)]

    # Create edges and nodes coordinates
    N = network.Nt*3

    x_edges = np.zeros(N)
    x_edges[np.arange(0, N, 3)] = network.coords[network.conns[:, 0]][:, 0]
    x_edges[np.arange(1, N, 3)] = network.coords[network.conns[:, 1]][:, 0]
    x_edges[np.arange(2, N, 3)] = np.nan

    y_edges = np.zeros(network.Nt*3)
    y_edges[np.arange(0, N, 3)] = network.coords[network.conns[:, 0]][:, 1]
    y_edges[np.arange(1, N, 3)] = network.coords[network.conns[:, 1]][:, 1]
    y_edges[np.arange(2, N, 3)] = np.nan

    z_edges = np.zeros(network.Nt*3)
    z_edges[np.arange(0, N, 3)] = network.coords[network.conns[:, 0]][:, 2]
    z_edges[np.arange(1, N, 3)] = network.coords[network.conns[:, 1]][:, 2]
    z_edges[np.arange(2, N, 3)] = np.nan

    # Create plotly's Scatter3d object for pores and throats
    trace_edges = go.Scatter3d(x=x_edges,
                               y=y_edges,
                               z=z_edges,
                               mode='lines',
                               line=dict(color=edge_color,
                                         width=edge_scale,
                                         colorscale=colormap),
                               text=edge_labels, hoverinfo='text')

    trace_nodes = go.Scatter3d(x=x_nodes,
                               y=y_nodes,
                               z=z_nodes,
                               mode='markers',
                               marker=dict(symbol='circle',
                                           size=node_size*node_scale,
                                           color=node_color,
                                           colorscale=colormap,
                                           line=dict(color='black', width=0.5)),
                               text=node_labels, hoverinfo='text')

    axis = dict(showbackground=False,
                showline=False,
                zeroline=False,
                showgrid=False,
                showticklabels=False,
                title='')

    layout = go.Layout(width=650,
                       height=625,
                       showlegend=False,
                       scene=dict(xaxis=dict(axis),
                                  yaxis=dict(axis),
                                  zaxis=dict(axis),),
                       margin=dict(t=100),
                       hovermode='closest')

    data = [trace_edges, trace_nodes]
    fig = go.Figure(data=data, layout=layout)
    return fig


def _generate_voxel_image(network, pore_shape, throat_shape, max_dim=200):
    r"""
    Generates a 3d numpy array from an OpenPNM network

    Parameters
    ----------
    network : OpenPNM Network
        Network from which voxel image is to be generated
    pore_shape : str
        Shape of pores in the network, valid choices are "sphere", "cube"
    throat_shape : str
        Shape of throats in the network, valid choices are "cylinder", "cuboid"
    max_dim : int
        Number of voxels in the largest dimension of the network

    Returns
    -------
    im : ndarray
        Voxelated image corresponding to the given pore network model

    Notes
    -----
    (1) The generated voxel image is labeled with 0s, 1s and 2s signifying
    solid phase, pores, and throats respectively.

    """
    from porespy.tools import insert_cylinder, overlay
    from skimage.morphology import ball, cube
    xyz = network["pore.coords"]
    cn = network["throat.conns"]

    # Distance bounding box from the network by a fixed amount
    delta = network["pore.diameter"].mean() / 2
    if isinstance(network, op.network.Cubic):
        try:
            delta = op.topotools.get_spacing(network).mean() / 2
        except AttributeError:
            delta = network.spacing.mean() / 2

    # Shift everything to avoid out-of-bounds
    extra_clearance = int(max_dim * 0.05)

    # Transform points to satisfy origin at (0, 0, 0)
    xyz0 = xyz.min(axis=0) - delta
    xyz += -xyz0
    res = (np.ptp(xyz, axis=0).max() + 2 * delta) / max_dim
    shape = np.rint((xyz.max(axis=0) + delta) / res).astype(int) + 2 * extra_clearance

    # Transforming from real coords to matrix coords
    xyz = np.rint(xyz / res).astype(int) + extra_clearance
    pore_radi = np.rint(network["pore.diameter"] * 0.5 / res).astype(int)
    throat_radi = np.rint(network["throat.diameter"] * 0.5 / res).astype(int)

    im_pores = np.zeros(shape, dtype=np.uint8)
    im_throats = np.zeros_like(im_pores)

    if pore_shape == "cube":
        pore_elem = cube
        rp = pore_radi * 2 + 1  # +1 since num_voxel must be odd
        rp_max = int(2 * round(delta / res)) + 1
    if pore_shape == "sphere":
        pore_elem = ball
        rp = pore_radi
        rp_max = int(round(delta / res))
    if throat_shape == "cuboid":
        raise Exception("Not yet implemented, try 'cylinder'.")

    # Generating voxels for pores
    for i, pore in enumerate(tqdm(network.Ps)):
        elem = pore_elem(rp[i])
        try:
            im_pores = overlay(im1=im_pores, im2=elem, c=xyz[i])
        except ValueError:
            elem = pore_elem(rp_max)
            im_pores = overlay(im1=im_pores, im2=elem, c=xyz[i])
    # Get rid of pore overlaps
    im_pores[im_pores > 0] = 1

    # Generating voxels for throats
    for i, throat in enumerate(tqdm(network.Ts)):
        try:
            im_throats = insert_cylinder(
                im_throats, r=throat_radi[i], xyz0=xyz[cn[i, 0]], xyz1=xyz[cn[i, 1]])
        except ValueError:
            im_throats = insert_cylinder(
                im_throats, r=rp_max, xyz0=xyz[cn[i, 0]], xyz1=xyz[cn[i, 1]])
    # Get rid of throat overlaps
    im_throats[im_throats > 0] = 1

    # Subtract pore-throat overlap from throats
    im_throats = (im_throats.astype(bool) * ~im_pores.astype(bool)).astype(np.uint8)
    im = im_pores * 1 + im_throats * 2

    return im[extra_clearance:-extra_clearance,
              extra_clearance:-extra_clearance,
              extra_clearance:-extra_clearance]


def generate_voxel_image(network, pore_shape="sphere", throat_shape="cylinder",
                         max_dim=None, rtol=0.1):
    r"""
    Generate a voxel image from a Network

    Parameters
    ----------
    network : OpenPNM Network
        Network from which voxel image is to be generated
    pore_shape : str
        Shape of pores in the network, valid choices are "sphere", "cube"
    throat_shape : str
        Shape of throats in the network, valid choices are "cylinder", "cuboid"
    max_dim : int
        Number of voxels in the largest dimension of the network
    rtol : float
        Stopping criteria for finding the smallest voxel image such that
        further increasing the number of voxels in each dimension by 25% would
        improve the predicted porosity of the image by less that ``rtol``

    Returns
    -------
    im : ndarray
        Voxelated image corresponding to the given pore network model

    Notes
    -----
    (1) The generated voxelated image is labeled with 0s, 1s and 2s signifying
    solid phase, pores, and throats respectively.

    (2) If max_dim is not provided, the method calculates it such that the
    further increasing it doesn't change porosity by much.

    """
    # If max_dim is provided, generate voxel image using max_dim
    if max_dim is not None:
        return _generate_voxel_image(
            network, pore_shape, throat_shape, max_dim=max_dim)
    max_dim = 200
    # If max_dim is not provided, find best max_dim that predicts porosity
    err = 100
    eps_old = 200
    while err > rtol:
        im = _generate_voxel_image(
            network, pore_shape, throat_shape, max_dim=max_dim)
        eps = im.astype(bool).sum() / np.prod(im.shape)
        err = abs(1 - eps / eps_old)
        eps_old = eps
        max_dim = int(max_dim * 1.25)
    return im


def create_pore_colors_from_array(a, cmap='viridis'):
    colormap = cm.get_cmap(cmap)
    return colormap(a/a.max())


def create_throat_colors_from_array(a, cmap='viridis'):
    colormap = cm.get_cmap(cmap)
    return np.repeat(colormap(a/a.max()), 2, axis=0)


def plot_vispy(
    network,
    pore_color=None,
    pore_size=None,
    throat_color=None,
    throat_size=None,
    bgcolor='grey',
):
    r"""
    Creates a pretty network plot using VisPy.

    Parameters
    ----------
    network
    """
    try:
        from vispy import scene
    except ModuleNotFoundError:
        raise Exception("vispy must be installed to use this function")
    canvas = scene.SceneCanvas(keys='interactive', show=True, bgcolor=bgcolor)
    view = canvas.central_widget.add_view()
    view.camera = 'turntable'
    view.camera.fov = 30
    view.camera.distance = 3*np.max(network['pore.coords'])

    if pore_color is None:
        pore_color = create_pore_colors_from_array(network['pore.diameter'],
                                                   cmap='viridis')
    else:
        pore_color = create_pore_colors_from_array(pore_color,
                                                   cmap='viridis')
    if throat_color is None:
        throat_color = create_throat_colors_from_array(network['throat.diameter'],
                                                       cmap='viridis')
    else:
        throat_color = create_throat_colors_from_array(throat_color,
                                                       cmap='viridis')
    if pore_size is None:
        pore_size = network['pore.diameter']
    if throat_size is None:
        throat_size = 2
    else:
        throat_size = np.max(throat_size)  # Arrays not supported here
    # plot spheres
    vis = scene.visuals.Markers(
        pos=network['pore.coords'],
        size=pore_size,
        antialias=0,
        face_color=pore_color,
        edge_width=0,
        scaling=True,
        spherical=True,
    )
    vis.parent = view.scene
    # plot axis
    # vispy.scene.visuals.XYZAxis(parent=view.scene)
    # set camera center
    view.camera.center = np.array((network['pore.coords'][:, 0].max()/2,
                                   network['pore.coords'][:, 1].max()/2,
                                   network['pore.coords'][:, 2].max()/2))
    # data preparation
    lines = np.zeros((len(network['throat.conns']), 2, 3))
    for i in range(len(network['throat.conns'])):
        pair = network['throat.conns'][i]
        line = np.array([[network['pore.coords'][pair[0]],
                          network['pore.coords'][pair[1]]]])
        lines[i, :, :] = line
    # plot throats
    vis2 = scene.visuals.Line(lines,
                              width=throat_size,
                              color=throat_color,
                              connect='segments',
                              antialias=True,)
    vis2.parent = view.scene


def set_mpl_style():  # pragma: no cover
    r"""
    Prettifies matplotlib's output by adjusting fonts, markersize etc.
    """
    sfont = 12
    mfont = 12
    lfont = 12

    image_props = {'interpolation': 'none',
                   'cmap': 'viridis'}
    line_props = {'linewidth': 2,
                  'markersize': 7,
                  'markerfacecolor': 'w'}
    font_props = {'size': sfont}
    axes_props = {'titlesize': lfont,
                  'labelsize': mfont,
                  'linewidth': 2,
                  'labelpad': 8}
    xtick_props = {'labelsize': sfont,
                   'top': True,
                   'direction': 'in',
                   'major.size': 6,
                   'major.width': 2}
    ytick_props = {'labelsize': sfont,
                   'right': True,
                   'direction': 'in',
                   'major.size': 6,
                   'major.width': 2}
    legend_props = {'fontsize': mfont,
                    'frameon': False}
    figure_props = {'titlesize': sfont,
                    'autolayout': True}

    plt.rc('font', **font_props)
    plt.rc('lines', **line_props)
    plt.rc('axes', **axes_props)
    plt.rc('xtick', **xtick_props)
    plt.rc('ytick', **ytick_props)
    plt.rc('legend', **legend_props)
    plt.rc('figure', **figure_props)
    plt.rc('image', **image_props)

    try:
        import IPython
        IPython.display.set_matplotlib_formats('png2x')
    except ModuleNotFoundError:
        pass
