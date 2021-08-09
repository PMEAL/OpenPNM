import numpy as np
import openpnm as op


def plot_connections(network,
                     throats=None,
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
    Produce a 3D plot of the network topology.

    This shows how throats connect for quick visualization without having
    to export data to veiw in Paraview.

    Parameters
    ----------
    network : GenericNetwork
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
        An ND-array of throat values (e.g. alg['throat.rate']).  These
        values are used to scale the ``linewidth``, so if the lines are too
        thin, then increase ``linewidth``.
    color_by : str or array_like (optional)
        An ND-array of throat values (e.g. alg['throat.rate']).
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
    >>> Ts = pn.throats('*boundary', mode='not')        # find internal throats
    >>> fig, ax = plt.subplots()                        # create empty figure
    >>> _ = op.topotools.plot_connections(network=pn,
    ...                                   throats=Ts)   # plot internal throats
    >>> Ts = pn.throats('*boundary')                    # find boundary throats
    >>> _ = op.topotools.plot_connections(network=pn,
    ...                                   throats=Ts,
    ...                                   ax=ax,
    ...                                   color='r')    # plot boundary throats in red

    """
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib import colors as mcolors
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib.collections import LineCollection
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
    # Override colors with color_by if given
    if color_by is not None:
        color = cm.get_cmap(name=cmap)(color_by / color_by.max())
        color[:, 3] = alpha
    if size_by is not None:
        linewidth = size_by / size_by.max() * linewidth

    if ThreeD:
        lc = Line3DCollection(throat_pos, colors=color, cmap=cmap,
                              linestyles=linestyle, linewidths=linewidth,
                              antialiaseds=np.ones_like(network.Ts), **kwargs)
    else:
        lc = LineCollection(throat_pos, colors=color, cmap=cmap,
                            linestyles=linestyle, linewidths=linewidth,
                            antialiaseds=np.ones_like(network.Ts), **kwargs)
    ax.add_collection(lc)

    _scale_axes(ax=ax, X=X, Y=Y, Z=Z)
    _label_axes(ax=ax, X=X, Y=Y, Z=Z)
    fig.tight_layout()

    return lc


def plot_coordinates(network,
                     pores=None,
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
    Produce a 3D plot showing specified pore coordinates as markers.

    Parameters
    ----------
    network : OpenPNM Network Object
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
        An ND-array of pore values (e.g. alg['pore.concentration']). These
        values are normalized by scaled by ``markersize``.
    color_by : str or array_like
        An ND-array of pore values (e.g. alg['pore.concentration']).
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
    >>> pn.add_boundary_pores()
    >>> Ps = pn.pores('internal')                       # find internal pores
    >>> fig, ax = plt.subplots()                        # create empty figure
    >>> _ = op.topotools.plot_coordinates(network=pn,
    ...                                   pores=Ps,
    ...                                   color='b',
    ...                                   ax=ax)        # plot internal pores
    >>> Ps = pn.pores('*boundary')                      # find boundary pores
    >>> _ = op.topotools.plot_coordinates(network=pn,
    ...                                   pores=Ps,
    ...                                   color='r',
    ...                                   ax=ax)        # plot boundary pores in red

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
    if color_by is not None:
        color = cm.get_cmap(name=cmap)(color_by / color_by.max())
    if size_by is not None:
        markersize = size_by / size_by.max() * markersize

    if ThreeD:
        sc = ax.scatter(X, Y, Z, c=color, s=markersize, marker=marker, alpha=alpha, **kwargs)
        _scale_axes(ax=ax, X=Xl, Y=Yl, Z=Zl)
    else:
        _X, _Y = np.column_stack((X, Y, Z))[:, dim].T
        sc = ax.scatter(_X, _Y, c=color, s=markersize, marker=marker, alpha=alpha, **kwargs)
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
    Creates a pretty 2d plot for 2d OpenPNM networks.

    Parameters
    ----------
    network : GenericNetwork
    plot_throats : boolean, optional
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
    from networkx import Graph, draw_networkx_nodes, draw_networkx_edges
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
        raise Exception('nan/inf values found in network["pore.diameter"]')
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
                                alpha=alpha, node_color="w", edgecolors=node_color,
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


def plot_vpython(network,
                 Psize='pore.diameter',
                 Tsize='throat.diameter',
                 Pcolor=None,
                 Tcolor=None,
                 cmap='jet',
                 **kwargs):  # pragma: no cover
    r"""
    Quickly visualize a network in 3D using VPython.

    Parameters
    ----------
    network : GenericNetwork
        The network to visualize.
    Psize : str (default = 'pore.diameter')
        The dictionary key pointing to the pore property by which sphere
        diameters should be scaled
    Tsize : str (default = 'throat.diameter')
        The dictionary key pointing to the throat property by which cylinder
        diameters should be scaled
    Pcolor : str
        The dictionary key pointing to the pore property which will control
        the sphere colors.  The default is None, which results in a bright
        red for all pores.
    Tcolor : str
        The dictionary key pointing to the throat property which will control
        the cylinder colors.  The default is None, which results in a unform
        pale blue for all throats.
    cmap : str or Matplotlib colormap object (default is 'jet')
        The color map to use when converting pore and throat properties to
        RGB colors.  Can either be a string indicating which color map to
        fetch from matplotlib.cmap, or an actual cmap object.
    kwargs : dict
        Any additional kwargs that are received are passed to the VPython
        ``canvas`` object.  Default options are:

        *'height' = 500* - Height of canvas

        *'width' = 800* - Width of canvas

        *'background' = [0, 0, 0]* - Sets the background color of canvas

        *'ambient' = [0.2, 0.2, 0.3]* - Sets the brightness of lighting

    Returns
    -------
    canvas : VPython Canvas object
        The canvas object containing the generated scene. The object has
        several useful methods.

    Notes
    -----
    **Important**

    a) This does not work in Spyder. It should only be called from a
    Jupyter Notebook.

    b) This is only meant for relatively small networks. For proper
    visualization use Paraview.

    """
    import matplotlib.pyplot as plt

    try:
        from vpython import canvas, vec, sphere, cylinder
    except ModuleNotFoundError:
        raise Exception('VPython must be installed to use this function')

    if isinstance(cmap, str):
        cmap = getattr(plt.cm, cmap)

    if Pcolor is None:
        Pcolor = [vec(230/255, 57/255, 0/255)]*network.Np
    else:
        a = cmap(network[Pcolor]/network[Pcolor].max())
        Pcolor = [vec(row[0], row[1], row[2]) for row in a]

    if Tcolor is None:
        Tcolor = [vec(51/255, 153/255, 255/255)]*network.Nt
    else:
        a = cmap(network[Tcolor]/network[Tcolor].max())
        Tcolor = [vec(row[0], row[1], row[2]) for row in a]

    # Set default values for canvas properties
    if 'background' not in kwargs.keys():
        kwargs['background'] = vec(1.0, 1.0, 1.0)
    if 'height' not in kwargs.keys():
        kwargs['height'] = 500
    if 'width' not in kwargs.keys():
        kwargs['width'] = 800
    # Parse any given values for canvas properties
    for item in kwargs.keys():
        try:
            kwargs[item] = vec(*kwargs[item])
        except TypeError:
            pass
    scene = canvas(title=network.name, **kwargs)

    for p in network.Ps:
        r = network[Psize][p]/2
        xyz = network['pore.coords'][p]
        c = Pcolor[p]
        sphere(pos=vec(*xyz), radius=r, color=c,
               shininess=.5)

    for t in network.Ts:
        head = network['throat.endpoints.head'][t]
        tail = network['throat.endpoints.tail'][t]
        v = tail - head
        r = network[Tsize][t]
        L = np.sqrt(np.sum((head-tail)**2))
        c = Tcolor[t]
        cylinder(pos=vec(*head), axis=vec(*v), opacity=1, size=vec(L, r, r),
                 color=c)

    return scene


def plot_tutorial(network,
                  font_size=12,
                  line_width=2,
                  node_color='b',
                  edge_color='r',
                  node_size=500):  # pragma: no cover
    r"""
    Generate a network plot suitable for tutorials and explanations.

    Parameters
    ----------
    network : GenericNetwork
        The network to plot, should be 2D, since the z-coordinate will be
        ignored.
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
    import networkx as nx
    import matplotlib.pyplot as plt
    from openpnm.io import NetworkX

    G = NetworkX.to_networkx(network=network)
    pos = {i: network['pore.coords'][i, 0:2] for i in network.Ps}
    labels = {i: i for i in network.Ps}
    edge_labels = {tuple(network['throat.conns'][i, :]): i for i in network.Ts}

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
    xy_range = network.coords.ptp(axis=0)[dims]
    aspect_ratio = xy_range[0] / xy_range[1]
    fig.set_size_inches(5, 5 / aspect_ratio)

    return gplot
