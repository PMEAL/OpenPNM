import numpy as np
import networkx as nx
import matplotlib.pyplot as plt


def plot_connections(network, throats=None, fig=None, size_by=None,
                     color_by=None, cmap='jet', color='b', alpha=1.0,
                     linestyle='solid', linewidth=1, **kwargs):
    r"""
    Produce a 3D plot of the network topology.

    This shows how throats connect for quick visualization without having
    to export data to veiw in Paraview.

    Parameters
    ----------
    network : OpenPNM Network Object
        The network whose topological connections to plot
    throats : array_like (optional)
        The list of throats to plot if only a sub-sample is desired.  This is
        useful for inspecting a small region of the network.  If no throats are
        specified then all throats are shown.
    fig : Matplotlib figure handle and line property arguments
        If a ``fig`` is supplied, then the topology will be overlaid on this
        plot.  This makes it possible to combine coordinates and connections,
        and to color throats differently
    size_by : array_like
        An ND-array of throat values (e.g. alg['throat.rate']).  These
        values are normalized by scaled by ``markersize``.
    color_by : str or array_like
        An ND-array of throat values (e.g. alg['throat.rate']).
    cmap : str or cmap object
        The matplotlib colormap to use if specfying a throat property
        for ``color_by``
    color : str
        A matplotlib named color (e.g. 'r' for red).
    alpha : float
        The transparency of the lines, with 1 being solid and 0 being invisible
    linestyle : str
        Can be one of {'solid', 'dashed', 'dashdot', 'dotted'}
    linewidth : float
        Controls the thickness of drawn lines.  Is used to scale the thickness
        if ``size_by`` is given.

    Notes
    -----
    The figure handle returned by this method can be passed into
    ``plot_coordinates`` to create a plot that combines pore coordinates and
    throat connections, and vice versa.

    See Also
    --------
    plot_coordinates

    Examples
    --------
    >>> import openpnm as op
    >>> import matplotlib as mpl
    >>> mpl.use('Agg')
    >>> pn = op.network.Cubic(shape=[10, 10, 3])
    >>> pn.add_boundary_pores()
    >>> Ts = pn.throats('*boundary', mode='nor')
    >>> # Create figure showing boundary throats
    >>> fig = op.topotools.plot_connections(network=pn, throats=Ts)
    >>> Ts = pn.throats('*boundary')
    >>> # Pass existing fig back into function to plot additional throats
    >>> fig = op.topotools.plot_connections(network=pn, throats=Ts,
    ...                                     fig=fig, colors='r')

    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib.collections import LineCollection
    from mpl_toolkits.mplot3d.art3d import Line3DCollection
    from matplotlib import colors as mcolors
    from matplotlib import cm
    from openpnm.topotools import dimensionality

    Ts = network.Ts if throats is None else network._parse_indices(throats)
    dim = dimensionality(network)
    ThreeD = True if dim.sum() == 3 else False
    # Add a dummy axis for 1D networks
    if dim.sum() == 1:
        dim[np.argwhere(dim == False)[0]] = True

    fig = plt.figure() if fig is None else fig
    ax = fig.gca()
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
        if not size_by.startswith('throat.'):
            size_by = 'throat.' + size_by
        linewidth = size_by / size_by.max() * linewidth

    if ThreeD:
        lc = Line3DCollection(throat_pos, colors=color, cmap=cmap,
                              linestyles=linestyle, linewidths=linewidth,
                              antialiaseds=np.ones_like(network.Ts))
    else:
        lc = LineCollection(throat_pos, colors=color, cmap=cmap,
                            linestyles=linestyle, linewidths=linewidth,
                            antialiaseds=np.ones_like(network.Ts))
    ax.add_collection(lc)

    _scale_3d_axes(ax=ax, X=X, Y=Y, Z=Z)
    _label_axes(ax=ax, X=X, Y=Y, Z=Z)

    return fig


def plot_coordinates(network, pores=None, fig=None, size_by=None,
                     color_by=None, cmap='jet', color='r', alpha=1.0,
                     marker='o', markersize=1, **kwargs):
    r"""
    Produce a 3D plot showing specified pore coordinates as markers.

    Parameters
    ----------
    network : OpenPNM Network Object
        The network whose topological connections to plot
    pores : array_like (optional)
        The list of pores to plot if only a sub-sample is desired.  This is
        useful for inspecting a small region of the network.  If no pores are
        specified then all are shown.
    fig : Matplotlib figure handle
        If a ``fig`` is supplied, then the coordinates will be overlaid.  This
        enables the plotting of multiple different sets of pores as well as
        throat connections from ``plot_connections``.
    size_by : str or array_like
        An ND-array of pore values (e.g. alg['pore.concentration']).  These
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

    Notes
    -----
    The figure handle returned by this method can be passed into
    ``plot_connections`` to create a plot that combines pore coordinates and
    throat connections, and vice versa.

    See Also
    --------
    plot_connections

    Examples
    --------
    >>> import openpnm as op
    >>> import matplotlib as mpl
    >>> mpl.use('Agg')
    >>> pn = op.network.Cubic(shape=[10, 10, 3])
    >>> pn.add_boundary_pores()
    >>> Ps = pn.pores('internal')
    >>> # Create figure showing internal pores
    >>> fig = op.topotools.plot_coordinates(pn, pores=Ps, c='b')
    >>> Ps = pn.pores('*boundary')
    >>> # Pass existing fig back into function to plot boundary pores
    >>> fig = op.topotools.plot_coordinates(pn, pores=Ps, fig=fig, c='r')

    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    from openpnm.topotools import dimensionality

    Ps = network.Ps if pores is None else network._parse_indices(pores)

    dim = dimensionality(network)
    ThreeD = True if dim.sum() == 3 else False
    # Add a dummy axis for 1D networks
    if dim.sum() == 1:
        dim[np.argwhere(dim == False)[0]] = True
    # Add 2 dummy axes for 0D networks (1 pore only)
    if dim.sum() == 0:
        dim[[0, 1]] = True

    fig = plt.figure() if fig is None else fig
    ax = fig.gca()
    if ThreeD and ax.name != '3d':
        fig.delaxes(ax)
        ax = fig.add_subplot(111, projection='3d')

    # Collect specified coordinates
    X, Y, Z = network['pore.coords'][Ps].T

    # Parse formating args
    if 'c' in kwargs.keys():
        color = kwargs.pop('c')
    if 's' in kwargs.keys():
        markersize = kwargs.pop('s')
    if color_by is not None:
        color = cm.get_cmap(name=cmap)(color_by / color_by.max())
    if size_by is not None:
        markersize = size_by / size_by.max() * markersize

    if ThreeD:
        ax.scatter(X, Y, Z, c=color, s=markersize,
                   marker=marker, alpha=alpha)
        _scale_3d_axes(ax=ax, X=X, Y=Y, Z=Z)
    else:
        X_temp, Y_temp = np.column_stack((X, Y, Z))[:, dim].T
        ax.scatter(X_temp, Y_temp, c=color, s=markersize,
                   marker=marker, alpha=alpha)
        _scale_3d_axes(ax=ax, X=X, Y=Y, Z=np.zeros_like(Y))

    _label_axes(ax=ax, X=X, Y=Y, Z=Z)

    return fig


def _label_axes(ax, X, Y, Z):
    labels = ["X", "Y", "Z"]
    dim = np.zeros(3, dtype=bool)
    for i, arr in enumerate([X, Y, Z]):
        if np.unique(arr).size > 1:
            dim[i] = True
    # Add a dummy axis for 1D networks
    if dim.sum() == 1:
        dim[np.argwhere(dim == False)[0]] = True
    # Add 2 dummy axes for 0D networks (1 pore only)
    if dim.sum() == 0:
        dim[[0, 1]] = True
    dim_idx = np.argwhere(dim == True).squeeze()
    ax.set_xlabel(labels[dim_idx[0]])
    ax.set_ylabel(labels[dim_idx[1]])
    if hasattr(ax, "set_zlim"):
        ax.set_zlabel("Z")


def _scale_3d_axes(ax, X, Y, Z):
    if not hasattr(ax, '_scaled'):
        ax._scaled = True
        if not hasattr(ax, "set_zlim"):
            ax.axis("equal")
        max_range = np.ptp([X, Y, Z]).max() / 2
        mid_x = (X.max() + X.min()) * 0.5
        mid_y = (Y.max() + Y.min()) * 0.5
        mid_z = (Z.max() + Z.min()) * 0.5
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        if hasattr(ax, "set_zlim"):
            ax.set_zlim(mid_z - max_range, mid_z + max_range)


def plot_networkx(network, plot_throats=True, labels=None, colors=None,
                  scale=1, ax=None, alpha=1.0):
    r"""
    Create a pretty 2d plot for 2d OpenPNM networks.

    Parameters
    ----------
    network : OpenPNM Network object

    plot_throats : boolean
        Plots throats as well as pores, if True.

    labels : list
        List of OpenPNM labels

    colors : list
        List of corresponding colors to the given `labels`.

    scale : float
        Scale factor for size of pores.
    """
    from networkx import Graph, draw_networkx_nodes, draw_networkx_edges
    from matplotlib.collections import PathCollection
    import matplotlib.pyplot as plt
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

    G = Graph()
    pos = {network.Ps[i]: [x[i], y[i]] for i in range(network.Np)}
    try:
        node_size = scale * network['pore.diameter']
    except KeyError:
        node_size = np.ones_like(x) * scale * 0.5
    if not np.isfinite(node_size).all():
        raise Exception('nan/inf values found in network["pore.diameter"]')
    node_color = np.array(['k'] * len(network.Ps))

    if labels:
        if type(labels) is not list:
            labels = [labels]
        if type(colors) is not list:
            colors = [colors]
        if len(labels) != len(colors):
            raise('len(colors) must be equal to len(labels)!')
        for label, color in zip(labels, colors):
            node_color[network.pores(label)] = color

    if ax is None:
        fig, ax = plt.subplots()
    ax.set_aspect('equal', adjustable='box')
    offset = node_size.max() * 0.25
    ax.set_xlim((x.min() - offset, x.max() + offset))
    ax.set_ylim((y.min() - offset, y.max() + offset))
    ax.axis("off")

    # Keep track of already plotted nodes
    temp = [id(item) for item in ax.collections if type(item) == PathCollection]

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
        if type(item) == PathCollection and id(item) not in temp:
            item.set_sizes(markersize)

    return gplot


def plot_vpython(network,
                 Psize='pore.diameter',
                 Tsize='throat.diameter',
                 Pcolor=None,
                 Tcolor=None,
                 cmap='jet', **kwargs):
    r"""
    Quickly visualize a network in 3D using VPython.

    Parameters
    ----------
    network : OpenPNM Network Object
        The network to visualize
    Psize : string (default = 'pore.diameter')
        The dictionary key pointing to the pore property by which sphere
        diameters should be scaled
    Tsize : string (default = 'throat.diameter')
        The dictionary key pointing to the throat property by which cylinder
        diameters should be scaled
    Pcolor : string
        The dictionary key pointing to the pore property which will control
        the sphere colors.  The default is None, which results in a bright
        red for all pores.
    Tcolor : string
        The dictionary key pointing to the throat property which will control
        the cylinder colors.  The default is None, which results in a unform
        pale blue for all throats.
    cmap : string or Matplotlib colormap object (default is 'jet')
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

    a) This does not work in Spyder.  It should only be called from a Jupyter
    Notebook.

    b) This is only meant for relatively small networks.  For proper
    visualization use Paraview.

    """
    import matplotlib.pyplot as plt

    try:
        from vpython import canvas, vec, sphere, cylinder
    except ModuleNotFoundError:
        raise Exception('VPython must be installed to use this function')

    if type(cmap) == str:
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


def plot_tutorial(network, font_size=24, line_width=3,
                  node_color='b', edge_color='r', node_size=2000):
    r"""
    Generate a network plot suitable for tutorials and explanations.

    Parameters
    ----------
    network : OpenPNM Network object
        The network to plot, should be 2D, since the z-coordinate will be
        ignored.
    font_size : int
        Size of font to use for labels
    line_width : int
        Thickness of edge lines and node borders
    node_color : str
        Color of node border
    edge_color : str
        Color of edge lines
    node_size : int
        Size of node circle

    Returns
    -------
    g : NetworkX plot object

    """
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
    nx.draw_networkx_edges(G, pos,
                           width=line_width,
                           edge_color=edge_color)
    nx.draw_networkx_labels(G, pos,
                            labels=labels,
                            font_size=font_size,
                            font_color='k')
    nx.draw_networkx_edge_labels(G, pos,
                                 edge_labels=edge_labels,
                                 font_size=font_size,
                                 font_color='k')
    plt.axis('off')
    return gplot
