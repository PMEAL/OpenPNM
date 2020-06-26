import networkx as nx
import matplotlib.pyplot as plt
from openpnm.io import NetworkX


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
