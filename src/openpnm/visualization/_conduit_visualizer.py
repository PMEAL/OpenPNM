import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle, Polygon, Rectangle

from openpnm._skgraph.tools import rotate_coords
from openpnm.models.geometry import throat_endpoints, throat_length
from openpnm.network import Network
from openpnm.topotools import find_connected_sites

__all__ = [
    'draw_conduit',
]


def draw_conduit(network, throat):
    """Draws a subset of a network given throat numbers."""
    pn = network
    P1, P2 = find_connected_sites(g=pn, bonds=throat)
    new_net = Network(coords=pn.coords[[P1, P2], :], conns=np.atleast_2d([0, 1]))
    new_net.regenerate_models()
    new_net['pore.diameter'] = pn['pore.diameter'][[P1, P2]]
    new_net['throat.diameter'] = [pn['throat.diameter'][throat]]
    new_net['throat.length'] = throat_length.circles_and_rectangles(new_net)
    new_net['throat.endpoints'] = throat_endpoints.spheres_and_cylinders(new_net)

    Pcrds = np.vstack((pn.coords[0, :], pn.coords[0, :]))
    Pcrds[1, :] = Pcrds[1, :] + np.array([new_net['throat.spacing'], 0, 0])
    fig, ax = plt.subplots(figsize=[10, 5])
    patches = []
    for i, P in enumerate([0, 1]):
        x, y, r = Pcrds[P, 0], Pcrds[P, 1], new_net['pore.diameter'][P]/2
        circle = Circle((x, y), r)
        patches.append(circle)
    Tcrds = np.vstack(list(new_net['throat.endpoints'].values()))
    Tcrds = rotate_coords(Tcrds, b=270)
    H = new_net['throat.diameter'][0]
    W = new_net['throat.length'][0]
    R1 = new_net['pore.diameter'][0]/2
    R2 = new_net['pore.diameter'][1]/2
    rect = Rectangle(xy=(Tcrds[0, 0], Tcrds[0, 1]-H/2), height=H, width=W)
    patches.append(rect)
    p = PatchCollection(patches, alpha=0.5, edgecolor='k', linewidth=3)
    p.set_array([1, 1, 2])
    p.cmap = plt.cm.bwr
    ax.add_collection(p)
    ax.scatter(*Pcrds[0, :2], marker='+', c='k', s=100)
    ax.scatter(*Pcrds[1, :2], marker='+', c='k', s=100)
    left = min(Pcrds[:, 0])-R1
    right = max(Pcrds[:, 0])+R2
    ax.set_xlim([left - abs(left-right)/32, right + abs(left-right)/32])
    temp = (Pcrds[:, 1]).mean()
    ax.set_ylim([temp - abs(left-right)/4, temp + abs(left-right)/4])
    ax.annotate(text='',
                xy=Pcrds[0, :2],
                xytext=(Pcrds[0, 0], Pcrds[0, 1]-R1),
                arrowprops=dict(arrowstyle='<-', lw=2))
    ax.annotate(text=f"{R1}",
                xy=Pcrds[0, :2],
                xytext=(Pcrds[0, 0], Pcrds[0, 1]-R1/2),
                textcoords='data')
    return ax
