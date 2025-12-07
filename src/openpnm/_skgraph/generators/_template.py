import numpy as np
from openpnm._skgraph.generators import cubic
from openpnm._skgraph.operations import trim_nodes


def cubic_template(template, spacing=1, connectivity=6,
                   node_prefix='node', edge_prefix='edge'):
    r"""
    Generate a simple cubic lattice matching the shape of the provided tempate

    Parameters
    ----------
    templte : ndarray
        Each ``True`` value will be treated as a vertex while all others
        will be trimmed.
    spacing : array_like or float
        The size of a unit cell in each direction. If an scalar is given it is
        applied in all 3 directions.

    Returns
    -------
    network : dict
        A dictionary containing 'node.coords' and 'edge.conns'

    """
    template = np.atleast_3d(template).astype(bool)
    # Generate a full cubic network
    temp = cubic(shape=template.shape, spacing=spacing,
                 connectivity=connectivity,
                 node_prefix=node_prefix, edge_prefix=edge_prefix)
    # Store some info about template
    coords = np.unravel_index(range(template.size), template.shape)
    coords = np.vstack(coords).T
    Np = coords.shape[0]
    temp[node_prefix+'.template_coords'] = coords
    temp[node_prefix+'.template_indices'] = np.arange(Np)
    # Trim pores not present in template
    temp = trim_nodes(network=temp, inds=~template.flatten())
    return temp


if __name__ == '__main__':
    im = np.ones([50, 50], dtype=bool)
    im[25:, ...] = False
    net = cubic_template(template=im)
    print(net.keys())
