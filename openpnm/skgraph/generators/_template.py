import numpy as np
from skgraph.generators import cubic, tools


def cubic_template(template, spacing=1):
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
        A dictionary containing 'vert.coords' and 'edge.conns'

    """
    template = np.atleast_3d(template).astype(bool)
    # Generate a full cubic network
    temp = cubic(shape=template.shape, spacing=spacing)
    # Store some info about template
    coords = np.unravel_index(range(template.size), template.shape)
    coords = np.vstack(coords).T
    Np = coords.shape[0]
    temp['vert.template_coords'] = coords
    temp['vert.template_indices'] = np.arange(Np)
    # Trim pores not present in template
    temp = tools.trim(network=temp, vert_ids=~template.flatten())
    return temp


if __name__ == '__main__':
    im = np.ones([50, 50], dtype=bool)
    im[25:, ...] = False
    net = cubic_template(template=im)
    print(net.keys())
