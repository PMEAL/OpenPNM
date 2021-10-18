import numpy as np
from openpnm.topotools.generators import cubic, tools


def cubic_template(template, spacing=1):
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
