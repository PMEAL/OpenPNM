import numpy as np
from openpnm.network.generators import cubic, tools


def cubic_template(template, spacing=1):
    template = np.atleast_3d(template)
    # Generate a full cubic network
    temp = cubic(shape=template.shape, spacing=spacing)
    # Store some info about template
    coords = np.unravel_index(range(template.size), template.shape)
    temp['pore.template_coords'] = np.vstack(coords).T
    temp['pore.template_indices'] = np.arange(Np)
    # Trim pores not present in template
    temp = tools.trim(network=temp, pores=template.flatten())
    return temp
