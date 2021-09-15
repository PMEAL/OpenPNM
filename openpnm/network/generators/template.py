import numpy as np
from openpnm.network.generators import cubic


def trim(network, pores=None, throats=None):
    if (pores is None) and (throats is None):
        raise Exception('Cannot trim pores and throats at the same time')
    if throats is not None:
        throats = np.atleast_1d(throats)
        keep = np.ones(network['throat.conns'].shape[0], dtype=bool)
        keep[throats] = False
        for item in network.keys():
            if item.startswith('throat'):
                network[item] = network[item][keep]
    elif pores is not None:
        pores = np.atleast_1d(pores)
        if pores.dtype == bool:
            pores = np.where(pores)[0]
        keep = np.ones(network['pore.coords'].shape[0], dtype=bool)
        keep[pores] = False
        for item in network.keys():
            if item.startswith('pore'):
                network[item] = network[item][keep]
        # Remove throats
        throats = np.any(np.isin(network['throat.conns'], pores), axis=1)
        network = trim(network, throats=throats)
        # Renumber throat conns
        remapping = np.cumsum(keep) - 1
        network['throat.conns'] = remapping[network['throat.conns']]
    return network


def cubic_template(template, spacing=1):
    template = np.atleast_3d(template)
    # Generate a full cubic network
    temp = cubic(shape=template.shape, spacing=spacing)
    # Store some info about template
    coords = np.unravel_index(range(template.size), template.shape)
    temp['pore.template_coords'] = np.vstack(coords).T
    temp['pore.template_indices'] = np.arange(Np)
    # Trim pores not present in template
    temp = trim(network=temp, pores=template.flatten())
    return temp
