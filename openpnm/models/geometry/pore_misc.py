import numpy as np


def neighbor(target, throat_prop='throat.seed', mode='min'):
    r"""
    Adopt a value from the values found in neighboring throats

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    throat_prop : string
        The dictionary key of the array containing the throat property to be
        used in the calculation.  The default is 'throat.seed'.

    mode : string
        Controls how the pore property is calculated.  Options are 'min',
        'max' and 'mean'.
    """
    network = target.project.network
    Ps = target.pores()
    data = target[throat_prop]
    neighborTs = network.find_neighbor_throats(pores=Ps,
                                               flatten=False,
                                               mode='intersection')
    values = np.ones((np.shape(Ps)[0],))*np.nan
    if mode == 'min':
        for pore in Ps:
            values[pore] = np.amin(data[neighborTs[pore]])
    if mode == 'max':
        for pore in Ps:
            values[pore] = np.amax(data[neighborTs[pore]])
    if mode == 'mean':
        for pore in Ps:
            values[pore] = np.mean(data[neighborTs[pore]])
    return values
