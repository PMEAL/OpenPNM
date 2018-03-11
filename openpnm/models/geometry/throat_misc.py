import np as np


def neighbor(target, pore_prop='pore.seed', mode='min'):
    r"""
    Adopt a value based on the values in neighboring pores

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_prop : string
        The dictionary key to the array containing the pore property to be
        used in the calculation.  Default is 'pore.seed'.

    mode : string
        Controls how the throat property is calculated.  Options are 'min',
        'max' and 'mean'.

    """
    network = target.project.network
    throats = network.throats(target.name)
    P12 = network.find_connected_pores(throats)
    pvalues = network[pore_prop][P12]
    if mode == 'min':
        value = np.amin(pvalues, axis=1)
    if mode == 'max':
        value = np.amax(pvalues, axis=1)
    if mode == 'mean':
        value = np.mean(pvalues, axis=1)
    return value
