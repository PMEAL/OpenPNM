r"""

.. autofunction:: openpnm.models.misc.neighbor_lookups.from_neighbor_throats
.. autofunction:: openpnm.models.misc.neighbor_lookups.from_neighbor_pores

"""

import numpy as np
from openpnm.utils import logging
logger = logging.getLogger(__name__)


def from_neighbor_throats(target, throat_prop='throat.seed', mode='min'):
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

    Returns
    -------
    value : NumPy ndarray
        Array containing customized values based on those of adjacent throats.

    """
    prj = target.project
    network = prj.network
    lookup = prj.find_full_domain(target)
    Ps = lookup.map_pores(target.pores(), target)
    data = lookup[throat_prop]
    neighborTs = network.find_neighbor_throats(pores=Ps,
                                               flatten=False,
                                               mode='or')
    values = np.ones((np.shape(Ps)[0],))*np.nan
    if mode == 'min':
        for pore in range(len(Ps)):
            values[pore] = np.amin(data[neighborTs[pore]])
    if mode == 'max':
        for pore in range(len(Ps)):
            values[pore] = np.amax(data[neighborTs[pore]])
    if mode == 'mean':
        for pore in range(len(Ps)):
            values[pore] = np.mean(data[neighborTs[pore]])
    return values


def from_neighbor_pores(target, pore_prop='pore.seed', mode='min'):
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

    Returns
    -------
    value : NumPy ndarray
        Array containing customized values based on those of adjacent pores.

    """
    prj = target.project
    network = prj.network
    throats = network.map_throats(target.throats(), target)
    P12 = network.find_connected_pores(throats)
    lookup = prj.find_full_domain(target)
    pvalues = lookup[pore_prop][P12]
    if mode == 'min':
        value = np.amin(pvalues, axis=1)
    if mode == 'max':
        value = np.amax(pvalues, axis=1)
    if mode == 'mean':
        value = np.mean(pvalues, axis=1)
    return value
