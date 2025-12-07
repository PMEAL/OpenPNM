"""
Neighbor Lookups
================

"""
import logging
import numpy as np
logger = logging.getLogger(__name__)


__all__ = ['from_neighbor_throats', 'from_neighbor_pores']


def from_neighbor_throats(target, prop, mode='min', ignore_nans=True):
    r"""
    Adopt a value from the values found in neighboring throats

    Parameters
    ----------
    target : Base
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.
    prop : str
        The dictionary key of the array containing the throat property to be
        used in the calculation.
    mode : str
        Controls how the pore property is calculated. The default value is
        'min'. Options are:

            ===========  =====================================================
            mode         meaning
            ===========  =====================================================
            'min'        Returns the value of the minimum property of the
                         neighboring throats
            'max'        Returns the value of the maximum property of the
                         neighboring throats
            'mean'       Returns the value of the mean property of the
                         neighboring throats
            'sum'        Returns the sum of the property of the neighboring
                         throats
            ===========  =====================================================

    Returns
    -------
    value : ndarray
        Array containing customized values based on those of adjacent throats.

    """
    network = target.network
    data = target[prop]
    nans = np.isnan(data)
    im = network.create_incidence_matrix()
    if mode == 'min':
        if ignore_nans:
            data[nans] = np.inf
        values = np.ones((network.Np, ))*np.inf
        np.minimum.at(values, im.row, data[im.col])
    if mode == 'max':
        if ignore_nans:
            data[nans] = -np.inf
        values = np.ones((network.Np, ))*-np.inf
        np.maximum.at(values, im.row, data[im.col])
    if mode == 'mean':
        if ignore_nans:
            data[nans] = 0
        values = np.zeros((network.Np, ))
        np.add.at(values, im.row, data[im.col])
        counts = np.zeros((network.Np, ))
        np.add.at(counts, im.row, np.ones((network.Nt, ))[im.col])
        if ignore_nans:
            np.subtract.at(counts, im.row, nans[im.col])
        values = values/counts
    if mode == 'sum':
        if ignore_nans:
            data[nans] = 0
        values = np.zeros((network.Np, ))
        np.add.at(values, im.row, data[im.col])
    return values


def from_neighbor_pores(target, prop, mode='min', ignore_nans=True):
    r"""
    Adopt a value based on the values in neighboring pores

    Parameters
    ----------
    target : Base
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.
    prop : str
        The dictionary key to the array containing the pore property to be
        used in the calculation.
    mode : str
        Controls how the pore property is calculated. The default value is
        'min'. Options are:

            ===========  =====================================================
            mode         meaning
            ===========  =====================================================
            'min'        Returns the value of the minimum property of the
                         neighboring pores
            'max'        Returns the value of the maximum property of the
                         neighboring pores
            'mean'       Returns the value of the mean property of the
                         neighboring pores
            'sum'        Returns the sum of the property of the neighrboring
                         pores
            ===========  =====================================================

    ignore_nans : bool (default is ``True``)
        If ``True`` the result will ignore ``nans`` in the neighbors

    Returns
    -------
    value : ndarray
        Array containing customized values based on those of adjacent pores.

    """
    network = target.network
    throats = target.Ts
    P12 = network.find_connected_pores(throats)
    pvalues = target[prop][P12]
    if ignore_nans:
        pvalues = np.ma.MaskedArray(data=pvalues, mask=np.isnan(pvalues))
    try:  # If pvalues is not empty
        if mode == 'min':
            value = np.amin(pvalues, axis=1)
        if mode == 'max':
            value = np.amax(pvalues, axis=1)
        if mode == 'mean':
            value = np.mean(pvalues, axis=1)
        if mode == 'sum':
            value = np.sum(pvalues, axis=1)
    except np.AxisError:  # Handle case of empty pvalues
        value = []
    return np.array(value)
