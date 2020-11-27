r"""
"""
import numpy as np
from openpnm.utils import logging
logger = logging.getLogger(__name__)


def from_neighbor_throats(target, prop, mode='min', ignore_nans=True):
    r"""
    Adopt a value from the values found in neighboring throats

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.
    prop : string
        The dictionary key of the array containing the throat property to be
        used in the calculation.
    mode : string
        Controls how the pore property is calculated.  Options are 'min',
        'max' and 'mean'.

    Returns
    -------
    value : ND-array
        Array containing customized values based on those of adjacent throats.

    """
    prj = target.project
    network = prj.network
    boss = prj.find_full_domain(target)
    data = boss[prop]
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
    Ps = boss.map_pores(target.pores(), target)
    return np.array(values)[Ps]


def from_neighbor_pores(target, prop, mode='min', ignore_nans=True):
    r"""
    Adopt a value based on the values in neighboring pores

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.
    prop : string
        The dictionary key to the array containing the pore property to be
        used in the calculation.
    mode : string
        Controls how the throat property is calculated.  Options are 'min',
        'max' and 'mean'.
    ignore_nans : boolean (default is ``True``)
        If ``True`` the result will ignore ``nans`` in the neighbors

    Returns
    -------
    value : ND-array
        Array containing customized values based on those of adjacent pores.

    """
    prj = target.project
    lookup = prj.find_full_domain(target)
    network = prj.network
    throats = network.map_throats(target.throats(), target)
    P12 = network.find_connected_pores(throats)
    pvalues = lookup[prop][P12]
    if ignore_nans:
        pvalues = np.ma.MaskedArray(data=pvalues, mask=np.isnan(pvalues))
    try:  # If pvalues is not empty
        if mode == 'min':
            value = np.amin(pvalues, axis=1)
        if mode == 'max':
            value = np.amax(pvalues, axis=1)
        if mode == 'mean':
            value = np.mean(pvalues, axis=1)
    except np.AxisError:  # Handle case of empty pvalues
        value = []
    return np.array(value)
