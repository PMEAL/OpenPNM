r"""
===============================================================================
pore_misc -- miscillaneous and generic functions to apply to pores
===============================================================================

"""
import scipy as _sp
from OpenPNM.Base import logging as _logging
_logger = _logging.getLogger(__name__)


def constant(geometry, value, **kwargs):
    r"""
    Assign specified constant value to all pores on the specifed Geometry
    object.

    Parameters
    ----------
    geometry : OpenPNM Geometry object
        The Geometry object to which this model will apply.  This is necessary
        to determine the length of the array to generate.

    value : scalar
        The scalar value to be applied

    Notes
    -----
    This function is redundant and could be accomplished with
    geometry['pore.prop'] = value.  This model was added to the code before
    this type of assignment was possible, and is only kept for backwards
    compatibility.
    """
    _logger.warning('This pore-scale model is deprecated, use direct ' +
                    'assignment of scalars to desired propname instead.')
    if _sp.shape(value):
        raise Exception('Only scalar values can be assigned')
    value = _sp.ones(geometry.num_pores(),)*value
    return value


def random(geometry, seed=None, num_range=[0, 1], **kwargs):
    r"""
    Assign random number to pore bodies
    note: should this be called 'poisson'?
    """
    range_size = num_range[1]-num_range[0]
    range_min = num_range[0]
    _sp.random.seed(seed=seed)
    value = _sp.random.rand(geometry.num_pores(),)
    value = value*range_size + range_min
    return value


def neighbor(geometry, throat_prop='throat.seed', mode='min', **kwargs):
    r"""
    Adopt a value from the values found in neighboring throats

    Parameters
    ----------
    geometry : OpenPNM Geometry Object
        The object containing the ``throat_prop`` to be used.

    throat_prop : string
        The dictionary key of the array containing the throat property to be
        used in the calculation.  The default is 'throat.seed'.

    mode : string
        Controls how the pore property is calculated.  Options are 'min',
        'max' and 'mean'.
    """
    network = geometry._net
    Ps = geometry.pores()
    data = geometry[throat_prop]
    neighborTs = network.find_neighbor_throats(pores=Ps,
                                               flatten=False,
                                               mode='intersection')
    values = _sp.ones((_sp.shape(Ps)[0],))*_sp.nan
    if mode == 'min':
        for pore in Ps:
            values[pore] = _sp.amin(data[neighborTs[pore]])
    if mode == 'max':
        for pore in Ps:
            values[pore] = _sp.amax(data[neighborTs[pore]])
    if mode == 'mean':
        for pore in Ps:
            values[pore] = _sp.mean(data[neighborTs[pore]])
    return values
