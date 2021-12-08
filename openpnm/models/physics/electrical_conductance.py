from openpnm.models.physics._utils import _poisson_conductance
from openpnm.models import _doctxt


__all__ = ["generic_electrical", "series_resistors"]


@_doctxt
def generic_electrical(
    target,
    pore_conductivity='pore.electrical_conductivity',
    throat_conductivity='throat.electrical_conductivity',
    size_factors='throat.diffusive_size_factors'
):
    r"""
    Calculate the electrical conductance of conduits in network, where a
    conduit is ( 1/2 pore - full throat - 1/2 pore ). See the notes section.

    Parameters
    ----------
    %(target_blurb)s
    pore_conductivity : str
        %(dict_blurb)s pore electrical conductivity
    throat_conductivity : str
        %(dict_blurb)s throat electrical conductivity
    size_factors: str
        %(dict_blurb)s conduit diffusive size factors

    Returns
    -------
    %(return_arr)s electrical conductance

    """
    return _poisson_conductance(target=target,
                                pore_conductivity=pore_conductivity,
                                throat_conductivity=throat_conductivity,
                                size_factors=size_factors)


def series_resistors(
    target,
    pore_conductivity='pore.electrical_conductivity',
    throat_conductivity='throat.electrical_conductivity',
    size_factors='throat.diffusive_size_factors'
):
    r"""
    Calculate the electrical conductance of conduits in network, where a
    conduit is ( 1/2 pore - full throat - 1/2 pore ). See the notes section.

    Parameters
    ----------
    %(target_blurb)s
    pore_conductivity : str
        %(dict_blurb)s pore electrical conductivity
    throat_conductivity : str
        %(dict_blurb)s throat electrical conductivity
    size_factors: str
        %(dict_blurb)s conduit diffusive size factors

    Returns
    -------
    %(return_arr)s electrical conductance

    """
    return _poisson_conductance(target=target,
                                pore_conductivity=pore_conductivity,
                                throat_conductivity=throat_conductivity,
                                size_factors=size_factors)
