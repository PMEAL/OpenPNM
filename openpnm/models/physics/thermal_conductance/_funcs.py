from openpnm.models.physics._utils import _poisson_conductance
from openpnm.models import _doctxt


__all__ = ["generic_thermal", "series_resistors"]


@_doctxt
def generic_thermal(target,
                    pore_conductivity='pore.thermal_conductivity',
                    throat_conductivity='throat.thermal_conductivity',
                    size_factors='throat.diffusive_size_factors'):
    r"""
    Calculate the thermal conductance of conduits in network.

    Parameters
    ----------
    %(target_blurb)s
    pore_conductivity : str
        %(dict_blurb)s thermal conductivity
    throat_conductivity : str
        %(dict_blurb)s throat thermal conductivity
    size_factors : str
        %(dict_blurb)s conduit diffusive size factors

    Returns
    -------
    %(return_arr)s thermal conductance

    """
    return _poisson_conductance(target=target,
                                pore_conductivity=pore_conductivity,
                                throat_conductivity=throat_conductivity,
                                size_factors=size_factors)


def series_resistors(target,
                     pore_thermal_conductivity='pore.thermal_conductivity',
                     throat_thermal_conductivity='throat.thermal_conductivity',
                     size_factors='throat.diffusive_size_factors'):
    r"""
    Calculate the thermal conductance of conduits in network.

    Parameters
    ----------
    %(target_blurb)s
    pore_conductivity : str
        %(dict_blurb)s thermal conductivity
    throat_conductivity : str
        %(dict_blurb)s throat thermal conductivity
    size_factors : str
        %(dict_blurb)s conduit diffusive size factors

    Returns
    -------
    %(return_arr)s thermal conductance

    """
    return _poisson_conductance(target=target,
                                pore_conductivity=pore_thermal_conductivity,
                                throat_conductivity=throat_thermal_conductivity,
                                size_factors=size_factors)
