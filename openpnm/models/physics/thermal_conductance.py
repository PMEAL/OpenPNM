from openpnm.models.physics._utils import _poisson_conductance
from openpnm.utils import Docorator


docstr = Docorator()
__all__ = ["generic_thermal", "series_resistors"]


@docstr.get_sections(base=['models.physics.thermal_conductance'],
                     sections=['Returns'])
@docstr.dedent
def generic_thermal(target,
                    pore_conductivity='pore.thermal_conductivity',
                    throat_conductivity='throat.thermal_conductivity',
                    size_factors='throat.diffusive_size_factors'):
    r"""
    Calculate the thermal conductance of conduits in network, where a
    conduit is (1/2 pore - full throat - 1/2 pore). See the notes section.

    Parameters
    ----------
    %(models.target.parameters)s
    pore_conductivity : str
        Dictionary key of the pore thermal conductivity values
    throat_conductivity : str
        Dictionary key of the throat thermal conductivity values
    size_factors: str
        Dictionary key of the conduit diffusive size factors' values.

    Returns
    -------
    g_therm : ndarray
        A numpy ndarray containing thermal conductance values

    Notes
    -----
    This function requires that all the necessary phase properties already
    be calculated.

    """
    return _poisson_conductance(target=target,
                                pore_conductivity=pore_conductivity,
                                throat_conductivity=throat_conductivity,
                                size_factors=size_factors)


@docstr.dedent
def series_resistors(target,
                     pore_thermal_conductivity='pore.thermal_conductivity',
                     throat_thermal_conductivity='throat.thermal_conductivity',
                     size_factors='throat.diffusive_size_factors'):
    r"""
    Calculate the thermal conductance of conduits in network, where a
    conduit is (1/2 pore - full throat - 1/2 pore). See the notes section.

    Parameters
    ----------
    %(models.target.parameters)s
    pore_conductivity : str
        Dictionary key of the pore thermal conductivity values
    throat_conductivity : str
        Dictionary key of the throat thermal conductivity values
    size_factors: str
        Dictionary key of the conduit diffusive size factors' values.

    Returns
    -------
    %(models.physics.thermal_conductance.returns)s

    Notes
    -----
    This function requires that all the necessary phase properties already
    be calculated.

    """
    return _poisson_conductance(target=target,
                                pore_conductivity=pore_thermal_conductivity,
                                throat_conductivity=throat_thermal_conductivity,
                                size_factors=size_factors)
