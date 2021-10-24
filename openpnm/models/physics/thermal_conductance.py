r"""
Pore-scale models for calculating the thermal conductance of conduits.
"""
from openpnm.models.physics.utils import _poisson_conductance

__all__ = ["generic_thermal", "series_resistors"]


def generic_thermal(target,
                    pore_conductivity='pore.thermal_conductivity',
                    throat_conductivity='throat.thermal_conductivity',
                    size_factors='throat.diffusive_size_factors'):
    r"""
    Calculate the thermal conductance of conduits in network, where a
    conduit is (1/2 pore - full throat - 1/2 pore). See the notes section.

    Parameters
    ----------
    target : GenericPhysics
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.
    pore_conductivity : str
        Dictionary key of the pore thermal conductivity values
    throat_conductivity : str
        Dictionary key of the throat thermal conductivity values
    size_factors: str
        Dictionary key of the conduit diffusive size factors' values.

    Returns
    -------
    g : ndarray
        Array containing thermal conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    This function requires that all the necessary phase properties already
    be calculated.

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
    Calculate the thermal conductance of conduits in network, where a
    conduit is (1/2 pore - full throat - 1/2 pore). See the notes section.

    Parameters
    ----------
    target : GenericPhysics
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.
    pore_conductivity : str
        Dictionary key of the pore thermal conductivity values
    throat_conductivity : str
        Dictionary key of the throat thermal conductivity values
    size_factors: str
        Dictionary key of the conduit diffusive size factors' values.

    Returns
    -------
    g : ndarray
        Array containing thermal conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    This function requires that all the necessary phase properties already
    be calculated.

    """
    return _poisson_conductance(target=target,
                                pore_conductivity=pore_thermal_conductivity,
                                throat_conductivity=throat_thermal_conductivity,
                                size_factors=size_factors)
