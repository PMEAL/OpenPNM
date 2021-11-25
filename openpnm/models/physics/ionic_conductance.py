r"""
Pore-scale models for calculating ionic conductance of conduits.
"""
import numpy as _np
from openpnm.utils import logging
from openpnm.models.physics.utils import _poisson_conductance
logger = logging.getLogger(__name__)

__all__ = ["poisson", "electroneutrality"]


def poisson(target,
            pore_permittivity='pore.permittivity',
            throat_permittivity='throat.permittivity',
            size_factors='throat.diffusive_size_factors'):
    r"""
    Calculate the ionic conductance of conduits in network (using the
    Poisson equation for charge conservation), where a conduit is
    (1/2 pore - full throat - 1/2 pore). See the notes section.

    Parameters
    ----------
    target : GenericPhysics
        Physics object with which this model is associated.
    pore_permittivity : str
        Dictionary key of the pore permittivity values.
    throat_permittivity : str
        Dictionary key of the throat permittivity values.
    size_factors: str
        Dictionary key of the conduit size factors' values.

    Returns
    -------
    g : ndarray
        Array containing ionic conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    This function requires that all the necessary phase properties already
    be calculated.

    """
    epsilon0 = 8.854187817e-12
    tmp = _poisson_conductance(target=target,
                               pore_conductivity=pore_permittivity,
                               throat_conductivity=throat_permittivity,
                               size_factors=size_factors)
    return tmp * epsilon0


def electroneutrality(target,
                      pore_diffusivity='pore.diffusivity',
                      throat_diffusivity='throat.diffusivity',
                      size_factors="throat.diffusive_size_factors",
                      pore_volume='pore.volume',
                      pore_temperature='pore.temperature',
                      throat_temperature='throat.temperature',
                      pore_valence='pore.valence',
                      throat_valence='throat.valence',
                      pore_concentration='pore.concentration',
                      ions=[]):
    r"""
    Calculates the ionic conductance of conduits in network (assuming
    electroneutrality for charge conservation), where a conduit is
    (1/2 pore - full throat - 1/2 pore). See the notes section.

    Parameters
    ----------
    target : GenericPhysics
        Physics object with which this model is associated.
    pore_diffusivity : str
        Dictionary key of the pore diffusivity values
    throat_diffusivity : str
        Dictionary key of the throat diffusivity values
    size_factors: str
        Dictionary key of the conduit size factors' values.
    pore_volume : str
        Dictionary key of the pore volume values
    pore_temperature : str
        Dictionary key of the pore temperature values
    throat_temperature : str
        Dictionary key of the throat temperature values
    pore_valence : str
       Dictionary key of the pore ionic species valence values
    throat_valence : str
       Dictionary key of the throat ionic species valence values
    pore_concentration : str
       Dictionary key of the pore ionic species concentration values

    Returns
    -------
    g : ndarray
        Array containing ionic conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    This function requires that all the necessary phase properties already
    be calculated.

    """
    if ions == []:
        raise Exception('List of ions must be provided')

    # Fetch GenericPhysicss
    network = target.network
    domain = target._domain
    throats = domain.throats(target.name)
    phase = target.project.find_phase(target)
    cn = network.conns[throats]

    # Fetch model parameters
    R = 8.3145
    F = 96485.3329
    SF = network[size_factors]
    V1, V2 = network[pore_volume][cn].T
    T1, T2 = phase[pore_temperature][cn].T
    Tt = phase[throat_temperature][throats]
    # Pre-allocate g vectors
    g1, g2, gt = _np.zeros((3, len(cn)))

    # Add the contribution of each ion present to the overall conductance
    for ion in ions:
        try:
            c1, c2 = phase[f"{pore_concentration}.{ion}"][cn].T
        except KeyError:
            c1, c2 = _np.zeros(network.Np)[cn].T
        ct = (c1*V1 + c2*V2)/(V1 + V2)
        D1, D2 = phase[f"{pore_diffusivity}.{ion}"][cn].T
        Dt = phase[f"{throat_diffusivity}.{ion}"][throats]
        v1, v2 = phase[f"{pore_valence}.{ion}"][cn].T
        vt = phase[f"{throat_valence}.{ion}"][throats]
        # Add the contribute of the ion to g
        g1 += F**2 * v1**2 * D1*c1 / (R*T1)
        g2 += F**2 * v2**2 * D2*c2 / (R*T2)
        gt += F**2 * vt**2 * Dt*ct / (R*Tt)

    # Apply size factors and calculate the final conductance
    if isinstance(SF, dict):
        g1 *= SF[f"{size_factors}.pore1"][throats]
        gt *= SF[f"{size_factors}.throat"][throats]
        g2 *= SF[f"{size_factors}.pore2"][throats]
        return 1 / (1/g1 + 1/gt + 1/g2)
    return gt * SF
