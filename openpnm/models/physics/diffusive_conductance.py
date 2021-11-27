import numpy as _np
import scipy.constants as _const
from openpnm.utils import Docorator
from openpnm.models.physics._utils import _poisson_conductance


docstr = Docorator()
__all__ = [
    "generic_diffusive",
    "ordinary_diffusion",
    "mixed_diffusion",
    "taylor_aris_diffusion",
    "multiphase_diffusion",
]


@docstr.get_sections(base='models.physics.diffusive_conductance',
                     sections=['Returns'])
@docstr.dedent
def generic_diffusive(target,
                      pore_diffusivity="pore.diffusivity",
                      throat_diffusivity="throat.diffusivity",
                      size_factors="throat.diffusive_size_factors"):
    r"""
    Calculates the diffusive conductance of conduits in network.

    A conduit is defined as (1/2 pore - full throat - 1/2 pore).

    Parameters
    ----------
    %(models.target.parameters)s
    %(models.phys.pore_diffusivity.parameters)s
    %(models.phys.throat_diffusivity.parameters)s
    size_factors: str
        Dictionary key of the conduit diffusive size factors' values

    Returns
    -------
    g_diffusive : ndarray
        A numpy ndarray containing diffusive conductance values for conduits

    """
    return _poisson_conductance(target=target,
                                pore_conductivity=pore_diffusivity,
                                throat_conductivity=throat_diffusivity,
                                size_factors=size_factors)


@docstr.dedent
def ordinary_diffusion(target,
                       pore_diffusivity="pore.diffusivity",
                       throat_diffusivity="throat.diffusivity",
                       size_factors="throat.diffusive_size_factors"):
    r"""
    Calculates the diffusive conductance of conduits in network.

    A conduit is defined as (1/2 pore - full throat - 1/2 pore).

    Parameters
    ----------
    %(models.target.parameters)s
    %(models.phys.pore_diffusivity.parameters)s
    %(models.phys.throat_diffusivity.parameters)s
    size_factors: str
        Dictionary key of the conduit diffusive size factors' values.

    Returns
    -------
    %(models.physics.diffusive_conductance.returns)s

    """
    return _poisson_conductance(target=target,
                                pore_conductivity=pore_diffusivity,
                                throat_conductivity=throat_diffusivity,
                                size_factors=size_factors)

@docstr.dedent
def mixed_diffusion(target,
                    pore_diameter="pore.diameter",
                    throat_diameter="throat.diameter",
                    pore_diffusivity="pore.diffusivity",
                    throat_diffusivity="throat.diffusivity",
                    pore_temperature="pore.temperature",
                    throat_temperature="throat.temperature",
                    molecular_weight="pore.molecular_weight",
                    size_factors="throat.diffusive_size_factors"):
    r"""
    Calculates the diffusive conductance of conduits in network with
    Knudsen correction.

    A conduit is defined as (1/2 pore - full throat - 1/2 pore). See
    Notes section for the limitations of this method.

    Parameters
    ----------
    %(models.target.parameters)s
    %(models.geometry.Pdia)s
    %(models.geometery.Tdia)s
    %(models.phys.pore_diffusivity.parameters)s
    %(models.phys.throat_diffusivity.parameters)s
    %(models.phys.pore_temperature.parameters)s
    %(models.phys.throat_temperature.parameters)s
    molecular_weigth : str
        Dictionary key of the pore molecular weight values.
    size_factors: str
        Dictionary key of the conduit diffusive size factors' values.

    Returns
    -------
    %(models.physics.diffusive_conductance.returns)s

    Notes
    -----
    This function is ONLY suitable for dilute mixtures and NOT those with
    concentrated species.

    This function requires that all the necessary phase properties already
    be calculated.

    """
    # Fetch GenericPhysicss
    network = target.network
    phase = target.project.find_phase(target)

    # Fetch model parameters
    Dp = phase[pore_diffusivity]
    Dt = phase[throat_diffusivity]
    dp = network[pore_diameter]
    dt = network[throat_diameter]
    MWp = phase[molecular_weight]
    MWt = phase.interpolate_data(propname=molecular_weight)
    Tp = phase[pore_temperature]
    Tt = phase[throat_temperature]

    # Calculate Knudsen contribution
    DKp = dp/3 * (8*_const.R*Tp / (_const.pi*MWp)) ** 0.5
    DKt = dt/3 * (8*_const.R*Tt / (_const.pi*MWt)) ** 0.5

    # Calculate mixed diffusivity
    Dp_eff = (1/DKp + 1/Dp) ** -1
    Dt_eff = (1/DKt + 1/Dt) ** -1

    return _poisson_conductance(target=target,
                                pore_conductivity=Dp_eff,
                                throat_conductivity=Dt_eff,
                                size_factors=size_factors)


@docstr.dedent
def taylor_aris_diffusion(target,
                          pore_area="pore.area",
                          throat_area="throat.area",
                          pore_diffusivity="pore.diffusivity",
                          throat_diffusivity="throat.diffusivity",
                          pore_pressure="pore.pressure",
                          throat_hydraulic_conductance="throat.hydraulic_conductance",
                          size_factors="throat.diffusive_size_factors"):
    r"""
    Calculates the diffusive conductance of conduits in network
    considering the Taylor-Aris effect (effect of flow on diffusion).

    A conduit is defined as (1/2 pore - full throat - 1/2 pore)

    Parameters
    ----------
    %(models.target.parameters)s
    %(models.geometry.Parea)s
    %(models.geometry.Tarea)s
    %(models.phys.pore_diffusivity.parameters)s
    %(models.phys.throat_diffusivity.parameters)s
    %(models.phys.pore_pressure.parameters)s
    throat_hydraulic_conductance : str
        Dictionary key of the throat hydraulic_conductance values.
    size_factors: str
        Dictionary key of the conduit size factors' values.

    Returns
    -------
    %(models.physics.diffusive_conductance.returns)s

    Notes
    -----
    This function requires that all the necessary phase properties are
    already calculated.

    """
    # Fetch GenericPhysicss
    network = target.network
    domain = target._domain
    throats = domain.throats(target.name)
    phase = target.project.find_phase(target)
    cn = network['throat.conns'][throats]
    F = network[size_factors]

    # Fetch model parameters
    A1, A2 = network[pore_area][cn].T
    At = network[throat_area][throats]
    D1, D2 = phase[pore_diffusivity][cn].T
    Dt = phase[throat_diffusivity][throats]
    P = phase[pore_pressure]
    gh = phase[throat_hydraulic_conductance]

    # Calculate Peclet number
    Qt = -gh * _np.diff(P[cn], axis=1).squeeze()
    u1, ut, u2 = [Qt/Ai for Ai in [A1, At, A2]]
    Pe1 = u1 * ((4 * A1 / _np.pi) ** 0.5) / D1
    Pe2 = u2 * ((4 * A2 / _np.pi) ** 0.5) / D2
    Pet = ut * ((4 * At / _np.pi) ** 0.5) / Dt

    # Calculate diffusive conductance
    if isinstance(F, dict):
        g1 = D1 * (1 + Pe1**2 / 192) * F[f"{size_factors}.pore1"][throats]
        gt = Dt * (1 + Pet**2 / 192) * F[f"{size_factors}.throat"][throats]
        g2 = D2 * (1 + Pe2**2 / 192) * F[f"{size_factors}.pore2"][throats]
        return 1 / (1/g1 + 1/gt + 1/g2)
    return Dt * (1 + Pet**2 / 192) * F[throats]


@docstr.dedent
def multiphase_diffusion(target,
                         pore_diffusivity="pore.diffusivity",
                         throat_diffusivity="throat.diffusivity",
                         size_factors="throat.diffusive_size_factors",
                         partition_coef_global="throat.partition_coef.all"):
    r"""
    Calculates the diffusive conductance of conduits in network.

    Parameters
    ----------
    %(models.target.parameters)s
    %(models.phys.pore_diffusivity.parameters)s
    %(models.phys.throat_diffusivity.parameters)s
    size_factors: str
        Dictionary key of the conduit size factors' values.

    Returns
    -------
    %(models.physics.diffusive_conductance.returns)s

    Notes
    -----
    This method assumes that phase["partition_coef"] contains information on
    binary phase partitioning. See ``MultiPhase`` class documentation for more
    information.

    """
    # Fetch GenericPhysicss
    network = target.network
    throats = target.throats(to_global=True)
    phase = target.project.find_phase(target)
    cn = network.conns[throats]
    F = network[size_factors]

    # Fetch model parameters
    D1, D2 = phase[pore_diffusivity][cn].T
    Dt = phase[throat_diffusivity][throats]
    g1 = D1 * F[f"{size_factors}.pore1"][throats]
    gt = Dt * F[f"{size_factors}.throat"][throats]
    g2 = D2 * F[f"{size_factors}.pore2"][throats]

    # Apply Henry's partitioning coefficient
    # Note: m12 = G21*c1 - G12*c2 NOT G12*c1 - G21*c2
    K12 = phase[partition_coef_global][throats]
    G21 = (1/g1 + 0.5/gt + K12 * (1/g2 + 0.5/gt)) ** -1
    G12 = K12 * G21

    return _np.vstack((G12, G21)).T
