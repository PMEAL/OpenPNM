import numpy as _np
import scipy.constants as _const
from openpnm.models import _doctxt
from openpnm.models.physics._utils import _poisson_conductance


__all__ = [
    "generic_diffusive",
    "ordinary_diffusion",
    "mixed_diffusion",
    "taylor_aris_diffusion",
]


@_doctxt
def generic_diffusive(phase,
                      pore_diffusivity="pore.diffusivity",
                      throat_diffusivity="throat.diffusivity",
                      size_factors="throat.diffusive_size_factors"):
    r"""
    Calculates the diffusive conductance of conduits in network.

    Parameters
    ----------
    %(phase)s
    pore_diffusivity : str
        %(dict_blurb)s pore diffusivity
    throat_diffusivity : str
        %(dict_blurb)s throat diffusivity
    size_factors : str
        %(dict_blurb)s conduit diffusive size factors

    Returns
    -------
    %(return_arr)s diffusive conductance

    """
    return _poisson_conductance(phase=phase,
                                pore_conductivity=pore_diffusivity,
                                throat_conductivity=throat_diffusivity,
                                size_factors=size_factors)


@_doctxt
def ordinary_diffusion(phase,
                       pore_diffusivity="pore.diffusivity",
                       throat_diffusivity="throat.diffusivity",
                       size_factors="throat.diffusive_size_factors"):
    r"""
    Calculates the diffusive conductance of conduits in network.

    Parameters
    ----------
    %(phase)s
    pore_diffusivity : str
        %(dict_blurb)s pore diffusivity
    throat_diffusivity : str
        %(dict_blurb)s throat diffusivity
    size_factors: str
        %(dict_blurb)s conduit diffusive size factors

    Returns
    -------
    %(return_arr)s diffusive conductance

    """
    return _poisson_conductance(phase=phase,
                                pore_conductivity=pore_diffusivity,
                                throat_conductivity=throat_diffusivity,
                                size_factors=size_factors)


@_doctxt
def mixed_diffusion(phase,
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

    Parameters
    ----------
    %(phase)s
    pore_diameter : str
        %(dict_blurb)s pore diameter
    throat_diameter : str
        %(dict_blurb)s throat diameter
    pore_diffusivity : str
        %(dict_blurb)s pore diffusivity
    throat_diffusivity : str
        %(dict_blurb)s pore diffusivity
    pore_temperature : str
        %(dict_blurb)s pore temperature
    molecular_weigth : str
        %(dict_blurb)s pore molecular weight
    size_factors : str
        %(dict_blurb)s size factor

    Returns
    -------
    %(return_arr)s diffusive conductance

    """
    # Fetch model parameters
    Dp = phase[pore_diffusivity]
    Dt = phase[throat_diffusivity]
    dp = phase.network[pore_diameter]
    dt = phase.network[throat_diameter]
    MWp = phase[molecular_weight]
    MWt = phase.interpolate_data(propname='throat.'+molecular_weight.split('.', 1)[-1])
    Tp = phase[pore_temperature]
    Tt = phase[throat_temperature]

    # Calculate Knudsen contribution
    DKp = dp/3 * (8*_const.R*Tp / (_const.pi*MWp)) ** 0.5
    DKt = dt/3 * (8*_const.R*Tt / (_const.pi*MWt)) ** 0.5

    # Calculate mixed diffusivity
    Dp_eff = (1/DKp + 1/Dp) ** -1
    Dt_eff = (1/DKt + 1/Dt) ** -1

    return _poisson_conductance(phase=phase,
                                pore_conductivity=Dp_eff,
                                throat_conductivity=Dt_eff,
                                size_factors=size_factors)


@_doctxt
def taylor_aris_diffusion(phase,
                          pore_area="pore.area",
                          throat_area="throat.cross_sectional_area",
                          pore_diffusivity="pore.diffusivity",
                          throat_diffusivity="throat.diffusivity",
                          pore_pressure="pore.pressure",
                          throat_hydraulic_conductance="throat.hydraulic_conductance",
                          size_factors="throat.diffusive_size_factors"):
    r"""
    Calculates the diffusive conductance of conduits in network
    considering the Taylor-Aris effect (effect of flow on diffusion).

    Parameters
    ----------
    %(phase)s
    pore_area : str
        %(dict_blurb)s pore cross-sectional area
    throat_area : str
        %(dict_blurb)s throat cross-sectional area
    pore_diffusivity : str
        %(dict_blurb)s pore diffusivity
    throat_diffusivity : str
        %(dict_blurb)s pore diffusivity
    pore_pressure : str
        %(dict_blurb)s pore pressure
    throat_hydraulic_conductance : str
        %(dict_blurb)s throat hydraulic conductance
    size_factors : str
        %(dict_blurb)s conduit size factors

    Returns
    -------
    %(return_arr)s diffusive conductance

    """
    network = phase.network
    cn = network['throat.conns']
    F = network[size_factors]

    # Fetch model parameters
    A1, A2 = network[pore_area][cn].T
    At = network[throat_area]
    D1, D2 = phase[pore_diffusivity][cn].T
    Dt = phase[throat_diffusivity]
    P = phase[pore_pressure]
    gh = phase[throat_hydraulic_conductance]

    # Calculate Peclet number
    Qt = -gh * _np.diff(P[cn], axis=1).squeeze()
    u1, ut, u2 = [Qt/Ai for Ai in [A1, At, A2]]
    Pe1 = u1 * ((4 * A1 / _np.pi) ** 0.5) / D1
    Pe2 = u2 * ((4 * A2 / _np.pi) ** 0.5) / D2
    Pet = ut * ((4 * At / _np.pi) ** 0.5) / Dt

    if F.ndim > 1:
        g1 = D1 * (1 + Pe1**2 / 192) * F[:, 0]
        gt = Dt * (1 + Pet**2 / 192) * F[:, 1]
        g2 = D2 * (1 + Pe2**2 / 192) * F[:, 2]
        gtot = 1 / (1/g1 + 1/gt + 1/g2)
    else:
        gtot = Dt * (1 + Pet**2 / 192) * F
    return gtot
