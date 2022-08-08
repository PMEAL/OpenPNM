from openpnm.utils import Docorator
import numpy as np


docstr = Docorator()


__all__ = [
    "liquid_mixture_tc",
    "gas_mixture_ce",
    "gas_mixture_fesg",
    # "gas_mixture_wilke_fuller",
]


def liquid_mixture_tc(
    phase,
    T='pore.temperature',
    mu='pore.viscosity',
    Vms_at_Tb='param.molar_volume_Tb.*',
    sigmas_at_Tb='param.surface_tension_Tb.*',
):
    r"""
    Uses Tyn-Calus model to estimate diffusion coefficient in a dilute liquid
    solution of A in B from first principles at conditions of interest

    Parameters
    ----------
    %(models.target.parameters)s
    %(models.phase.T)s
    Vms_at_Tb : str or list of scalars
        Molar volumes of each component at its boiling temperature (m3/mol).
        Can either be a string to a parameter on each component of a list
        of scalar values which are used directly
    sigmas_at_Tb : float
        Surface tension of component A at boiling temperature (N/m)

    Returns
    -------
    %(models.phase.diffusivity.returns)s
    """
    T = phase[T]
    mu = phase[mu]
    sigma_A, sigma_B = phase.get_comp_vals(sigmas_at_Tb).values()
    VA, VB = phase.get_comp_vals(Vms_at_Tb).values()
    A = 8.93e-8*(VB*1e6)**0.267/(VA*1e6)**0.433*T
    B = (sigma_B/sigma_A)**0.15/(mu*1e3)
    value = A*B*1e-4
    return value


def gas_mixture_ce(
    phase,
    T='pore.temperature',
    P='pore.pressure',
    Tcs='param.critical_temperature.*',
    Pcs='param.critical_pressure.*',
    omegas='param.acentric_factor.*',
    MWs='param.molecular_weight.*',
    epsilons='param.LJ_energy.*',
    sigmas='param.LJ_diameter.*',
):
    r"""
    Calculate gas phase diffusion coefficient using Chapman-Enskog equation.

    Parameters
    ----------
    %(models.target.parameters)s
    MA, MB : scalar
        The molecular mass of species A and B in units of kg/mol
    sigma_AB : scalar
        The collision diameter in units of Angstroms
    %(models.phase.T)s
    %(models.phase.P)s

    Returns
    -------
    %(models.phase.diffusivity.returns)s
    """
    # Fetch values from components
    T = phase[T]
    P = phase[P]
    try:
        MA, MB = phase.get_comp_vals(MWs).values()
    except ValueError:
        raise Exception('This function only works on binary mixtures')
    MWAB = 2/(1/MA + 1/MB)
    omega = phase[omegas].values()
    Tc = phase[Tcs].values()
    Pc = phase[Pcs].values()
    k = 1.380649e-23  # Boltzmann constant
    try:
        eA, eB = phase.get_comp_vals(epsilons).values()
    except Exception:
        # Use correlation of Tee, Gotoh, & Stewart
        eA, eB = (0.7915 + 0.1693*omega)*Tc*k
    eAB = (eA*eB)**0.5
    # Compute collision integral using Neufeld's correlation (RPP Eq.(11-3.6))
    A, B, C, D, E, F, G, H = (1.06036, 0.15610, 0.19300, 0.47635, 1.03587,
                              1.52996, 1.76474, 3.89411)
    Tstar = k*T/eAB
    Omega = Omega = A/(Tstar**B) + C/np.exp(D*Tstar) + E/np.exp(F*Tstar) \
        + G/np.exp(H*Tstar)
    # Now apply RPP Eq.(11-3.2)
    try:
        sA, sB = phase.get_comp_vals(sigmas).values()
    except:
        # Use correlation of Tee, Gotoh, & Stewart
        sA, sB = (2.3551 - 0.08847*omega)*(Tc/Pc)**(1/3)
    sAB = (sA + sB)/2
    DAB = 0.00266 * (T**1.5) / (P/101325 * MWAB**0.5 * sAB**2 * Omega) * 1e-4
    return DAB


def gas_mixture_fesg(
    phase,
    T='pore.temperature',
    P='pore.pressure',
    MWs='param.molecular_weight.*',
    Vdms='param.molar_diffusion_volume.*',
):
    r"""
    Estimates the diffusion coefficient of both species in a binary gas
    mixture using the Fuller et al correlation [1]

    Parameters
    ----------
    %(models.phase.parameters)s
    molecular_weight : string
        Dictionary key containing the molecular weight of each species. The
        default is 'pore.molecular_weight'
    molar_diffusion_volume : string
        Dictionary key containing the molar diffusion volume of each species.
        This is used by the Fuller correlation. The default is
        'pore.molar_diffusion_volume'
    %(models.phase.T)s
    %(models.phase.P)s

    Returns
    -------
    Dij : dict containing ND-arrys
        The dict contains one array for each component, containing the
        diffusion coefficient of that component at each location.

    References
    ----------
    [1] Fuller, E. N., and J. C. Giddings: J. Gas Chromatogr., 3: 222 (1965).
    [2] Fuller, E. N., P. D. Schettler, and J. C. Giddings: Ind. Eng. Chem.,
        58(5): 18 (1966).
    [3] Fuller, E. N., K. Ensley, and J. C. Giddings: J. Phys. Chem.,
        73: 3679 (1969).
    """
    T = phase[T]
    P = phase[P]
    # The following is to accomodate proper mixtures AND normal Phase objects
    try:
        MA, MB = phase[MWs].values()
    except AttributeError:
        MA, MB = phase[MWs]
    try:
        vA, vB = phase[Vdms].values()
    except AttributeError:
        vA, vB = phase[Vdms]
    MAB = 2/(1/MA + 1/MB)
    P = P/101325
    DAB = 0.00143*T**1.75/(P*(MAB**0.5)*(vA**(1./3) + vB**(1./3))**2)*1e-4
    return DAB


def gas_mixture_fw(
    phase,
    T='pore.temperature',
    P='pore.pressure',
    MWs='param.molecular_weight.*',
    Vdms='pore.molar_diffusion_volume.*',
):
    r"""
    Estimates the diffusion coeffient of each species in a gas mixture

    Uses the Fuller equation to estimate binary diffusivity between pairs,
    then uses the correction of Fairbanks and Wilke [1] to account for the
    composition of the gas mixture.

    Parameters
    ----------
    %(models.target.parameters)s
    molecular_weight : str
        Dictionary key containing the molecular weight of each species.  The
        default is 'pore.molecular_weight'
    molar_diffusion_volume : str
        Dictionary key containing the molar diffusion volume of each species.
        This is used by the Fuller correlation.  The default is
        'pore.molar_diffusion_volume'
    %(models.phase.T)s
    %(models.phase.P)s

    Returns
    -------
    Dij : dict containing ND-arrys
        The dict contains one array for each component, containing the
        diffusion coefficient of that component at each location.

    Reference
    ---------
    [1] Fairbanks DF and CR Wilke, Diffusion Coefficients in Multicomponent
    Gas Mixtures. Industrial & Engineering Chemistry, 42(3), p471–475 (1950).
    `DOI: 10.1021/ie50483a022 <http://doi.org/10.1021/ie50483a022>`_
    """
    raise Exception('This function is not ready yet')
    comps = list(phase.components.values())
    values = {}
    for i in range(len(comps)):
        A = comps[i]
        denom = 0.0
        for j in range(len(comps)):
            if i != j:
                B = comps[j]
                D = gas_mixture_fesg(phase=phase,
                                     molecular_weight=MW_pair,
                                     molar_diffusion_volume=Vdm_pair,
                                     temperature=T,
                                     pressure=P)
                yB = phase['pore.mole_fraction.' + B.name]
                denom += yB/D
        yA = phase['pore.mole_fraction.' + A.name]
        values[A.name] = (1 - yA)/denom
    return values
