from openpnm.phases.mixtures import GenericSpecies
from openpnm.utils import PrintableDict
import chemicals as chem
from chemicals.utils import R, k
from openpnm.utils import logging
import numpy as np
logger = logging.getLogger(__name__)


class SpeciesByName(GenericSpecies):
    r"""
    Creates Phase object that represents a single species in a multicomponent
    mixture system.

    Parameters
    ----------
    network : OpenPNM Network object
        The network to which this phase object will be attached.
    species : str
        The name of the species to generate.  This is used to lookup tabulated
        constants in the ``chemicals`` package, which attempts to find a match. For
        instance, 'water', 'Water', and 'H2O' all work.
    project : OpenPNM Project object, optional
        The Project with which this phase should be associted.  If a
        ``network`` is given then this is ignored and the Network's project
        is used.  If a ``network`` is not given then this is mandatory.
    name : string, optional
        The name of the phase.  This is useful to keep track of the objects
        throughout the simulation.  The name must be unique to the project.
        If no name is given, one is generated.

    """

    def __init__(self, species, **kwargs):
        super().__init__(**kwargs)
        CAS = chem.CAS_from_any(species)
        self.settings['CAS'] = CAS
        self.parameters = PrintableDict()
        self.parameters['molecular_weight'] = chem.MW(CAS)
        self.parameters['critical_temperature'] = chem.critical.Tc(CAS)
        self.parameters['critical_pressure'] = chem.critical.Pc(CAS)
        self.parameters['critical_volume'] = chem.critical.Vc(CAS)
        self.parameters['critical_compressibilty_factor'] = chem.critical.Zc(CAS)
        self.parameters['boiling_temperature'] = chem.Tb(CAS)
        self.parameters['melting_temperature'] = chem.Tm(CAS)
        self.parameters['acentric_factor'] = chem.acentric.omega(CAS)
        self.parameters['dipole_moment'] = chem.dipole.dipole_moment(CAS)
        e_k = chem.lennard_jones.Stockmayer(CAS)
        if e_k is not None:
            self.parameters['lennard_jones_epsilon'] = e_k*k
        s = chem.lennard_jones.molecular_diameter(CAS)
        if s is not None:
            self.parameters['lennard_jones_sigma'] = s


class GasByName(SpeciesByName):
    def __init__(self, species, **kwargs):
        super().__init__(species=species, **kwargs)
        self.add_model(propname='pore.heat_capacity',
                       model=gas_heat_capacity)
        self.add_model(propname='pore.thermal_conductivity',
                       model=gas_thermal_conductivity)
        self.add_model(propname='pore.viscosity',
                       model=gas_viscosity)
        self.add_model(propname='pore.vapor_pressure',
                       model=vapor_pressure)


class LiquidByName(SpeciesByName):
    def __init__(self, species, **kwargs):
        super().__init__(species=species, **kwargs)
        self.add_model(propname='pore.gas_heat_capacity',
                       model=gas_heat_capacity)
        self.add_model(propname='pore.heat_capacity',
                       model=liquid_heat_capacity)
        self.add_model(propname='pore.thermal_conductivity',
                       model=liquid_thermal_conductivity)
        self.add_model(propname='pore.viscosity',
                       model=liquid_viscosity)
        self.add_model(propname='pore.density',
                       model=liquid_density)
        self.add_model(propname='pore.vapor_pressure',
                       model=vapor_pressure)


def vapor_pressure(target, temperature='pore.temperature'):
    T = target[temperature]
    CAS = target.settings['CAS']
    Tc = target.parameters['critical_temperature']
    try:
        coeffs = chem.vapor_pressure.Psat_data_AntoineExtended.loc[CAS]
        _, A, B, C, Tc, to, n, E, F, Tmin, Tmax = coeffs
        x = (T - to - 273.15)/Tc
        x = np.clip(x, 0, np.inf)
        x4 = x*x*x*x
        PV = 10.**(A - B/(T+C) + 0.43429*x**n + x4*x4*(E + F*x4))
    except KeyError:
        coeffs = chem.vapor_pressure.Psat_data_AntoinePoling.loc[CAS]
        _, A, B, C, Tmin, Tmax = coeffs
        PV = chem.vapor_pressure.Antoine(T=T, A=A, B=B, C=C)
    return PV


def liquid_density(target, temperature='pore.temperature'):
    from chemicals import Vm_to_rho
    T = target[temperature]
    MW = target.parameters['molecular_weight']
    Tc = target.parameters['critical_temperature']
    Vc = target.parameters['critical_volume']
    omega = target.parameters['acentric_factor']
    T[T > Tc] = Tc
    Tr = T/Tc
    tau = 1.0 - Tr
    tau_cbrt = (tau)**(1.0/3.)
    a = 0.296123
    b = 0.0480645
    c = 0.0427258
    d = 0.386914
    e = 0.190454
    f = 0.81446
    g = 1.43907
    h = 1.52816
    V_delta = (-a + Tr*(Tr*(-b*Tr - c) + d))/(Tr - 1.00001)
    V_0 = tau_cbrt*(tau_cbrt*(tau_cbrt*(e*tau_cbrt - f) + g) - h) + 1.0
    Vm = Vc*V_0*(1.0 - omega*V_delta)
    rhoL = Vm_to_rho(Vm=Vm, MW=MW)
    return rhoL


def liquid_viscosity(target, temperature='pore.temperature'):
    T = target[temperature]
    MW = target.parameters['molecular_weight']
    Tc = target.parameters['critical_temperature']
    Pc = target.parameters['critical_pressure']
    omega = target.parameters['acentric_factor']
    muL = chem.viscosity.Letsou_Stiel(T=T, MW=MW, Tc=Tc, Pc=Pc, omega=omega)
    return muL


def gas_viscosity(target, temperature='pore.temperature'):
    T = target[temperature]
    MW = target.parameters['molecular_weight']
    Tc = target.parameters['critical_temperature']
    Pc = target.parameters['critical_pressure']
    Tr = T/Tc
    mask = Tr < 0.2
    Tr[mask] = 0.2
    T[mask] = 0.2*Tc
    muG = 1E-5*Pc*Tr + (0.091 - 0.477/MW)*T + \
        MW*(1E-5*Pc - 8.0*MW*MW/(T*T))*(10.7639/Tc - 4.1929/T)
    muG = 1e-7*muG
    return muG


def liquid_thermal_conductivity(target, temperature='pore.temperature'):
    T = target[temperature]
    MW = target.parameters['molecular_weight']
    Tb = target.parameters['boiling_temperature']
    Pc = target.parameters['critical_pressure']
    omega = target.parameters['acentric_factor']
    kL = chem.thermal_conductivity.Gharagheizi_liquid(T=T, MW=MW, Tb=Tb,
                                                      Pc=Pc, omega=omega)
    return kL


def gas_thermal_conductivity(target, temperature='pore.temperature'):
    T = target[temperature]
    MW = target.parameters['molecular_weight']
    Tb = target.parameters['boiling_temperature']
    Pc = target.parameters['critical_pressure']
    omega = target.parameters['acentric_factor']
    kG = chem.thermal_conductivity.Gharagheizi_gas(T=T, MW=MW, Tb=Tb,
                                                   Pc=Pc, omega=omega)
    return kG


def liquid_heat_capacity(target,
                         Cpgm='pore.gas_heat_capacity',
                         T='pore.temperature'):
    T = target['pore.temperature']
    Cpgm = target['pore.gas_heat_capacity']
    Tc = target.parameters['critical_temperature']
    omega = target.parameters['acentric_factor']
    Tr = T/Tc
    Cplm = Cpgm + R*(1.586 + 0.49/(1.0-Tr) +
                     omega*(4.2775 + 6.3*(1.0-Tr)**(1/3.)/Tr + 0.4355/(1.0-Tr)))
    return Cplm


def gas_heat_capacity(target, temperature='pore.temperature'):
    from numpy import exp
    T = target[temperature]
    props = chem.heat_capacity.TRC_gas_data.loc[target.settings['CAS']]
    _, Tmin, Tmax, a0, a1, a2, a3, a4, a5, a6, a7, I, J, Hfg = props
    if np.any(T <= a7):
        y = 0.
    else:
        y = (T - a7)/(T + a6)
    T_inv = 1.0/T
    y2 = y*y
    T_m_a7 = T - a7
    Cp = R*(a0 + (a1*T_inv*T_inv)*exp(-a2*T_inv) +
            y2*(a3 + (a4 - a5/(T_m_a7*T_m_a7))*y2*y2*y2))
    return Cp
