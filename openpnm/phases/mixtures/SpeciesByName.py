from openpnm.phases.mixtures import GenericSpecies
from openpnm.utils import PrintableDict
import chemicals as chem
from chemicals import numba_vectorized
from chemicals import Vm_to_rho
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
        if CAS in extra_LJ.keys():
            s, e_k = extra_LJ[CAS]
        else:
            e_k = chem.lennard_jones.Stockmayer(CAS)
            s = chem.lennard_jones.molecular_diameter(CAS)
        if e_k is not None:
            self.parameters['lennard_jones_epsilon'] = e_k*k
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


extra_LJ = {}
extra_LJ['7727-37-9'] = (3.788, 71.4)


def vapor_pressure(target, temperature='pore.temperature'):
    T = target[temperature]
    CAS = target.settings['CAS']
    Tc = target.parameters['critical_temperature']
    try:
        coeffs = chem.vapor_pressure.Psat_data_AntoineExtended.loc[CAS]
        _, A, B, C, Tc, to, n, E, F, Tmin, Tmax = coeffs
        PV = numba_vectorized.TRC_Antoine_extended(T, A, B, C, n, E, F)
    except KeyError:
        coeffs = chem.vapor_pressure.Psat_data_AntoinePoling.loc[CAS]
        _, A, B, C, Tmin, Tmax = coeffs
        PV = 10**(A - B/(T + C))
        # PV = numba_vectorized.vapor_pressure.Antoine(T=T, A=A, B=B, C=C)
    return PV


def liquid_density(target, temperature='pore.temperature'):
    T = target[temperature]
    MW = target.parameters['molecular_weight']
    Tc = target.parameters['critical_temperature']
    Vc = target.parameters['critical_volume']
    omega = target.parameters['acentric_factor']
    Vm = numba_vectorized.COSTALD(T, Tc, Vc, omega)
    rhoL = Vm_to_rho(Vm, MW)
    return rhoL


def liquid_viscosity(target, temperature='pore.temperature'):
    T = target[temperature]
    MW = target.parameters['molecular_weight']
    Tc = target.parameters['critical_temperature']
    Pc = target.parameters['critical_pressure']
    omega = target.parameters['acentric_factor']
    muL = numba_vectorized.Letsou_Stiel(T, MW, Tc, Pc, omega)
    return muL


def gas_viscosity(target, temperature='pore.temperature'):
    T = target[temperature]
    MW = target.parameters['molecular_weight']
    Tc = target.parameters['critical_temperature']
    Pc = target.parameters['critical_pressure']
    muG = numba_vectorized.viscosity_gas_Gharagheizi(T, Tc, Pc, MW)
    return muG


def liquid_thermal_conductivity(target, temperature='pore.temperature'):
    T = target[temperature]
    MW = target.parameters['molecular_weight']
    Tb = target.parameters['boiling_temperature']
    Pc = target.parameters['critical_pressure']
    omega = target.parameters['acentric_factor']
    kL = numba_vectorized.Gharagheizi_liquid(T, MW, Tb, Pc, omega)
    return kL


def gas_thermal_conductivity(target, temperature='pore.temperature'):
    T = target[temperature]
    MW = target.parameters['molecular_weight']
    Tb = target.parameters['boiling_temperature']
    Pc = target.parameters['critical_pressure']
    omega = target.parameters['acentric_factor']
    kG = numba_vectorized.Gharagheizi_gas(T, MW, Tb, Pc, omega)
    return kG


def liquid_heat_capacity(target,
                         Cpgm='pore.gas_heat_capacity',
                         T='pore.temperature'):
    T = target['pore.temperature']
    Cpgm = target['pore.gas_heat_capacity']
    Tc = target.parameters['critical_temperature']
    omega = target.parameters['acentric_factor']
    Cplm = numba_vectorized.Rowlinson_Poling(T, Tc, omega, Cpgm)
    return Cplm


def gas_heat_capacity(target, temperature='pore.temperature'):
    T = target[temperature]
    props = chem.heat_capacity.TRC_gas_data.loc[target.settings['CAS']]
    _, Tmin, Tmax, a0, a1, a2, a3, a4, a5, a6, a7, I, J, Hfg = props
    Cp = numba_vectorized.TRCCp(T, a0, a1, a2, a3, a4, a5, a6, a7)
    return Cp
