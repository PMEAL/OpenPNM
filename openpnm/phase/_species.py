from openpnm.phase import Phase
import chemicals as chem
from chemicals.utils import k
import openpnm.models.phase as mods
import logging


logger = logging.getLogger(__name__)


__all__ = [
    'Species',
    'SpeciesByName',
    'GasByName',
    'LiquidByName',
]


class Species(Phase):
    r"""
    Creates Phase object that represents a single species in a multicomponent
    mixture system.

    Parameters
    ----------
    network : Network
        The network to which this phase object will be attached.
    name : str, optional
        The name of the phase.  This is useful to keep track of the objects
        throughout the simulation.  The name must be unique to the project.
        If no name is given, one is generated.

    """
    @property
    def mixture(self):
        for item in self.project:
            if hasattr(item, 'components'):
                for comp in item.components.values():
                    if self is comp:
                        return item


class SpeciesByName(Species):
    r"""
    Creates Phase object that represents a single species in a multicomponent
    mixture system.

    Parameters
    ----------
    network : GenericNetwork
        The network to which this phase object will be attached.
    species : str
        The name of the species to generate.  This is used to lookup tabulated
        constants in the ``chemicals`` package, which attempts to find a match.
        For instance, 'water', 'Water', and 'H2O' all work.
    name : str, optional
        The name of the phase. This is useful to keep track of the objects
        throughout the simulation.  The name must be unique to the project.
        If no name is given, one is generated.

    """

    def __init__(self, species, **kwargs):
        super().__init__(**kwargs)
        CAS = chem.CAS_from_any(species)
        self['param.CAS'] = CAS
        a = chem.identifiers.search_chemical(CAS)
        self['param.common_name'] = a.common_name
        self['param.molecular_weight'] = a.MW/1000  # Convert to kg/mol
        self['param.critical_temperature'] = chem.critical.Tc(CAS)
        self['param.critical_pressure'] = chem.critical.Pc(CAS)
        self['param.critical_volume'] = chem.critical.Vc(CAS)
        self['param.critical_compressibilty_factor'] = chem.critical.Zc(CAS)
        self['param.boiling_temperature'] = chem.Tb(CAS)
        self['param.melting_temperature'] = chem.Tm(CAS)
        self['param.acentric_factor'] = chem.acentric.omega(CAS)
        self['param.dipole_moment'] = chem.dipole.dipole_moment(CAS)
        self['param.lennard_jones_epsilon'] = k*chem.lennard_jones.Stockmayer(CAS)
        self['param.lennard_jones_sigma'] = chem.lennard_jones.molecular_diameter(CAS)


class GasByName(SpeciesByName):
    r"""
    Creates a phase object based on given chemical name, including some
    additional properties for a gas
    """
    def __init__(self, species, **kwargs):
        super().__init__(species=species, **kwargs)
        self.add_model(propname='pore.heat_capacity',
                       model=mods.heat_capacity.gas_heat_capacity)
        self.add_model(propname='pore.thermal_conductivity',
                       model=mods.thermal_conductivity.gas_thermal_conductivity)
        self.add_model(propname='pore.viscosity',
                       model=mods.viscosity.gas_viscosity)


class LiquidByName(SpeciesByName):
    r"""
        Creates a phase object based on given chemical name, including some
        additional properties for a liquid
    """
    def __init__(self, species, **kwargs):
        super().__init__(species=species, **kwargs)
        self.add_model(propname='pore.heat_capacity',
                       model=mods.heat_capacity.liquid_heat_capacity)
        self.add_model(propname='pore.thermal_conductivity',
                       model=mods.thermal_conductivity.liquid_thermal_conductivity)
        self.add_model(propname='pore.viscosity',
                       model=mods.viscosity.liquid_viscosity)
        self.add_model(propname='pore.density',
                       model=mods.density.liquid_density)
        self.add_model(propname='pore.vapor_pressure',
                       model=mods.vapor_pressure.vapor_pressure)
