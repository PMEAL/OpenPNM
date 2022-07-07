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
        logger.warn("No mixture phase found for this species")


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
        # Create temp first to ensure all look-ups pass before initializing obj
        temp = {}
        CAS = chem.CAS_from_any(species)
        temp['param.CAS'] = CAS
        a = chem.identifiers.search_chemical(CAS)
        temp['param.common_name'] = a.common_name
        temp['param.molecular_weight'] = a.MW/1000  # Convert to kg/mol
        temp['param.critical_temperature'] = chem.critical.Tc(CAS)
        temp['param.critical_pressure'] = chem.critical.Pc(CAS)
        temp['param.critical_volume'] = chem.critical.Vc(CAS)
        temp['param.critical_compressibilty_factor'] = chem.critical.Zc(CAS)
        temp['param.boiling_temperature'] = chem.Tb(CAS)
        temp['param.melting_temperature'] = chem.Tm(CAS)
        temp['param.acentric_factor'] = chem.acentric.omega(CAS)
        temp['param.dipole_moment'] = chem.dipole.dipole_moment(CAS)
        temp['param.lennard_jones_epsilon'] = k*chem.lennard_jones.Stockmayer(CAS)
        temp['param.lennard_jones_sigma'] = chem.lennard_jones.molecular_diameter(CAS)
        super().__init__(**kwargs)
        self._params.update(temp)


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
