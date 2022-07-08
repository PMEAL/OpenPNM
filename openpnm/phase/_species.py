from openpnm.phase import Phase
import chemicals as chem
from chemicals.utils import k
import openpnm.models as mods
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
        temp['CAS'] = CAS
        a = chem.identifiers.search_chemical(CAS)
        temp['common_name'] = a.common_name
        temp['molecular_weight'] = a.MW/1000  # Convert to kg/mol
        temp['critical_temperature'] = chem.critical.Tc(CAS)
        temp['critical_pressure'] = chem.critical.Pc(CAS)
        temp['critical_volume'] = chem.critical.Vc(CAS)
        temp['critical_compressibilty_factor'] = chem.critical.Zc(CAS)
        temp['boiling_temperature'] = chem.Tb(CAS)
        temp['melting_temperature'] = chem.Tm(CAS)
        temp['acentric_factor'] = chem.acentric.omega(CAS)
        temp['dipole_moment'] = chem.dipole.dipole_moment(CAS)
        temp['lennard_jones_epsilon'] = k*chem.lennard_jones.Stockmayer(CAS)
        temp['lennard_jones_sigma'] = chem.lennard_jones.molecular_diameter(CAS)
        super().__init__(**kwargs)
        self.params.update(temp)


class GasByName(SpeciesByName):
    r"""
    Creates a phase object based on given chemical name, including some
    additional properties for a gas
    """

    def __init__(self, species, **kwargs):
        super().__init__(species=species, **kwargs)
        self.add_model_collection(mods.collections.phase.generic_gas())
        self.regenerate_models()


class LiquidByName(SpeciesByName):
    r"""
        Creates a phase object based on given chemical name, including some
        additional properties for a liquid
    """

    def __init__(self, species, **kwargs):
        super().__init__(species=species, **kwargs)
        self.add_model_collection(mods.collections.phase.generic_liquid())
        self.regenerate_models()
