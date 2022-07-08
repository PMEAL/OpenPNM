from openpnm.phase import Phase
import chemicals as chem
from chemicals.utils import k
import logging


logger = logging.getLogger(__name__)


__all__ = [
    'Species',
]


class Species(Phase):
    r"""
    A special Phase object that represents a single species in a mixture

    This class provides a ``mixture`` attribute which allows one to lookup
    which mixture the species is associated with.

    Parameters
    ----------
    network : GenericNetwork
        The network to which this phase object will be attached.
    species : str, optional
        If provided, this is used to lookup tabulated constants from the
        ``chemicals`` package, which attempts to find a match. For instance,
        'water', 'Water', and 'H2O' all work. The contants are stored in
        ``species.params`` like ``species.params['molecular_weight']``.
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
        temp['charge'] = a.charge
        temp['formula'] = a.formula
        temp['molecular_weight'] = a.MW
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

    @property
    def mixture(self):
        for item in self.project:
            if hasattr(item, 'components'):
                for comp in item.components.values():
                    if self is comp:
                        return item
        logger.warn("No mixture phase found for this species")
