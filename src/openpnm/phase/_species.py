import logging

from openpnm.models.collections.phase import standard_gas, standard_liquid
from openpnm.phase import Phase, _fetch_chemical_props

logger = logging.getLogger(__name__)


__all__ = [
    'Species',
    'StandardGas',
    'StandardLiquid',
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
        ``species.params`` like ``species.params['molecular_weight']``. They
        can also be accessed using ``species['param.molecular_weight']`` using
        some behind the scenes python magic.
    name : str, optional
        The name of the phase. This is useful to keep track of the objects
        throughout the simulation.  The name must be unique to the project.
        If no name is given, one is generated.

    """

    def __init__(self, species, **kwargs):
        # Create temp first to ensure all look-ups pass before initializing obj
        from thermo import Chemical
        temp = _fetch_chemical_props(Chemical(species))
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


class StandardLiquid(Species):
    """A ``Species`` object with built-in standard liquid-phase models."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.add_model_collection(standard_liquid)
        self.regenerate_models()


class StandardGas(Species):
    """A ``Species`` object with built-in standard gas-phase models."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.add_model_collection(standard_gas)
        self.regenerate_models()
