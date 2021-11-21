from openpnm.physics import GenericPhysics
from openpnm.models import physics as mods
from openpnm.utils import Docorator

docstr = Docorator()
__all__ = ["Basic"]


@docstr.dedent
class Basic(GenericPhysics):
    r"""
    Minimal subclass of GenericPhysics for performing diffusion and/or
    flow simulations.

    Parameters
    ----------
    %(GenericPhysics.parameters)s

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.add_model(propname='throat.hydraulic_conductance',
                       model=mods.hydraulic_conductance.generic_hydraulic)
        self.add_model(propname='throat.diffusive_conductance',
                       model=mods.diffusive_conductance.generic_diffusive)
        self.add_model(propname='throat.entry_pressure',
                       model=mods.capillary_pressure.washburn)
