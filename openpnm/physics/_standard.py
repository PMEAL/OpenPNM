from openpnm.physics import GenericPhysics
from openpnm.models import physics as mods
from openpnm.utils import logging, Docorator


logger = logging.getLogger(__name__)
docstr = Docorator()


@docstr.dedent
class Standard(GenericPhysics):
    r"""
    Generic class to generate Physics objects

    Parameters
    ----------
    %(GenericPhysics.parameters)s
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.add_model(propname='throat.hydraulic_conductance',
                       model=mods.hydraulic_conductance.hagen_poiseuille)
        self.add_model(propname='throat.diffusive_conductance',
                       model=mods.diffusive_conductance.mixed_diffusion)
        self.add_model(propname='throat.ad_dif_conductance',
                       model=mods.ad_dif_conductance.ad_dif)
        self.add_model(propname='throat.entry_pressure',
                       model=mods.capillary_pressure.washburn)
        self.add_model(propname='throat.thermal_conductance',
                       model=mods.thermal_conductance.generic_thermal)
        self.add_model(propname='throat.electrical_conductance',
                       model=mods.electrical_conductance.generic_electrical)
