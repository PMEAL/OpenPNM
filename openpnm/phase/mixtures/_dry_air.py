import logging
from openpnm.phase.mixtures import IdealGas, species
import openpnm.models as mods
logger = logging.getLogger(__name__)


class DryAir(IdealGas):
    r"""
    """

    def __init__(self, network, **kwargs):
        super().__init__(network=network, components=[], **kwargs)

        N2 = species.gases.N2(network=network, name='N2_'+self.name)
        O2 = species.gases.O2(network=network, name='O2_'+self.name)
        self.set_component([O2, N2])
        self.set_mole_fraction(component=N2, values=0.791)
        self.set_mole_fraction(component=O2, values=0.209)
        self.add_model(propname='pore.diffusivity.N2',
                       model=mods.phase.mixtures.fuller_diffusivity)
        self.add_model(propname='pore.diffusivity.O2',
                       model=mods.phase.mixtures.fuller_diffusivity)
