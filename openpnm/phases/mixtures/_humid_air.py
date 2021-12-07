from openpnm.phases.mixtures import IdealGas, species
import openpnm.models as mods
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class HumidAir(IdealGas):
    r"""
    """

    def __init__(self, network, **kwargs):
        super().__init__(network=network, components=[], **kwargs)

        N2 = species.gases.N2(network=network, name='N2_'+self.name)
        O2 = species.gases.O2(network=network, name='O2_'+self.name)
        H2O = species.liquids.H2O(network=network, name='H2O_'+self.name)
        self.settings['components'] = [O2.name, N2.name, H2O.name]
        self.set_mole_fraction(component=N2, values=0.791)
        self.set_mole_fraction(component=O2, values=0.209)
        self.set_mole_fraction(component=H2O, values=0.000)
        self.add_model(propname='pore.vapor_pressure',
                       model=mods.phases.vapor_pressure.water)
#        self.add_model(propname='pore.mole_fraction.'+H2O.name,
#                       model=mods.misc.fraction,
#                       numerator='pore.vapor_pressure',
#                       denominator='pore.pressure')
        self.add_model(propname='pore.diffusivity',
                       model=mods.phases.mixtures.wilke_fuller_diffusivity)
