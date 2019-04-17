import numpy as np
from openpnm.phases import GenericMixture, Water
from openpnm.phases.components.ions import Cl, Na
from openpnm.models.phases.mixtures import mole_weighted_average
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class SalineWater(GenericMixture):
    r"""
    """
    def __init__(self, network, components=None, **kwargs):
        if components is not None:
            logger.warn('Ignoring received components')
        super().__init__(network=network, components=[], **kwargs)

        C = Cl(network=network, name='Cl')
        N = Na(network=network, name='Na')
        W = Water(network=network, name='H2O')
        self.settings['components'] = [C.name, N.name, W.name]
        self.set_mole_fraction(component=W, values=1.0)
        self.set_mole_fraction(component=N, values=0.0)
        self.set_mole_fraction(component=C, values=0.0)
        self.add_model(propname='pore.molar_mass',
                       model=mole_weighted_average,
                       prop='pore.molecular_weight')
