import numpy as np
from openpnm.phases.mixtures import GenericMixture, species
from openpnm.models.phases.mixtures import mole_weighted_average
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class Air(GenericMixture):
    r"""
    """
    def __init__(self, network, **kwargs):
        super().__init__(network=network, components=[], **kwargs)

        N2 = species.gases.N2(network=network, name='pure_N2')
        O2 = species.gases.O2(network=network, name='pure_O2')
        self.settings['components'] = [O2.name, N2.name]
        self.set_mole_fraction(component=N2, values=0.791)
        self.set_mole_fraction(component=O2, values=0.209)
        self.add_model(propname='pore.molar_mass',
                       model=mole_weighted_average,
                       prop='pore.molecular_weight')
