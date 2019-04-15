import numpy as np
from openpnm.phases import GenericMixture, Water
from openpnm.phases.components.ions import Cl, Na
from openpnm.utils import logging
logger = logging.getLogger(__name__)


def mixture_mole_average(target, prop):
    # temp = set(target.components.keys()).difference({solvent})
    # solutes = []
    # for name in temp:
    #     solutes.append(target.project[name])
    # solvent = target.project[solvent]
    element = prop.split('.')[0]
    vals = np.zeros(target._count(element))
    for item in target.components.keys():
        frac = target[element + '.mole_fraction.' + item]
        temp = target.project[item][prop]
        vals += temp*frac
    return vals


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
                       model=mixture_mole_average,
                       prop='pore.molecular_weight')
