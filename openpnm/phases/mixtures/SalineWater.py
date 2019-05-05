from openpnm.phases.mixtures import GenericMixture, species
from openpnm import models as mods
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class SalineWater(GenericMixture):
    r"""
    """
    def __init__(self, network, components=None, **kwargs):
        if components is not None:
            logger.warn('Ignoring received components')
        super().__init__(network=network, components=[], **kwargs)

        Cl = species.ions.Cl(network=network, name='Cl')
        Na = species.ions.Na(network=network, name='Na')
        W = species.liquids.H2O(network=network, name='H2O')
        self.settings['components'] = [Cl.name, Na.name, W.name]
        self.set_concentration(component=W, values=998/0.018)
        self.set_concentration(component=Na, values=0.0)
        self.set_concentration(component=Cl, values=0.0)
        self.add_model(propname='pore.salt_concentration',
                       model=mods.misc.summation,
                       props=['pore.concentration.Na',
                              'pore.concentration.Cl'])
        self.add_model(propname='pore.salinity',
                       model=mods.phases.mixtures.salinity,
                       concentration='pore.concentration.Na')
        self.add_model(propname='pore.mass_density',
                       model=mods.phases.density.water,
                       salinity='pore.salinity')
        self.add_model(propname='pore.viscosity',
                       model=mods.phases.viscosity.water,
                       salinity='pore.salinity')
        self.update_mole_fractions()
