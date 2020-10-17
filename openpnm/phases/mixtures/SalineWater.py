from openpnm.phases.mixtures import GenericMixture, species
from openpnm import models as mods
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class SalineWater(GenericMixture):

    def __init__(self, network, components=None, **kwargs):
        if components is not None:
            logger.warn('Ignoring received components')
        super().__init__(network=network, components=[], **kwargs)

        Cl = species.ions.Cl(network=network, name='Cl_'+self.name)
        Na = species.ions.Na(network=network, name='Na_'+self.name)
        W = species.liquids.H2O(network=network, name='H2O_'+self.name)
        self.set_component([Cl, Na, W])
        self.set_concentration(component=W, values=998/0.018)
        self.set_concentration(component=Na, values=0.0)
        self.set_concentration(component=Cl, values=0.0)
        self.update_mole_fractions()
        self.add_model(propname='pore.salt_concentration',
                       model=mods.misc.summation,
                       props=['pore.concentration.Na_'+self.name,
                              'pore.concentration.Cl_'+self.name])
        self.add_model(propname='pore.salinity',
                       model=mods.phases.mixtures.salinity,
                       concentration='pore.concentration.Na_'+self.name)
        self.add_model(propname='pore.mass_density',
                       model=mods.phases.density.water,
                       salinity='pore.salinity')
        self.add_model(propname='pore.viscosity',
                       model=mods.phases.viscosity.water,
                       salinity='pore.salinity')
        self.add_model(propname='pore.permittivity',
                       model=mods.phases.permittivity.water)
