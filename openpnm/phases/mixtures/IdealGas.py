# from collections import ChainMap  # Might use eventually
from openpnm.phases.mixtures import GenericMixture
import openpnm.models as mods
from openpnm.utils import logging, Docorator
docstr = Docorator()
logger = logging.getLogger(__name__)


@docstr.dedent
class IdealGas(GenericMixture):
    r"""
    Creates Mixture object that represents a ideal gas system
    consisting of a given list of OpenPNM Phase objects as components.

    Parameters
    ----------
    %(GenericMixture.parameters)s

    """

    def __init__(self, settings={}, **kwargs):
        super().__init__(settings={'prefix': 'mix'}, **kwargs)
        self.settings.update(settings)

        self.add_model(propname='pore.molar_density',
                       model=mods.phases.molar_density.ideal_gas,
                       regen_mode='deferred')
