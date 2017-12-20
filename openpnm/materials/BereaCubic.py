from openpnm.core import logging
from openpnm.network import Cubic
from openpnm.geometry import models as gm
logger = logging.getLogger(__name__)


class BereaCubic(Cubic):

    def __init__(self, shape, **kwargs):
        super().__init__(shape=shape, spacing=1e-4, connectivity=26, **kwargs)

        self.add_model(propname='pore.seed',
                       model=gm.pore_seed.spatially_correlated,
                       weights=[1, 1, 1])
        self.add_model(propname='pore.diameter',
                       model=gm.pore_size.weibull,
                       shape=.5, scale=5e-4, loc=1e-6)
        self.add_model(propname='pore.volume',
                       model=gm.pore_volume.sphere)
        self.add_model(propname='throat.length',
                       model=gm.throat_length.straight)
        self.add_model(propname='throat.diameter',
                       model=gm.throat_size.neighbor,
                       pore_prop='pore.diameter', mode='min')
        self.add_model(propname='throat.volume',
                       model=gm.throat_volume.cylinder)

        self.regenerate_models()
