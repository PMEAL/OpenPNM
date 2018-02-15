from openpnm.core import logging, Simulation
from openpnm.network import Cubic
from openpnm.geometry import GenericGeometry
import openpnm.models.geometry as gm
logger = logging.getLogger(__name__)


class BereaCubic(Simulation):

    def __init__(self, shape, name='BereaCubic', **kwargs):
        super().__init__(name=name)

        net = Cubic(shape=shape, spacing=1e-4, connectivity=26,
                    simulation=self, **kwargs)

        geom = GenericGeometry(network=net, pores=net.Ps, throats=net.Ts)

        geom.add_model(propname='pore.seed',
                       model=gm.pore_seed.spatially_correlated,
                       weights=[1, 1, 1])
        geom.add_model(propname='pore.diameter',
                       model=gm.pore_size.weibull,
                       shape=2.2, scale=2e-5, loc=1e-6)
        geom.add_model(propname='pore.volume',
                       model=gm.pore_volume.sphere)
        geom.add_model(propname='throat.length',
                       model=gm.throat_length.straight)
        geom.add_model(propname='throat.diameter',
                       model=gm.throat_size.neighbor,
                       pore_prop='pore.diameter', mode='min')
        geom.add_model(propname='throat.volume',
                       model=gm.throat_volume.cylinder)

        geom.regenerate_models()
