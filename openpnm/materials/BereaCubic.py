from openpnm.utils import logging, Project
from openpnm.network import Cubic
from openpnm.geometry import GenericGeometry
import openpnm.models.geometry as gm
logger = logging.getLogger(__name__)


class BereaCubic(Project):
    r"""
    A traditional Berea Sandstone on a Cubic lattice

    Berea Sandstone is one of the standard materials used on geoscience
    studies due to it's importance in oil reservoir engineering as well as
    having well defined pore structure.  This class creates a Cubic Network
    with the appropriate lattice spacing and connectivity, then adds a Geometry
    object with the necessary pore-scale models and prescribed parameters.

    Parameters
    ----------
    shape : array_like
        The number of pores along each direction of the domain.  All other
        aspects of this model are prescribed by the code.

    name : string, optional
        The name to give the Project

    Notes
    -----
    The source code for this Material is relatively straight-forward, so is a
    good example starting point for creating custom materials.

    References
    ----------
    [1] ???

    Examples
    --------

    """

    def __init__(self, shape, name=None, **kwargs):
        super().__init__(name=name)

        net = Cubic(shape=shape, spacing=1e-4, connectivity=26,
                    project=self, **kwargs)

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
                       model=gm.throat_length.piecewise)
        geom.add_model(propname='throat.diameter',
                       model=gm.throat_size.from_neighbor_pores,
                       prop='pore.diameter', mode='min')
        geom.add_model(propname='throat.volume',
                       model=gm.throat_volume.cylinder)

        geom.regenerate_models()
