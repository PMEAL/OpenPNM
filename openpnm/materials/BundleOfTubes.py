import scipy as sp
import scipy.stats as spst
from openpnm.utils import logging, Project
from openpnm.network import Cubic
from openpnm.geometry import GenericGeometry
import openpnm.models.geometry as gm
logger = logging.getLogger(__name__)


class BundleOfTubes(Project):
    r"""
    A basic 'bundle-of-tubes' model.

    Parameters
    ----------
    shape : array_like
        The number of pores in the X and Y direction of the domain.  It will
        be 2 pores thick in the Z direction, as required to be a bundle of
        tubes model.

    psd_params : dictionary
        The parameters of the statistical distribution of the pore sizes.
        The dictionary must contain the type of distribution to use (specified
        as 'distribution'), selected from the ``scipy.stats`` module, and the
        parameters corresponding to the chosen distribution.
        The default is: ``{'distribution': 'norm', 'loc': 0.5, 'scale': 1}.
        Note that ``scipy.stats`` uses *loc* and *scale* to be consistent
        between different distribution types, instead of things like *mean*
        and *stddev*.

    name : string, optional
        The name to give the Project

    Examples
    --------

    """

    def __init__(self, shape, psd_params={'distribution': 'normal',
                                          'loc': 0.5, 'shape': 1},
                 name=None, **kwargs):
        super().__init__(name=name)

        if isinstance(shape, int):
            shape = sp.array([shape, shape, 1])
        elif len(shape) == 2:
            shape = sp.concatenate((sp.array(shape), [1]))
        else:
            raise Exception('shape not understood, must be int of list of 2 ints')

        net = Cubic(shape=shape, spacing=1, project=self, **kwargs)

        geom = GenericGeometry(network=net, pores=net.Ps, throats=net.Ts)

        geom.add_model(propname='pore.seed',
                       model=gm.pore_seed.random)
        func = getattr(spst, psd_params['distribution'])
        psd = func.freeze(loc=1, scale=0.5)
        geom.add_model(propname='pore.diameter',
                       seeds='pore.seed',
                       model=gm.pore_size.generic_distribution)
        geom.add_model(propname='pore.volume',
                       model=gm.pore_volume.sphere)
        geom.add_model(propname='throat.length',
                       model=gm.throat_length.piecewise)
        geom.add_model(propname='throat.diameter',
                       model=gm.throat_size.from_neighbor_pores,
                       pore_prop='pore.diameter', mode='min')
        geom.add_model(propname='throat.volume',
                       model=gm.throat_volume.cylinder)

        geom.regenerate_models()
