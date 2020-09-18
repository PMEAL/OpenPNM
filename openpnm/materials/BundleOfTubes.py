import numpy as np
from openpnm.network import Cubic
from openpnm.utils import logging, Project, GenericSettings
from openpnm.geometry import GenericGeometry
from openpnm.topotools import trim
import openpnm.models as mods
logger = logging.getLogger(__name__)
defsets = {'adjust_psd': 'clip'}


class BundleOfTubesSettings(GenericSettings):
    r"""
    The following parameters as used when generating the network and geometry

    Parameters
    ----------
    distribution : str
        The type of statistical distriubtion to use.  Must correspond to one
        of the Scipy.stats options.
    loc : float
        The loc to pass to the distribution used, so refer to its help for
        more information
    scale : float
        The scale to pass to the distribution used, so refer to its help for
        more information
    scale : float
        The scale to pass to the distribution used, so refer to its help for
        more information
    adjust_psd : string {'normalize' (default), 'clip'}
        How the statistical distribution should be adjust in the event that
        it produces values larger than the spacing.
    seed : int (default is ``None``)
        The seed value to use in the random number generator. Setting this to
        a constant value should reproduce the same model on each run.
    """
    distribution = "norm"
    loc = 0.0
    scale = 0.4
    shape = 2.2
    adjust_psd = 'normalize'
    seed = None


class BundleOfTubes(Project):
    r"""
    The materials class creats a network that matches the bundle-of-tubes model

    An OpenPNM project object is returned that contain a network with a
    bundle-of-tubes topology, and geometry object with the necessary pore
    size information.

    Parameters
    ----------
    shape : array_like or int
        The number of pores in the X and Y direction of the domain.  It will
        be 2 pores thick in the Z direction, as required to be a bundle of
        tubes model (1 pore on each end of each throat makes a tube).  If an
        ``int`` is given it will be applied to both the X and Y directions.
    spacing : array_like or float
        The spacing between the tubes in the X and Y direction.  If a ``float``
        is given it will be applied to both the X and Y directions.
    length : float
        The length of the tubes or thickness of the domain in the z-direction.
    name : string, optional
        The name to give the Project

    """

    def __init__(self, shape, spacing=1.0, length=1.0, settings={},
                 name=None, **kwargs):
        import scipy.stats as spst

        super().__init__()
        self.settings._update_settings_and_docs(BundleOfTubesSettings())
        self.settings.update(settings)

        if isinstance(shape, int):
            shape = np.array([shape, shape, 2])
        elif len(shape) == 2:
            shape = np.concatenate((np.array(shape), [2]))
        else:
            raise Exception('shape not understood, must be int '
                            + ' or list of 2 ints')

        if isinstance(spacing, float) or isinstance(spacing, int):
            spacing = float(spacing)
            self.settings['spacing'] = spacing
            spacing = np.array([spacing, spacing, length])
        else:
            raise Exception('spacing not understood, must be float')

        net = Cubic(shape=shape, spacing=spacing, project=self, **kwargs)
        Ps_top = net.pores('top')
        Ps_bot = net.pores('bottom')
        Ts = net.find_connecting_throat(P1=Ps_top, P2=Ps_bot)
        Ts = net.tomask(throats=Ts)
        trim(network=net, throats=~Ts)

        geom = GenericGeometry(network=net, pores=net.Ps, throats=net.Ts)

        geom.add_model(propname='throat.seed',
                       model=mods.geometry.throat_seed.random)
        if self.settings['loc'] is None:
            self.settings['loc'] = spacing[0]/2
        if self.settings['scale'] is None:
            self.settings['scale'] = spacing[0]/10
        if self.settings['distribution'] in ['norm', 'normal', 'gaussian']:
            geom.add_model(propname='throat.size_distribution',
                           seeds='throat.seed',
                           model=mods.geometry.throat_size.normal,
                           loc=self.settings['loc'],
                           scale=self.settings['scale'])
        elif self.settings['distribution'] in ['weibull']:
            geom.add_model(propname='throat.size_distribution',
                           seeds='throat.seed',
                           model=mods.geometry.throat_size.weibull,
                           loc=self.settings['loc'],
                           scale=self.settings['scale'],
                           shape=self.settings['shape'])
        else:
            temp = self.settings.copy()
            func = getattr(spst, temp.pop('distribution'))
            psd = func.freeze(**temp)
            geom.add_model(propname='throat.size_distribution',
                           seeds='throat.seed',
                           model=mods.geometry.throat_size.generic_distribution,
                           func=psd)

        if np.any(geom['throat.size_distribution'] < 0):
            logger.warning('Given size distribution produced negative '
                           + 'tube diameters...these will be set to 0')
        geom.add_model(propname='throat.diameter',
                       model=mods.misc.clip,
                       prop='throat.size_distribution',
                       xmin=1e-12, xmax=np.inf)

        if self.settings['adjust_psd'] is None:
            if geom['throat.size_distribution'].max() > spacing[0]:
                logger.warning('Given size distribution produced tubes '
                               + 'larger than the spacing.')

        elif self.settings['adjust_psd'] == 'clip':
            geom.add_model(propname='throat.diameter',
                           model=mods.misc.clip,
                           prop='throat.size_distribution',
                           xmin=1e-12, xmax=spacing[0])
            if geom['throat.size_distribution'].max() > spacing[0]:
                logger.warning('Given size distribution produced tubes '
                               + 'larger than the spacing...tube diameters '
                               + 'will be clipped between 0 and given spacing')

        elif self.settings['adjust_psd'] == 'normalize':
            tmin = max(1e-12, geom['throat.size_distribution'].min())
            geom.add_model(propname='throat.diameter',
                           model=mods.misc.normalize,
                           prop='throat.size_distribution',
                           xmin=tmin, xmax=spacing[0])
            if geom['throat.size_distribution'].max() > spacing[0]:
                logger.warning('Given size distribution produced tubes '
                               + 'larger than the spacing...tube diameters '
                               + 'will be normalized to fit given spacing')
        else:
            logger.warning('Settings not understood, ignoring')

        geom.add_model(propname='pore.diameter',
                       model=mods.geometry.pore_size.from_neighbor_throats,
                       throat_prop='throat.diameter', mode='max')
        geom.add_model(propname='pore.diameter',
                       model=mods.misc.constant, value=0.0)
        geom.add_model(propname='throat.length',
                       model=mods.geometry.throat_length.ctc)
        geom.add_model(propname='throat.area',
                       model=mods.geometry.throat_area.cylinder)
        geom.add_model(propname='pore.area',
                       model=mods.misc.from_neighbor_throats,
                       throat_prop='throat.area')
        geom.add_model(propname='pore.volume',
                       model=mods.misc.constant, value=0.0)
        geom.add_model(propname='throat.volume',
                       model=mods.geometry.throat_volume.cylinder)

        geom.regenerate_models()
