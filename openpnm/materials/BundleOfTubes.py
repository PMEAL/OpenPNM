import numpy as np
from openpnm.network import Cubic
from openpnm.utils import logging, Project
from openpnm.geometry import GenericGeometry
from openpnm.phases import GenericPhase
from openpnm.topotools import trim
import openpnm.models as mods
logger = logging.getLogger(__name__)
defsets = {'adjust_psd': 'clip'}


class BundleOfTubes(Project):
    r"""
    The materials class creats a network that matches the bundle-of-tubes model

    An OpenPNM project object is returned that contain a network with a
    bundle-of-tubes topology, and geometry object with the necessary pore
    size information, and a phase object with pre-defined pore-scale physics
    models attached.  Note that this phase object does not have any actual
    thermophysical properties which must be added by the user.

    Parameters
    ----------
    shape : array_like or int
        The number of pores in the X and Y direction of the domain.  It will
        be 2 pores thick in the Z direction, as required to be a bundle of
        tubes model.  If an ``int`` is given it will be applied to both the
        X and Y directions.

    spacing : array_like or float
        The spacing between the tubes in the X and Y direction.  If a ``float``
        is given it will be applied to both the X and Y directions.

    length : float
        The length of the tubes or thickness of the domain in the z-direction.

    psd_params : dictionary
        The parameters of the statistical distribution of the pore sizes.
        The dictionary must contain the type of distribution to use (specified
        as 'distribution'), selected from the ``scipy.stats`` module, and the
        parameters corresponding to the chosen distribution.
        The default is: ``{'distribution': 'norm', 'loc': 0.5, 'scale': 0.11}.
        Note that ``scipy.stats`` uses *loc* and *scale* to be consistent
        between different distribution types, instead of things like *mean*
        and *stddev*.

    name : string, optional
        The name to give the Project

    """
    def __init__(
        self,
        shape,
        spacing=1.0,
        length=1.0,
        psd_params={"distribution": "norm", "loc": None, "scale": None},
        name=None,
        settings={},
        **kwargs
    ):
        import scipy.stats as spst

        super().__init__(name=name)
        self.settings.update(defsets)
        self.settings.update(settings)

        if isinstance(shape, int):
            shape = np.array([shape, shape, 2])
        elif len(shape) == 2:
            shape = np.concatenate((np.array(shape), [2]))
        else:
            raise Exception('shape not understood, must be int '
                            + ' or list of 2 ints')

        if isinstance(spacing, (float, int)):
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
        if psd_params['loc'] is None:
            psd_params['loc'] = spacing[0]/2
        if psd_params['scale'] is None:
            psd_params['scale'] = spacing[0]/10
        if psd_params['distribution'] in ['norm', 'normal', 'gaussian']:
            geom.add_model(propname='throat.size_distribution',
                           seeds='throat.seed',
                           model=mods.geometry.throat_size.normal,
                           loc=psd_params['loc'], scale=psd_params['scale'])
        elif psd_params['distribution'] in ['weibull']:
            geom.add_model(propname='throat.size_distribution',
                           seeds='throat.seed',
                           model=mods.geometry.throat_size.weibull,
                           loc=psd_params['loc'],
                           scale=psd_params['scale'],
                           shape=psd_params['shape'])
        else:
            temp = psd_params.copy()
            func = getattr(spst, temp.pop('distribution'))
            psd = func.freeze(**temp)
            geom.add_model(propname='throat.size_distribution',
                           seeds='throat.seed',
                           model=mods.geometry.throat_size.generic_distribution,
                           func=psd)

        if np.any(geom['throat.size_distribution'] < 0):
            logger.warning('Given size distribution produced negative '
                           + 'throat diameters...these will be set to 0')
        geom.add_model(propname='throat.diameter',
                       model=mods.misc.clip,
                       prop='throat.size_distribution',
                       xmin=1e-12, xmax=np.inf)

        if self.settings['adjust_psd'] is None:
            if geom['throat.size_distribution'].max() > spacing[0]:
                logger.warning('Given size distribution produced throats '
                               + 'larger than the spacing.')

        elif self.settings['adjust_psd'] == 'clip':
            geom.add_model(propname='throat.diameter',
                           model=mods.misc.clip,
                           prop='throat.size_distribution',
                           xmin=1e-12, xmax=spacing[0])
            if geom['throat.size_distribution'].max() > spacing[0]:
                logger.warning('Given size distribution produced throats '
                               + 'larger than the spacing...tube diameters '
                               + 'will be clipped between 0 and given spacing')

        elif self.settings['adjust_psd'] == 'normalize':
            tmin = max(1e-12, geom['throat.size_distribution'].min())
            geom.add_model(propname='throat.diameter',
                           model=mods.misc.normalize,
                           prop='throat.size_distribution',
                           xmin=tmin, xmax=spacing[0])
            if geom['throat.size_distribution'].max() > spacing[0]:
                logger.warning('Given size distribution produced throats '
                               + 'larger than the spacing...tube diameters '
                               + 'will be normalized to fit given spacing')
        else:
            logger.warning('Settings not understood, ignoring')

        geom.add_model(propname='pore.diameter',
                       model=mods.geometry.pore_size.from_neighbor_throats,
                       prop='throat.diameter', mode='max')
        geom.add_model(propname='pore.diameter',
                       model=mods.misc.constant, value=0.0)
        geom.add_model(propname='throat.length',
                       model=mods.geometry.throat_length.ctc)
        geom.add_model(propname='throat.area',
                       model=mods.geometry.throat_cross_sectional_area.cylinder)
        geom.add_model(propname='pore.area',
                       model=mods.misc.from_neighbor_throats,
                       prop='throat.area')
        geom.add_model(propname='pore.volume',
                       model=mods.misc.constant, value=0.0)
        geom.add_model(propname='throat.volume',
                       model=mods.geometry.throat_volume.cylinder)

        geom.regenerate_models()

        # Now create a generic phase with physics models on it
        phase = GenericPhase(network=net)
        m = mods.physics.hydraulic_conductance.classic_hagen_poiseuille
        phase.add_model(propname='throat.hydraulic_conductance',
                        model=m, regen_mode='deferred')
        m = mods.physics.diffusive_conductance.classic_ordinary_diffusion
        phase.add_model(propname='throat.diffusive_conductance',
                        model=m, regen_mode='deferred')
        m = mods.physics.diffusive_conductance.classic_ordinary_diffusion
        phase.add_model(propname='throat.entry_pressure',
                        model=m, regen_mode='deferred')
