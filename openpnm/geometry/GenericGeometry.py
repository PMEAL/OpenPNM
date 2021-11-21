from copy import deepcopy
from openpnm.core import Subdomain, ModelsMixin, ParamMixin
from openpnm.utils import Docorator, SettingsAttr
from openpnm.utils import Workspace, logging
logger = logging.getLogger(__name__)
ws = Workspace()
docstr = Docorator()


@docstr.get_sections(base='GeometrySettings', sections=['Parameters'])
@docstr.dedent
class GeometrySettings:
    r"""

    Parameters
    ----------
    %(BaseSettings.parameters)s
    """
    prefix = 'geo'


@docstr.get_sections(base='GenericGeometry', sections=['Parameters',
                                                       'Examples'])
@docstr.dedent
class GenericGeometry(ParamMixin, Subdomain, ModelsMixin):
    r"""
    This generic class is meant as a starter for custom Geometry objects

    It has no pore-scale models assigned to it, so is a blank slate.  Note that
    all OpenPNM Geometry sub-classes are just GenericGeometry instances with a
    assortment of models added.

    Parameters
    ----------
    pores : array_like
        The list of pores where this geometry applies.
    throats : array_like
        The list of throats where this geometry applies.
    %(Base.parameters)s


    """

    def __init__(self, pores=[], throats=[], settings={}, **kwargs):
        self.settings = SettingsAttr(GeometrySettings, settings)
        super().__init__(settings=self.settings, **kwargs)

        network = self.project.network
        if network:
            network[f'pore.{self.name}'] = False
            network[f'throat.{self.name}'] = False
            try:
                self.set_locations(pores=pores, throats=throats, mode='add')
            except Exception as e:
                network.project.purge_object(self)
                raise e
                logger.error(f'{e}, instantiation cancelled')

    def _get_phys(self):
        return self.project.find_physics(geometry=self)
    physics = property(fget=_get_phys)
