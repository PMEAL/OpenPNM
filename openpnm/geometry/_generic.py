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
        The list of throats where this Geometry applies.
    name : str
        A unique name to apply to the object.  This name will also be used as a
        label to identify where this Geometry applies.

    Examples
    --------
    .. plot::

       import openpnm as op
       import matplotlib.pyplot as plt

       pn = op.network.Cubic(shape=[5, 5, 5])
       Ps = pn.pores('all')    # Get all pores
       Ts = pn.throats('all')  # Get all throats
       geom = op.geometry.GenericGeometry(network=pn, pores=Ps, throats=Ts)

       # Now assign pore-scale models to the empty object
       geom.add_model(propname='pore.size',
                      model=op.models.misc.random,
                      element='pore',
                      num_range=[0.01, 0.1])

       # Confirm that the object has one added model
       print(geom.models)

       # The results of the model can be seen using the ``show_hist`` function:
       geom.show_hist('pore.size')

       plt.show()

    """

    def __init__(self, pores=[], throats=[], settings=None, **kwargs):
        self.settings = SettingsAttr(GeometrySettings, settings)
        super().__init__(settings=self.settings, **kwargs)

        network = self.project.network
        if network:
            network[f'pore.{self.name}'] = False
            network[f'throat.{self.name}'] = False
            try:
                self.set_locations(pores=pores, throats=throats, mode='add')
            except Exception as e:
                logger.error(f'{e}, instantiation cancelled')
                network.project.purge_object(self)
                raise e

    def _get_phys(self):
        """A shortcut to get a handle to the associated physics."""
        return self.project.find_physics(geometry=self)
    physics = property(fget=_get_phys)
