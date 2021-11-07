from openpnm.core import Subdomain, ModelsMixin, ParamMixin
from openpnm.utils import Docorator
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


class GenericGeometry(ParamMixin, Subdomain, ModelsMixin):
    r"""
    This generic class is meant as a starter for custom Geometry objects

    It has no pore-scale models assigned to it, so a a blank slate.  Note that
    all OpenPNM Geometry sub-classes are just GenericGeometry instances with a
    number of models added.

    Parameters
    ----------
    network : GenericNetwork
        The Network object to which this Geometry applies.
    pores : array_like
        The list of pores where this Geometry applies.
    throats : array_like
        The list of throats where this Geometry applies.
    name : str
        A unique name to apply to the object.  This name will also be used as a
        label to identify where this Geometry applies.
    project : Project, optional
        A Project can be specified instead of ``network``.

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[5, 5, 5])
    >>> Ps = pn.pores('all')    # Get all pores
    >>> Ts = pn.throats('all')  # Get all throats
    >>> geom = op.geometry.GenericGeometry(network=pn, pores=Ps, throats=Ts)

    Now assign pore-scale models to the empty object:

    >>> geom.add_model(propname='pore.size',
    ...                model=op.models.misc.random,
    ...                element='pore',
    ...                num_range=[0.01, 0.1])

    Confirm that the object has one added model:

    >>> print(geom.models)
    ―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
    #   Property Name                       Parameter                 Value
    ―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
    1   pore.size                           model:                    random
                                            element:                  pore
                                            num_range:                [0.01, 0.1]
                                            seed:                     None
                                            regeneration mode:        normal
    ―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――

    The results of the model can be seen using the ``show_hist`` function:

    >>> import matplotlib as mpl
    >>> mpl.use('Agg')
    >>> geom.show_hist('pore.size')

    .. image:: /../docs/_static/images/generic_geometry_histogram.png
        :width: 500px
        :align: center

    """

    def __init__(self, pores=[], throats=[], settings={}, **kwargs):
        self.settings._update(GeometrySettings, docs=True)
        self.settings._update(settings)  # Add user supplied settings
        super().__init__(**kwargs)

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
