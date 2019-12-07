from openpnm.core import Subdomain, ModelsMixin
from openpnm.utils import Workspace, logging
logger = logging.getLogger(__name__)
ws = Workspace()


class GenericGeometry(Subdomain, ModelsMixin):
    r"""
    This generic class is meant as a starter for custom Geometry objects

    It has no pore-scale models assigned to it, so a a blank slate.  Note that
    all OpenPNM Geometry sub-classes are just GenericGeometry instances with a
    number of models added.

    Parameters
    ----------
    network : OpenPNM Network Object
        The Network object to which this Geometry applies.

    pores : array_like
        The list of pores where this Geometry applies.

    throats : array_like
        The list of throats where this Geometry applies.

    name : string
        A unique name to apply to the object.  This name will also be used as a
        label to identify where this Geometry applies.

    project : OpenPNM Project object (optional)
        A Project can be specified instead of ``network``.

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[5, 5, 5])
    >>> Ps = pn.pores('all')  # Get all pores
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

    .. image:: /../docs/static/images/generic_geometry_histogram.png
        :width: 500px
        :align: center

    """

    def __init__(self, network=None, project=None, pores=None, throats=None,
                 settings={}, **kwargs):
        # Define some default settings
        self.settings.update({'prefix': 'geo'})
        # Overwrite with user supplied settings, if any
        self.settings.update(settings)

        # Deal with network or project arguments
        if network is not None:
            if project is not None:
                assert network is project.network
            else:
                project = network.project

        super().__init__(project=project, **kwargs)

        network = self.project.network
        if network:
            network['pore.'+self.name] = False
            network['throat.'+self.name] = False
            if (pores is None) and (throats is None):
                logger.info('No pores and throats given, assigning '
                            + self.name + ' to entire domain')
                pores = network.Ps
                throats = network.Ts
            try:
                self._add_locations(pores=pores, throats=throats)
            except Exception as e:
                network.project.purge_object(self)
                logger.error(str(e) +  ', instantiation cancelled')
