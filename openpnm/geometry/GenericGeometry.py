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

    def __init__(self, network=None, project=None, pores=[], throats=[],
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
            network[f'pore.{self.name}'] = False
            network[f'throat.{self.name}'] = False
            try:
                self.add_locations(pores=pores, throats=throats)
            except Exception as e:
                network.project.purge_object(self)
                logger.error(f'{e}, instantiation cancelled')

    def add_locations(self, pores=[], throats=[]):
        r"""
        Adds associations between this geometry and the given pore and/or
        throat locations.

        Parameters
        ----------
        pores and throats : array_like
            The pore and/or throat locations for which the association should
            be added.  These indices are for the full domain.

        Notes
        -----
        If a physics object is associated with this geometry, then its
        pore and/or throat associations are also changed.
        """
        pores = self.network._parse_indices(pores)
        throats = self.network._parse_indices(throats)
        objects = self.project.find_physics(self)
        objects.append(self)
        for obj in objects:
            if len(pores) > 0:
                obj._set_locations(element='pore', indices=pores, mode='add')
            if len(throats) > 0:
                obj._set_locations(element='throat', indices=throats, mode='add')

    def drop_locations(self, pores=[], throats=[]):
        r"""
        Removes association between this geometry and the given pore and/or
        throat locations.

        Parameters
        ----------
        pores and throats : array_like
            The pore and/or throat locations from which the association should
            be removed.  These indices refer to the full domain.

        Notes
        -----
        If a physics object is associated with this geometry, then its
        pore and/or throat associations are also changed.

        """
        pores = self.network._parse_indices(pores)
        throats = self.network._parse_indices(throats)
        objects = self.project.find_physics(self)
        objects.append(self)
        for obj in objects:
            if len(pores) > 0:
                obj._set_locations(element='pore', indices=pores, mode='drop')
            if len(throats) > 0:
                obj._set_locations(element='throat', indices=throats, mode='drop')
