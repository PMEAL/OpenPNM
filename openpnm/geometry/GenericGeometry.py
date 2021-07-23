from openpnm.core import Subdomain, ModelsMixin
from openpnm.utils import Workspace, Docorator, GenericSettings, logging
ws = Workspace()
docstr = Docorator()
logger = logging.getLogger(__name__)


@docstr.get_sections(base='GenericGeometry',
                     sections=['Parameters'])
@docstr.dedent
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
    project : Project, optional
        A Project can be specified instead of ``network``.
    pores : array_like
        The list of pores where this Geometry applies.
    throats : array_like
        The list of throats where this Geometry applies.
    name : str
        A unique name to apply to the object.  This name will also be used as a
        label to identify where this Geometry applies.

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
