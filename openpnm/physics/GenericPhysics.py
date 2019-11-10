from openpnm.core import Subdomain, ModelsMixin
from openpnm.network import GenericNetwork
from openpnm.phases import GenericPhase
from openpnm.geometry import GenericGeometry
from openpnm.utils import Workspace, logging
logger = logging.getLogger(__name__)
ws = Workspace()


class GenericPhysics(Subdomain, ModelsMixin):
    r"""
    This generic class is meant as a starter for custom Physics objects

    It produces a blank object with no pore-scale models attached.  Users can
    add models from the ``models`` module (or create their own).

    Parameters
    ----------
    network : OpenPNM Network object
        The network to which this Physics should be attached

    phase : OpenPNM Phase object
        The Phase object to which this Physics applies.  If nothing is supplied
        then a new phase object is created and associated with this physics.

    geometry : OpenPNM Geometry object
        The Geometry object that defines the pores/throats where this Physics
        should be applied.  If nothing is supplied, then there are three cases:

        (1) If no geometries currently exist on the project, then an empty
        ``GenericGeometry`` will be created and assocated with this physics

        (2) If only one geometry is found on the project, then it will be
        associated with this physics

        (3) If more than one geometry is found then an Exception is raised
        since it's not clear which geometry to associate with this physics

    name : str, optional
        A unique string name to identify the Physics object, typically same as
        instance name but can be anything.  If left blank, and name will be
        generated that include the class name and a random string.

    """

    def __init__(self, project=None, network=None, phase=None,
                 geometry=None, settings={}, **kwargs):

        # Define some default settings
        self.settings.update({'prefix': 'phys'})
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
        if phase is None:
            logger.warning('No phase provided...a GenericPhase'
                           + ' will be created')
            phase = GenericPhase(network=network)
        self.set_phase(phase=phase)

        if geometry is None:
            g = list(self.project.geometries().values())
            if len(g) == 0:
                logger.warning('No geometry provided...a new GenericGeometry'
                               + ' will be created')
                geometry = GenericGeometry(network=network,
                                           pores=network.Ps,
                                           throats=network.Ts)
            elif len(g) == 1:
                if len(self.project.find_physics(geometry=g[0])):
                    raise Exception('No geometry provided, and the only'
                                    + ' existing geometry is already assigned'
                                    + ' to another physics')
                logger.warning('No geometry provided... so ' + self.name
                               + ' will be associated with the only available'
                               + ' geometry')
                geometry = g[0]
            elif len(g) > 1:
                raise Exception('No geometry provided, but multiple '
                                + ' geometries exist on project.  Please '
                                + ' specify which one.')
        elif isinstance(geometry, GenericNetwork):
            logger.warning('No geometry provided...a new GenericGeometry'
                           + ' will be created')
            geometry = GenericGeometry(network=network,
                                       pores=network.Ps,
                                       throats=network.Ts)
        self.set_geometry(geometry=geometry)

    def set_phase(self, phase, mode='add'):
        r"""
        """
        if phase not in self.project:
            raise Exception(self.name + ' not in same project as given phase')
        if mode == 'add':
            ok = False
            try:
                self.project.find_phase(self)
            except Exception:
                ok = True
            if ok:
                phase['pore.'+self.name] = False
                phase['throat.'+self.name] = False
            else:
                raise Exception('Physics is already associated with a phase'
                                + ' , use mode \'swap\' to change phases')
        elif mode == 'remove':
            current_phase = self.project.find_phase(self)
            if current_phase is not phase:
                raise Exception('Physics cannot be removed from a phase it'
                                + ' is not already associated with')
            phase.pop('pore.'+self.name, None)
            phase.pop('throat.'+self.name, None)
        elif mode == 'swap':
            current_phase = self.project.find_phase(self)
            self.set_phase(phase=current_phase, mode='remove')
            self.set_phase(phase=phase, mode='add')
        else:
            raise Exception('mode ' + mode + ' not understood')

    def set_geometry(self, geometry):
        r"""
        """
        if geometry not in self.project:
            raise Exception(self.name + ' not in same project as given geometry')
        network = self.network
        Ps = network.pores(geometry.name)
        Ts = network.throats(geometry.name)
        self._drop_locations(pores=Ps, throats=Ts)
        self._add_locations(pores=Ps, throats=Ts)
