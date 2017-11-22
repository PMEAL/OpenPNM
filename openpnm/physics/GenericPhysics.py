from openpnm.core import Base, Workspace, logging, utils
import scipy as sp
logger = logging.getLogger(__name__)
ws = Workspace()


class GenericPhysics(Base):
    r"""
    Generic class to generate Physics objects

    Parameters
    ----------
    network : OpenPNM Network object
        The network to which this Physics should be attached

    phase : OpenPNM Phase object
        The Phase object to which this Physics applies

    geometry : OpenPNM Geometry object
        The Geometry object that defines the pores/throats where this Physics
        should be applied.  If this argument is supplied, then pores and
        throats cannot be specified.

    name : str, optional
        A unique string name to identify the Physics object, typically same as
        instance name but can be anything.  If left blank, and name will be
        generated that include the class name and a random string.

    """

    def __init__(self, network=None, phase=None, geometry=None, name=None):
        super().__init__(name=name, simulation=network.simulation)
        logger.name = self.name

        # Initialize a label dictionary in the associated phase and network
        phase['pore.'+self.name] = False
        phase['pore.'+self.name][network.pores(geometry.name)] = True
        phase['throat.'+self.name] = False
        phase['throat.'+self.name][network.throats(geometry.name)] = True
        self['pore.all'] = sp.ones(shape=sp.size(geometry.Ps), dtype=bool)
        self['throat.all'] = sp.ones(shape=sp.size(geometry.Ts), dtype=bool)
        network.simulation.add_physics(physics=self, geometry=geometry,
                                       phase=phase)
