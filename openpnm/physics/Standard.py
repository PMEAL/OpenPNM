from openpnm.utils import Workspace, logging
from openpnm.physics import GenericPhysics
from openpnm.models import physics as mods
logger = logging.getLogger(__name__)
ws = Workspace()


class Standard(GenericPhysics):
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
        should be applied.

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

        super().__init__(project=project, phase=phase, geometry=geometry,
                         **kwargs)

        self.add_model(propname='throat.flow_shape_factors',
                       model=mods.flow_shape_factors.conical_frustum_and_stick)
        self.add_model(propname='throat.hydraulic_conductance',
                       model=mods.hydraulic_conductance.hagen_poiseuille)
        self.add_model(propname='throat.poisson_shape_factors',
                       model=mods.poisson_shape_factors.conical_frustum_and_stick)
        self.add_model(propname='throat.diffusive_conductance',
                       model=mods.diffusive_conductance.mixed_diffusion)
        self.add_model(propname='throat.ad_dif_conductance',
                       model=mods.ad_dif_conductance.ad_dif)
        self.add_model(propname='throat.entry_pressure',
                       model=mods.capillary_pressure.washburn)
        self.add_model(propname='throat.thermal_conductance',
                       model=mods.thermal_conductance.series_resistors)
        self.add_model(propname='throat.electrical_conductance',
                       model=mods.electrical_conductance.series_resistors)
