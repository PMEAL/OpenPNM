from openpnm.physics import GenericPhysics
from openpnm.models import physics as mods

__all__ = ["Basic"]


class Basic(GenericPhysics):
    r"""
    Minimal subclass of GenericPhysics for performing diffusion and/or
    flow simulations.

    Parameters
    ----------
    network : GenericNetwork
        The network to which this Physics should be attached/

    phase : GenericPhase
        The Phase object to which this Physics applies.

    geometry : GenericGeometry
        The Geometry object that defines the pores/throats where this
        Physics should be applied.

    name : str, optional
        A unique string name to identify the Physics object, typically
        same as instance name but can be anything.  If left blank, and
        name will be generated that include the class name and a random
        string.

    Returns
    -------
    GenericPhysics
        Basic subclass of GenericPhysics with minimal models attached
        required to perform diffusion and/or flow simulations.

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

        self.add_model(propname='throat.hydraulic_conductance',
                       model=mods.hydraulic_conductance.generic_hydraulic)
        self.add_model(propname='throat.diffusive_conductance',
                       model=mods.diffusive_conductance.generic_diffusive)
        self.add_model(propname='throat.entry_pressure',
                       model=mods.capillary_pressure.washburn)
