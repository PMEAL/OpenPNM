from openpnm.core import Base, ModelsMixin
from openpnm.utils import PrintableDict, Workspace, logging
import openpnm.models as mods
logger = logging.getLogger(__name__)
ws = Workspace()
from numpy import ones


class GenericPhase(Base, ModelsMixin):
    r"""
    This generic class is meant as a starter for custom Phase objects

    This class produces a blank-slate object with no pore-scale models for
    calculating any thermophysical properties.  Users must add models and
    specify parameters for all the properties they require.

    Parameters
    ----------
    network : openpnm Network object
        The network to which this Phase should be attached

    name : str, optional
        A unique string name to identify the Phase object, typically same as
        instance name but can be anything.

    Examples
    --------
    Create a new empty phase:

    >>> import openpnm as op
    >>> pn = op.network.Cubic([10, 10, 10])
    >>> phase = op.phases.GenericPhase(network=pn)

    And add a model:

    >>> phase.add_model(propname='pore.molar_density',
    ...                 model=op.models.phases.molar_density.ideal_gas)

    Now confirm that the model was added and data was calculated.  The
    ``models`` attribute can be printed:

    >>> print(phase.models)
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
    #   Property Name             Parameter                 Value
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
    1   pore.molar_density        model:                    ideal_gas
                                  pressure:                 pore.pressure
                                  temperature:              pore.temperature
                                  regeneration mode:        normal
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――

    And the Phase itself has a nice printout using ``print(phase)``.

    """

    def __init__(self, network=None, project=None, settings={}, **kwargs):
        # Define some default settings
        self.settings.update({'prefix': 'phase'})
        # Overwrite with user supplied settings, if any
        self.settings.update(settings)

        # Deal with network or project arguments
        if network is not None:
            if project is not None:
                assert network is project.network
            else:
                project = network.project

        super().__init__(project=project, **kwargs)

        # If project has a network object, adjust pore and throat sizes
        network = self.project.network
        if network:
            self['pore.all'] = ones((network.Np, ), dtype=bool)
            self['throat.all'] = ones((network.Nt, ), dtype=bool)

        # Set standard conditions on the fluid to get started
        self['pore.temperature'] = 298.0
        self['pore.pressure'] = 101325.0

    def __getitem__(self, key):
        element, prop = key.split('.', 1)
        # Deal with special keys first
        if prop == '_id':
            net = self.project.network
            return net[element+'._id']
        if prop == self.name:
            return self[element+'.all']
        # Now get values if present, or regenerate them
        vals = self.get(key)
        if vals is None:
            if element == 'throat' and 'pore.'+prop in self.keys():
                vals = self.interpolate_data(propname='pore.'+prop)
                logger.info(key + ', not found, interpolating from ' +
                            'pore.' + prop)
            elif element == 'pore' and 'throat.'+prop in self.keys():
                vals = self.interpolate_data(propname='throat.'+prop)
                logger.info(key + ', not found, interpolating from ' +
                            'throat.' + prop)
            else:
                vals = self.interleave_data(key)
        return vals
