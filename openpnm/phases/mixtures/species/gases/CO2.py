from openpnm.phases import GenericPhase
import openpnm.models as mods


class CO2(GenericPhase):
    r"""
    Creates Phase object with preset models and values for CO2 gas

    Parameters
    ----------
    network : OpenPNM Network object
        The network to which this phase object will be attached.

    project : OpenPNM Project object, optional
        The Project with which this phase should be associted.  If a
        ``network`` is given then this is ignored and the Network's project is
        used.  If a ``network`` is not given then this is mandatory.

    name : string, optional
        The name of the phase.  This is useful to keep track of the objects
        throughout the simulation.  The name must be unique to the project.  If
        no name is given, one is generated.

    Examples
    --------
    >>> import openpnm as op
    >>> import openpnm.phases.mixtures as mixtures
    >>> pn = op.network.Cubic(shape=[5, 5, 5])
    >>> CO2 = mixtures.species.gases.CO2(network=pn)

    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self['pore.molecular_weight'] = 0.04401  # kg/mol
        self['pore.molar_diffusion_volume'] = 17.9  # Wrong
