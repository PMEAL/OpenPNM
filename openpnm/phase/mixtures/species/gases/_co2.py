from openpnm.phase import Phase
import openpnm.models as mods


class CO2(Phase):
    r"""
    Creates Phase object with preset models and values for CO2 gas

    Parameters
    ----------
    network : GenericNetwork
        The network to which this phase object will be attached.
    name : str, optional
        The name of the phase.  This is useful to keep track of the objects
        throughout the simulation.  The name must be unique to the project.  If
        no name is given, one is generated.

    Examples
    --------
    >>> import openpnm as op
    >>> import openpnm.phase.mixtures as mixtures
    >>> pn = op.network.Cubic(shape=[5, 5, 5])
    >>> CO2 = mixtures.species.gases.CO2(network=pn)

    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self['pore.molecular_weight'] = 0.04401  # kg/mol
        self['pore.molar_diffusion_volume'] = 17.9  # Wrong
