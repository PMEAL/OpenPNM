from openpnm.phase import mixtures
import openpnm.models as mods


class Cl(mixtures.GenericSpecies):
    r"""
    Creates Phase object with preset models and values for Cl ions

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
    >>> Cl = mixtures.species.ions.Cl(network=pn)

    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self['pore.molecular_weight'] = 0.03545  # kg/mol
        self['pore.diffusivity'] = 2.03e-09  # m2/s
        self['throat.diffusivity'] = 2.03e-09  # m2/s
        self['pore.valence'] = -1
        self['throat.valence'] = -1
