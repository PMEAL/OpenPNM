from openpnm.phases import mixtures
import openpnm.models as mods


class H2(mixtures.GenericSpecies):
    r"""
    Creates Phase object with preset models and values for Na ions

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
    >>> H2 = mixtures.species.gases.H2(network=pn)

    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self['pore.molecular_weight'] = 0.002016  # kg/mol
        self['pore.molar_diffusion_volume'] = 17.9  # Wrong
