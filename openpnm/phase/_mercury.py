from openpnm.models.collections.phase import mercury
from openpnm.phases import GenericPhase
import openpnm.models as mods


class Mercury(GenericPhase):
    r"""
    Creates Phase object with a default name 'Hg' and preset values and
    pore-scale models for mercury.

    Parameters
    ----------
    %(GenericPhase.parameters)s

    References
    ----------
    The correlations and constants for this class were taken from:

    ::

        Thermophysical Properties of Materials for Nuclear Engineering:
        IAEA, Vienna, 2008. ISBN 978-92-0-106508-7:

    """
    def __init__(self, name=None, **kwargs):
        super().__init__(name=name, **kwargs)
        self.models.update(mercury)
        self.regenerate_models()
