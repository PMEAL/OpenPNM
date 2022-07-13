from openpnm.models.collections.phase import mercury
from openpnm.phase import Phase


class Mercury(Phase):
    r"""
    Creates Phase object with and preset values and pore-scale models for
    mercury

    Parameters
    ----------
    %(Phase.parameters)s

    References
    ----------
    The correlations and constants for this class were taken from:

    ::

        Thermophysical Properties of Materials for Nuclear Engineering:
        IAEA, Vienna, 2008. ISBN 978-92-0-106508-7:

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.add_model_collection(mercury())
        self.regenerate_models()
