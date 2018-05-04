from openpnm.algorithms import TransientReactiveTransport, StokesFlow
from openpnm.core import logging
logger = logging.getLogger(__name__)


class TransientStokesFlow(TransientReactiveTransport, StokesFlow):
    r"""
    A subclass of GenericTransport to simulate Stokes flow.

    """

    def __init__(self, settings={}, **kwargs):
        super().__init__(**kwargs)
        # Apply any received settings to overwrite defaults
        self.settings.update(settings)
