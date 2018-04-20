from openpnm.algorithms import StokesFlow
from openpnm.algorithms import TransientTransport
from openpnm.core import logging
logger = logging.getLogger(__name__)


class TransientStokesFlow(TransientTransport, StokesFlow):
    r"""
    A subclass of GenericTransport to simulate Stokes flow.

    """

    def __init__(self, settings={}, **kwargs):
        super().__init__(**kwargs)
        # Apply any received settings to overwrite defaults
        self.settings.update(settings)
        # Save A matrix of the steady sys of eqs
        self.A = self._build_A()
        self._coef = 1
