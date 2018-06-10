from openpnm.algorithms import TransientReactiveTransport, Dispersion
from openpnm.core import logging
logger = logging.getLogger(__name__)


class TransientDispersion(TransientReactiveTransport, Dispersion):
    r"""

    """

    def __init__(self, settings={}, **kwargs):
        super().__init__(**kwargs)
        # Apply any received settings to overwrite defaults
        self.settings.update(settings)
