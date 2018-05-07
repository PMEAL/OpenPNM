from openpnm.algorithms import OrdinaryPercolation
from openpnm.core import logging
logger = logging.getLogger()

default_settings = {'pore_volume': 'pore.volume',
                    'throat_volume': 'throat.volume'}


class Porosimetry(OrdinaryPercolation):

    def __init__(self, settings={}, **kwargs):
        super().__init__(**kwargs)
        self.settings.update(default_settings)
        # Apply user settings, if any
        self.settings.update(settings)
        # Use the reset method to initialize all arrays
        self.reset()

    def reset(self):
        super().reset()

    def setup(self, **kwargs):
        super().setup(**kwargs)
