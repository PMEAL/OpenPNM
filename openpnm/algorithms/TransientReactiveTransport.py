from openpnm.algorithms import ReactiveTransport, TransientTransport
from openpnm.core import logging
logger = logging.getLogger(__name__)


class TransientReactiveTransport(ReactiveTransport, TransientTransport):
    r"""
    """

    def __init__(self, **kwargs):
        self.settings.update({'dummy': None})
        super().__init__(**kwargs)

    def run(self, x=None, t=0):
        r"""

        """
        print('―'*80)
        print('Running TransientReactiveTransport')
        self.setup()
        x = self._run_transient_reactive(x=x, t=t)

    def _run_transient_reactive(self, x, t):
        x = self._run_reactive(x=x)
        self.update_A()
        self.update_b()
        if t < self.settings['t_final']:
            self._run_transient_reactive(x=x, t=t+self.settings['t_step'])
        return x
