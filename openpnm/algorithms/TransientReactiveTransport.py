import scipy as sp
import scipy.sparse as sprs
from openpnm.algorithms import ReactiveTransport, TransientTransport
from openpnm.core import logging
logger = logging.getLogger(__name__)


class TransientReactiveTransport(ReactiveTransport, TransientTransport):
    r"""
    """

    def __init__(self, **kwargs):
        self.settings.update({'dummy': None})
        super().__init__(**kwargs)

    def run(self, x, t):
        r"""

        """
        x = self._run_transient_reactive(x=x, t=t)

    def _run_transient_reactive(self, x, t):
        x = self._reactive_run(x0=0)
        self.update_A()
        self.update_b()
        if t < self.settings['time_final']:
            self._transient_reactive_run(x=x, t=t+self.settings['time_step'])
        return x
