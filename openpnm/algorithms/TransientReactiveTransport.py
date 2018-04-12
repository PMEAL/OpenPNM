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
        x = self._run_transient_reactive(x=x, t=t)

    def _run_transient_reactive(self, x, t):
        self._apply_BCs()
        self._update_A()
        self._update_b()
        x = self._run_reactive(x=x)
        if t < self.settings['t_final']:
            print('Current time step: '+str(t))
            self._run_transient_reactive(x=x, t=t+self.settings['t_step'])
        else:
            print('Maximum time step reached')
