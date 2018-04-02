from openpnm.algorithms import GenericTransport
from openpnm.core import logging
logger = logging.getLogger(__name__)


class TransientTransport(GenericTransport):
    r"""
    """

    def __init__(self, **kwargs):
        self.settings.update({'t_initial': 0,
                              't_final': None,
                              't_step': None})
        super().__init__(**kwargs)

    def set_IC(self, values):
        r"""
        """
        self[self.settings['quantity']] = values

    def _update_A(self):
        r"""
        """
        pass

    def _update_b(self):
        r"""
        """
        pass

    def run(self, t=0):
        r"""
        """
        print('―'*80)
        print('Running TransientTransport')
        self._run_transient(t=t)

    def _run_transient(self, t):
        self._apply_BCs()
        self._update_A()
        self._update_b()
        x_new = self._solve()
        self[self.settings['quantity']] = x_new
        if t < self.settings['t_final']:
            print('Current time step: '+str(t))
            self._run_transient(t=t + self.settings['t_step'])
        else:
            print('Maximum time step reached')
