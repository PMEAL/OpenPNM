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

    def update_A(self):
        r"""
        """
        pass

    def update_b(self):
        r"""
        """
        pass

    def run(self, t=0):
        r"""
        """
        self.setup()
        self._run_transient(t=t)

    def _run_transient(self, t):
        self.update_A()
        self.update_b()
        x_new = self._solve()
        self[self.settings['quantity']] = x_new
        if t < self.settings['t_final']:
            self._run_transient(t=t + self.settings['t_step'])
        else:
            print('Maximum time step reached')
