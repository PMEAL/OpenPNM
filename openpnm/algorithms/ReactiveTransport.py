import numpy as np
from openpnm.algorithms import GenericTransport
from openpnm.core import logging
logger = logging.getLogger(__name__)


class ReactiveTransport(GenericTransport):
    r"""
    """

    def __init__(self, **kwargs):
        self.settings.update({'sources': [],
                              'tolerance': 0.001})
        super().__init__(**kwargs)

    def set_source_term(self, propname, pores):
        r"""
        """
        self.settings['sources'].append(propname)
        self[propname] = self.tomask(pores=pores)

    def apply_sources(self):
        phase = self.project.phases()[self.settings['phase']]
        for item in self.settings['sources']:
            Ps = self.pores(item)
            # Regenerate models with new guess
            quantity = self.settings['quantity']
            # Put quantity on phase so physics finds it when regenerating
            phase[quantity] = self[quantity]
            phase.regenerate_models(propnames=item)
            # Add S1 to diagonal of A
            # TODO: We need this to NOT overwrite the A and b, but create
            # copy, otherwise we have to regenerate A and b on each loop
            datadiag = self.A.diagonal()
            datadiag[Ps] = datadiag[Ps] + phase[item+'_S1'][Ps]
            self.A.setdiag(datadiag)
            # Add S2 to b
            self.b[Ps] = self.b[Ps] - phase[item+'_S2'][Ps]

    def run(self, x=None):
        r"""
        Builds the A and b matrices, and calls the solver specified in the
        ``settings`` attribute.

        Parameters
        ----------
        x : ND-array
            Initial guess of unknown variable

        """
        print('―'*80)
        print('Running ReactiveTransport')
        self.setup()
        self._run_reactive(x=x)

    def _run_reactive(self, x):
        if self.settings['quantity'] not in self.keys():
            self[self.settings['quantity']] = 0
        self.setup()
        self.apply_sources()
        if x is None:
            x = np.zeros(shape=[self.Np, ], dtype=float)
        x_new = self._solve()
        self[self.settings['quantity']] = x_new
        res = np.sum(np.absolute(x**2 - x_new**2))
        if res < self.settings['tolerance']:
            print('Solution converged: ' + str(res))
            return
        else:
            print('Tolerance not met: ' + str(res))
            self._run_reactive(x=x_new)
