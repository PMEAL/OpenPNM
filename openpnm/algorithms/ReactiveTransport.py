import numpy as np
from openpnm.algorithms import GenericTransport
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class ReactiveTransport(GenericTransport):
    r"""
    """

    def __init__(self, **kwargs):
        self.settings.update({'sources': [],
                              'tolerance': 0.001,
                              'max_iter': 10000,
                              'relaxation_source': 1,
                              'relaxation_quantity': 1})
        super().__init__(**kwargs)

    def set_source(self, propname, pores):
        r"""
        Applies a given source term to the specified pores

        Parameters
        ----------
        propname : string
            The property name of the source term model to be applied

        pores : array_like
            The pore indices where the source term should be applied

        """
        self.settings['sources'].append(propname)
        self[propname] = self.tomask(pores=pores)

    def _update_physics(self):
        phase = self.project.phases()[self.settings['phase']]
        physics = self.project.find_physics(phase=phase)
        for item in self.settings['sources']:
            # Regenerate models with new guess
            quantity = self.settings['quantity']
            # Put quantity on phase so physics finds it when regenerating
            phase[quantity] = self[quantity]
            # Regenerate models, on either phase or physics
            phase.regenerate_models(propnames=item)
            for phys in physics:
                phys.regenerate_models(propnames=item)

    def _apply_sources(self):
        if self.settings['t_scheme'] == 'cranknicolson':
            f1 = 0.5
        else:
            f1 = 1
        phase = self.project.phases()[self.settings['phase']]
        relax = self.settings['relaxation_source']
        for item in self.settings['sources']:
            Ps = self.pores(item)
            # Add S1 to diagonal of A
            # TODO: We need this to NOT overwrite the A and b, but create
            # copy, otherwise we have to regenerate A and b on each loop
            datadiag = self._A.diagonal().copy()
            # Source term relaxation
            S1_old = phase[item+'.'+'S1'][Ps].copy()
            S2_old = phase[item+'.'+'S2'][Ps].copy()
            self._update_physics()
            S1 = phase[item+'.'+'S1'][Ps]
            S2 = phase[item+'.'+'S2'][Ps]
            S1 = relax*S1 + (1-relax)*S1_old
            S2 = relax*S2 + (1-relax)*S2_old
            phase[item+'.'+'S1'][Ps] = S1
            phase[item+'.'+'S2'][Ps] = S2
            datadiag[Ps] = datadiag[Ps] - f1*S1
            # Add S1 to A
            self._A.setdiag(datadiag)
            # Add S2 to b
            self._b[Ps] = self._b[Ps] + f1*S2

    def run(self, x=None):
        r"""
        Builds the A and b matrices, and calls the solver specified in the
        ``settings`` attribute.

        Parameters
        ----------
        x : ND-array
            Initial guess of unknown variable

        relaxation : float
                Under relaxation parameter, a float between 0 and 1.
                relaxation = 1, no relaxation (fast but may be unstable).
                relaxation < 1, improved stability but may be slow.

        """
        logger.info('―'*80)
        logger.info('Running ReactiveTransport')
        x = self._run_reactive(x=x)
        return x

    def _run_reactive(self, x):
        if x is None:
            x = np.zeros(shape=[self.Np, ], dtype=float)
        self[self.settings['quantity']] = x
        relax = self.settings['relaxation_quantity']
        res = 1e+06
        for itr in range(int(self.settings['max_iter'])):
            if res >= self.settings['tolerance']:
                self[self.settings['quantity']] = x
                self._build_A(force=True)
                self._build_b(force=True)
                self._apply_BCs()
                self._apply_sources()
                x_new = self._solve()
                # Relaxation
                x_new = relax*x_new + (1-relax)*self[self.settings['quantity']]
                self[self.settings['quantity']] = x_new
                res = np.sum(np.absolute(x**2 - x_new**2))
                x = x_new
                logger.info('Tolerance not met: ' + str(res))
            elif res < self.settings['tolerance']:
                logger.info('Solution converged: ' + str(res))
                break
        return x_new
