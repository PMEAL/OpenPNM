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
        self._update_physics()
        for item in self.settings['sources']:
            Ps = self.pores(item)
            # Add S1 to diagonal of A
            # TODO: We need this to NOT overwrite the A and b, but create
            # copy, otherwise we have to regenerate A and b on each loop
            datadiag = self._A.diagonal().copy()
            datadiag[Ps] = datadiag[Ps] - f1*phase[item+'.'+'S1'][Ps]
            self._A.setdiag(datadiag)
            # Add S2 to b
            self._b[Ps] = self._b[Ps] + f1*phase[item+'.'+'S2'][Ps]

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
        x = self._run_reactive(x=x)
        return x

    def _run_reactive(self, x):
        if x is None:
            x = np.zeros(shape=[self.Np, ], dtype=float)
        self[self.settings['quantity']] = x
        self._build_A()
        self._build_b()
        self._apply_BCs()
        self._apply_sources()
        x_new = self._solve()
        self[self.settings['quantity']] = x_new
        res = np.sum(np.absolute(x**2 - x_new**2))
        if res < self.settings['tolerance']:
            logger.info('Solution converged: ' + str(res))
            return x_new
        else:
            logger.info('Tolerance not met: ' + str(res))
            self._run_reactive(x=x_new)
