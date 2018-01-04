import scipy as sp
from scipy.optimize import newton_krylov, anderson, broyden1, broyden2
from openpnm.core import ModelsMixin
from openpnm.utils.misc import PrintableDict
from openpnm.algorithms import GenericAlgorithm


class GenericReaction(GenericAlgorithm, ModelsMixin):

    def __init__(self, network, algorithm, **kwargs):
        phase = algorithm.simulation.phases[algorithm.settings['phase']]
        super().__init__(network=network, phase=phase, **kwargs)
        self.update({'pore.all': sp.ones_like(phase.Ps, dtype=bool)})
        self.update({'throat.all': sp.ones_like(phase.Ts, dtype=bool)})
        self.settings = PrintableDict({'phase': phase.name,
                                       'algorithm': algorithm.name,
                                       'quantity': None})
        self['pore._id'] = network['pore._id']
        # Add self to algorithm's settings
        algorithm.settings['sources'].append(self.name)

    def setup(self, quantity=None):
        if quantity is not None:
            self.settings['quantity'] = quantity
        # Fetch algorithm object from simulation
        alg = self.simulation[self.settings['algorithm']]
        # Assign algorithm's present values to self
        self[self.settings['quantity']] = alg[self.settings['quantity']]
        # Fetch a copy of b from algorithm object
        self.b = alg.b.copy()

    def solve(self, x0=None):
        if x0 is None:
            x0 = self[self.settings['quantity']]
#        x = newton_krylov(F=self._residual, xin=x0, method='lgmres', verbose=1)
#        x = anderson(F=self._residual, xin=x0, verbose=1)
#        x = broyden1(F=self._residual, xin=x0, verbose=1)
        x = broyden2(F=self._residual, xin=x0, verbose=1)
        return x

    def _residual(self, x):
        # Fetch algorithm object from simulation
        alg = self.simulation[self.settings['algorithm']]
        self[self.settings['quantity']] = x
        # Regenerate models with new guess
        self.regenerate_models()
        # Adjust b vector with rate based on new guess of x
        alg.b = self.b + self['pore.rate']
        x_new = alg.solve()
        res = x_new - x
        return res
