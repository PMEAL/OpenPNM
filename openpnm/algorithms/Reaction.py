import scipy as sp
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
                                       'quantity': algorithm.settings['quantity']})
        self['pore._id'] = network['pore._id']
        # Add self to algorithm's settings
        algorithm.settings['sources'].append(self.name)

    def setup(self, quantity=None):
        if quantity is not None:
            self.settings['quantity'] = quantity
        # Fetch algorithm object from simulation
        alg = self.simulation[self.settings['algorithm']]
        # Fetch a copy of A and b from algorithm object
        self.b = alg.b.copy()
        self.A = alg.A.copy()

    def solve(self, x0=None):
        if x0 is None:
            x0 = 0.5*sp.ones(shape=(self.Np, ))
#        x = newton_krylov(F=self._residual, xin=x0, method='gmres', verbose=1)
        x = anderson(F=self._residual, xin=x0, verbose=1)
#        x = broyden1(F=self._residual, xin=x0, verbose=1)
#        x = broyden2(F=self._residual, xin=x0, verbose=1)
        return x

    def _residual(self, x):
        # Fetch algorithm object from simulation
        alg = self.simulation[self.settings['algorithm']]
        self[self.settings['quantity']] = x
        # Regenerate models with new guess
        self.regenerate_models()
        # Adjust b vector with rate based on new guess of x
        alg.b = self.b + self['pore.rate'][:, 0]
        x_new = alg.solve()
        res = x_new - x
        return res

    def run(self, x=None, tol=0.001):
        if x is None:
            x = sp.zeros_like(self.Ps)
        # Fetch algorithm object from simulation
        alg = self.simulation[self.settings['algorithm']]
        self[self.settings['quantity']] = x
        # Regenerate models with new guess
        self.regenerate_models()
        # Add S1 to diagonal of A
        datadiag = self.A.diagonal()
        datadiag = datadiag + self['pore.rate'][:, 1]
        alg.A.setdiag(datadiag)
        alg.b = self.b - self['pore.rate'][:, 2]
        x_new = alg.solve()
        res = sp.sum(sp.absolute(x_new - x))
        if res < tol:
            print('Tolerance met, current residual: '+str(res))
            return x
        else:
            print('Tolerance not met, current residual: '+str(res))
            self.run(x_new)
