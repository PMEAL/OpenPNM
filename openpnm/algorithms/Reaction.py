import scipy as sp
from scipy.optimize import newton_krylov, anderson, broyden1, broyden2
from openpnm.core import ModelsMixin
from openpnm.utils.misc import PrintableDict
from openpnm.algorithms import GenericAlgorithm


class GenericReaction(GenericAlgorithm, ModelsMixin):

    def __init__(self, network, algorithm, pores, **kwargs):
        phase = algorithm.simulation.phases[algorithm.settings['phase']]
        super().__init__(network=network, phase=phase, **kwargs)
        self.update({'pore.all': sp.ones_like(pores, dtype=bool)})
        self.update({'throat.all': sp.ones_like(phase.Ts, dtype=bool)})
        self.settings = PrintableDict({'phase': phase.name,
                                       'algorithm': algorithm.name,
                                       'quantity': None})
        self['pore._id'] = network['pore._id'][pores]
        # Add self to algorithm's settings
        algorithm.settings['sources'].append(self.name)

    def setup(self, quantity=None):
        if quantity is not None:
            self.settings['quantity'] = quantity
        # Fetch algorithm object from simulation
        alg = self.simulation[self.settings['algorithm']]
        # Fetch a copy of b from algorithm object
        self.b = alg.b.copy()

    def solve(self, x0=None):
        if x0 is None:
            x0 = 0.5*sp.ones(shape=(self.Np, ))
#        x = newton_krylov(F=self._residual, xin=x0, method='gmres', verbose=1)
        x = anderson(F=self._residual, xin=x0, verbose=1)
#        x = broyden1(F=self._residual, xin=x0, verbose=1)
#        x = broyden2(F=self._residual, xin=x0, verbose=1)
        return x

    def _residual(self, x):
        net = self.simulation.network
        Ps = net.map_pores(self['pore._id'])
        # Fetch algorithm object from simulation
        alg = self.simulation[self.settings['algorithm']]
        self[self.settings['quantity']] = x
        # Regenerate models with new guess
        self.regenerate_models()
        # Adjust b vector with rate based on new guess of x
        alg.b[Ps] = self.b[Ps] + self['pore.rate']
        x_new = alg.solve()[Ps]
        res = x_new - x
        return res
