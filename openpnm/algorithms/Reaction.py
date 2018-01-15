import scipy as sp
from scipy.optimize import anderson
from openpnm.core import ModelsMixin
from openpnm.utils.misc import PrintableDict
from openpnm.algorithms import GenericAlgorithm


class GenericReaction(GenericAlgorithm, ModelsMixin):

    def __init__(self, network, pores, **kwargs):
        super().__init__(network=network, **kwargs)
        self.settings = PrintableDict()
        self.update({'pore.all': sp.ones_like(pores, dtype=bool)})
        self.update({'throat.all': sp.ones(shape=(0, ), dtype=bool)})
        self['pore._id'] = network['pore._id'][pores]
        self.settings['rate_model'] = 'pore.rate'

    def setup(self, algorithm):
        self.settings['algorithm'] = algorithm.name
        self.settings['quantity'] = algorithm.settings['quantity']

    def solve(self, x0=None):
        if x0 is None:
            x0 = 0.5*sp.ones(shape=(self.Np, ))
        x = anderson(F=self._residual, xin=x0, verbose=1)
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

    def apply(self):
        net = self.simulation.network
        Ps = net.map_pores(self['pore._id'])
        # Fetch algorithm object from simulation
        alg = self.simulation[self.settings['algorithm']]
        quantity = alg.settings['quantity']
        x = alg[quantity].copy()
        self[quantity] = x[Ps]
        # Regenerate models with new guess
        self.regenerate_models()
        # Add S1 to diagonal of A
        datadiag = alg.A.diagonal()
        datadiag[Ps] = datadiag[Ps] +self[self.settings['rate_model']][:, 1]
        alg.A.setdiag(datadiag)
        alg.b[Ps] = alg.b[Ps] - self[self.settings['rate_model']][:, 2]
