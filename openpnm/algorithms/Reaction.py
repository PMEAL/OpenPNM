import scipy as sp
from openpnm.core import ModelsMixin
from openpnm.algorithms import GenericAlgorithm


class GenericReaction(GenericAlgorithm, ModelsMixin):

    _prefix = 'rxn'

    def __init__(self, network, pores, **kwargs):
        super().__init__(network=network, **kwargs)
        self.settings.update({'rate_model': 'pore.rate'})
        self.update({'pore.all': sp.ones_like(pores, dtype=bool)})
        self.update({'throat.all': sp.ones(shape=(0, ), dtype=bool)})
        self['pore._id'] = network['pore._id'][pores]

    def setup(self, algorithm):
        self.settings['algorithm'] = algorithm.name
        self.settings['quantity'] = algorithm.settings['quantity']

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
        datadiag[Ps] = datadiag[Ps] + self[self.settings['rate_model']][:, 1]
        alg.A.setdiag(datadiag)
        alg.b[Ps] = alg.b[Ps] - self[self.settings['rate_model']][:, 2]
