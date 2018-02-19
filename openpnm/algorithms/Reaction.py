import scipy as sp
from openpnm.core import ModelsMixin
from openpnm.algorithms import GenericAlgorithm


class GenericReaction(GenericAlgorithm, ModelsMixin):

    def __init__(self, network, pores, settings={}, **kwargs):
        self.settings.update({'rate_model': 'pore.rate',
                              'prefix': 'rxn'})
        self.settings.update(settings)
        super().__init__(Np=sp.size(pores), network=network, **kwargs)
        self['pore._id'] = network['pore._id'][pores]

    def setup(self, algorithm):
        self.settings['algorithm'] = algorithm.name
        self.settings['quantity'] = algorithm.settings['quantity']

    def apply(self, A=None, b=None):
        net = self.project.network
        Ps = net.map_pores(self['pore._id'])
        # Fetch algorithm object from project
        alg = self.project[self.settings['algorithm']]
        if A is None:
            A = alg.A
        if b is None:
            b = alg.b
        quantity = alg.settings['quantity']
        x = alg[quantity].copy()
        self[quantity] = x[Ps]
        # Regenerate models with new guess
        self.regenerate_models()
        # Add S1 to diagonal of A
        datadiag = A.diagonal()
        datadiag[Ps] = datadiag[Ps] + self[self.settings['rate_model']][:, 1]
        A.setdiag(datadiag)
        b[Ps] = b[Ps] - self[self.settings['rate_model']][:, 2]
