import scipy as sp
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
        # Fetch a copy of A and b from algorithm object
        self.b = alg.b.copy()
        self.A = alg.A.copy()

    def run(self, x=None, tol=0.001):
        if x is None:
            x = sp.zeros_like(self.Ps)
        net = self.simulation.network
        Ps = net.map_pores(self['pore._id'])
        # Fetch algorithm object from simulation
        alg = self.simulation[self.settings['algorithm']]
        self[self.settings['quantity']] = x
        # Regenerate models with new guess
        self.regenerate_models()
        # Add S1 to diagonal of A
        datadiag = self.A.diagonal()
        datadiag[Ps] = datadiag[Ps] + self['pore.rate'][:, 1]
        alg.A.setdiag(datadiag)
        # Add S2 to b
        alg.b[Ps] = self.b[Ps] - self['pore.rate'][:, 2]
        # Call the solver on the algorithm
        x_new = alg.solve()[Ps]
        # Check the residual difference between current x and previous value
        res = sp.sum(sp.absolute(x_new**2 - x**2))
        # Recursively call self.run until desired tolerance is met
        if res < tol:
            print('Tolerance met, current residual: '+str(res))
            return x
        else:
            print('Tolerance not met, current residual: '+str(res))
            self.run(x_new)
