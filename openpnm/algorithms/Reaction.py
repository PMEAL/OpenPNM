import scipy as sp
from openpnm.core import ModelsMixin
from openpnm.utils.misc import PrintableDict
from openpnm.algorithms import GenericAlgorithm
import scipy.sparse as sprs


class GenericReaction(GenericAlgorithm, ModelsMixin):

    def __init__(self, network, algorithm, pores, **kwargs):
        phase = algorithm.simulation.phases[algorithm.settings['phase']]
        super().__init__(network=network, phase=phase, **kwargs)
        pores = self._parse_indices(pores)
        self.update({'pore.all': sp.ones_like(pores, dtype=bool)})
        self.update({'throat.all': sp.ones(shape=(0, ), dtype=bool)})
        self.settings = PrintableDict({'phase': phase.name,
                                       'algorithm': algorithm.name,
                                       'quantity': None})
        self['pore._id'] = network['pore._id'][pores]

    def setup(self, quantity):
        self.settings['quantity'] = quantity

    def apply(self):
        # Fetch algorithm object from simulaiton
        alg = self.simulation[self.settings['algorithm']]
        # Identify which pores this object applies to
        net = self.simulation.network
        Ps = net.map_pores(self['pore._id'])
        # Retrieve info from settings and store on self
        quant = self.settings['quantity']
        if quant not in alg.keys():
            alg[quant] = .5
        self[quant] = alg[quant][Ps]
        # Run models to generate slope and intercept of Taylor series
        self.regenerate_models()
        # Fetch A and b from algorithm object
        A = alg.A
        b = alg.b
        # Determine which elements in the COO of A are the diagonal
        diag = sp.where(A.col == A.row)[0]
        # Add source term to the diagonal of the coefficient matrix
        A.data[diag[Ps]] -= self['pore.source_term'][:, 0]
        # Add source term to the RHS matrix
        b[Ps] += self['pore.source_term'][:, 1]
        alg.run = self.run

    def run(self, rel_tol=1e-5, max_iter=5):
        # Fetch algorithm object from simulaiton
        alg = self.simulation[self.settings['algorithm']]
        count = 0
        while count < max_iter:
            count += 1
            x_old = alg[self.settings['quantity']].copy()
            x = sprs.linalg.spsolve(A=alg.A.tocsr(), b=alg.b)
            print(sp.vstack((x_old, x)).T)
            if sp.allclose(x, x_old, rtol=rel_tol, atol=0):
                break
            alg[self.settings['quantity']] = x
            alg.build_A()
            alg.build_b()
            self.apply()
            print(count)
        return x
