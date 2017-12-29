import scipy as sp
from openpnm.core import ModelsMixin
from openpnm.utils.misc import PrintableDict
from openpnm.algorithms import GenericAlgorithm


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

    def setup(self, quantity):
        # Initialize an array to store the current value of the quantity
        self['pore.quantity'] = 0
        self.settings['quantity'] = quantity

    def apply(self):
        # Fetch algorithm object from simulaiton
        alg = self.simulation[self.settings['algorithm']]
        # Identify which pores this object applies to
        Ps = self.map_pores()
        # Retrieve info from settings and store on self
        self[self.settings['quantity']] = alg[self.settings['quantity']][Ps]
        # Run models to generate slope and intercept of Taylor series
        self.regenerate_models()
        # Fetch A and b from algorithm object
        A = alg.A
        b = alg.b
        # Determine which elements in the COO of A are the diagonal
        diag = sp.where(A.col == A.row)[0]
        # Add source term to the diagonal of the coefficient matrix
        A.data[diag][Ps] -= self['pore.S1']
        # Add source term to the RHS matrix
        b[Ps] += self['pore.S2']
