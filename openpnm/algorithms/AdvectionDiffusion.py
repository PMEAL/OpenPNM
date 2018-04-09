import numpy as np
import scipy.sparse as sprs
from openpnm.algorithms import GenericTransport
from openpnm.core import logging
logger = logging.getLogger(__name__)


class AdvectionDiffusion(GenericTransport):
    r"""
    A subclass of GenericAlgorithm to simulate advection diffusion.

    """

    def __init__(self, settings={}, **kwargs):
        super().__init__(**kwargs)
        # Set some default settings
        self.settings.update({'quantity': 'pore.mole_fraction',
                              'diffusive_conductance':
                                  'throat.diffusive_conductance',
                              'hydraulic_conductance':
                                  'throat.hydraulic_conductance',
                              'pressure': 'pore.pressure'})
        # Apply any received settings to overwrite defaults
        self.settings.update(settings)

    def build_A(self):
        r"""
        """
        network = self.project.network
        phase = self.project.phases()[self.settings['phase']]
        gh = phase[self.settings['hydraulic_conductance']]
        P = phase[self.settings['pressure']]
        D = phase[self.settings['diffusive_conductance']]
        # Pore neighbor throats
        nt = network.find_neighbor_throats(pores=phase['pore.all'],
                                           flatten=False,
                                           mode='not_intersection')
        A = np.zeros((network.Np, network.Np))  # Initialize A matrix
        for i in range(network.Np):
            q = gh[nt[i]]*(P[network.find_neighbor_pores(i)]-P[i])  # Flow rate
            qP = np.where(q > 0, q, 0)  # Throat positive flow rates
            qN = np.where(q < 0, q, 0)  # Throat negative flow rates
            A[i, i] = np.sum(qN - D[nt[i]])  # Diagonal
            j1 = network['throat.conns'][nt[i]]  # Find off diag cells to fill
            j2 = np.reshape(j1, np.size(j1))
            j = j2[j2 != i]
            A[i, j] = -(-qP - D[nt[i]])  # Off diagonal
        A = sprs.coo_matrix(A)
        self.A = A
        return A
