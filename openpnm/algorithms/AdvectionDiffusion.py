import numpy as np
import scipy.sparse as sprs
from scipy.sparse.csgraph import laplacian
from openpnm.algorithms import GenericTransport
from openpnm.core import logging
logger = logging.getLogger(__name__)


class AdvectionDiffusion(GenericTransport):
    r"""
    A subclass of GenericTransport to simulate advection diffusion.

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
        return self.build_A_vectorized()

    def build_A_vectorized(self):
        network = self.project.network
        phase = self.project.phases()[self.settings['phase']]
        P = phase[self.settings['pressure']]
        gh = phase[self.settings['hydraulic_conductance']]
        gh = np.tile(gh, 2)
        gd = phase[self.settings['diffusive_conductance']]
        gd = np.tile(gd, 2)

        pores_ij = network['throat.conns']
        conns = np.vstack((pores_ij, np.flip(pores_ij, axis=1)))
        Qij = gh*np.diff(P[conns], axis=1).squeeze()
        qP = np.where(Qij > 0, Qij, 0)  # Throat positive flow rates
        qN = np.where(Qij < 0, Qij, 0)
        # Put the flow rates in the coefficient matrix
        A = network.create_adjacency_matrix(weights=(qP + gd))
        # Overwrite the diagonal
        am = network.create_adjacency_matrix(weights=(qN - gd))
        A_diags = laplacian(am)
        A.setdiag(A_diags.diagonal())
        self.A = A
        return A

    def build_A_looped(self):
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
            q = gh[nt[i]]*(P[i]-P[network.find_neighbor_pores(i)])  # Flow rate
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
