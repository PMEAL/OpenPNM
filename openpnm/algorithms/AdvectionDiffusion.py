import numpy as np
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
                              'pressure': 'pore.pressure',
                              'molar_density': 'pore.molar_density'})
        # Apply any received settings to overwrite defaults
        self.settings.update(settings)

    def build_A(self):
        network = self.project.network
        phase = self.project.phases()[self.settings['phase']]
        P = phase[self.settings['pressure']]
        gh_0 = phase[self.settings['hydraulic_conductance']]
        gd = phase[self.settings['diffusive_conductance']]
        gd = np.tile(gd, 2)
        md = phase[self.settings['molar_density']][0]

        conns = np.flip(network['throat.conns'], axis=1)
        Qij = md*gh_0*np.diff(P[conns], axis=1).squeeze()
        Qij = np.append(Qij, -Qij)

        qP = np.where(Qij > 0, Qij, 0)  # Throat positive flow rates
        qN = np.where(Qij < 0, Qij, 0)

        # Put the flow rates in the coefficient matrices
        am1 = network.create_adjacency_matrix(weights=(qP + gd))
        am2 = network.create_adjacency_matrix(weights=(-qN + gd))
        A_diags = laplacian(am1)
        A = laplacian(am2)
        # Overwrite the diagonal
        A.setdiag(A_diags.diagonal())
        self.A = A
        return A
