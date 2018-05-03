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
                              'molar_density': 'pore.molar_density',
                              's_scheme': 'powerlaw'})
        # Apply any received settings to overwrite defaults
        self.settings.update(settings)

    def setup(self, phase=None, quantity='', diffusive_conductance='',
              hydraulic_conductance='', pressure='', molar_density='',
              s_scheme='', **kwargs):
        r"""

        """
        if phase:
            self.settings['phase'] = phase.name
        if quantity:
            self.settings['quantity'] = quantity
        if diffusive_conductance:
            self.settings['diffusive_conductance'] = diffusive_conductance
        if hydraulic_conductance:
            self.settings['hydraulic_conductance'] = hydraulic_conductance
        if pressure:
            self.settings['pressure'] = pressure
        if molar_density:
            self.settings['molar_density'] = molar_density
        if s_scheme:
            self.settings['s_scheme'] = s_scheme
        super().setup(**kwargs)

    def _build_A(self, force=False):
        s_dis = self.settings['s_scheme']
        network = self.project.network
        phase = self.project.phases()[self.settings['phase']]
        P = phase[self.settings['pressure']]
        gh_0 = phase[self.settings['hydraulic_conductance']]
        gd = phase[self.settings['diffusive_conductance']]
        gd = np.tile(gd, 2)
        md = phase[self.settings['molar_density']][0]

        conns1 = network['throat.conns']
        conns2 = np.flip(conns1, axis=1)

        Qij1 = md*gh_0*np.diff(P[conns1], axis=1).squeeze()
        Qij1 = np.append(Qij1, -Qij1)
        Qij2 = md*gh_0*np.diff(P[conns2], axis=1).squeeze()
        Qij2 = np.append(Qij2, -Qij2)

        qP1 = np.where(Qij1 > 0, Qij1, 0)  # Throat positive flow rates
        qN2 = np.where(Qij2 < 0, Qij2, 0)

        if force:
            self._pure_A = None
        if self._pure_A is None:
            if (s_dis == 'upwind'):
                # Put the flow rates in the coefficient matrices
                am1 = network.create_adjacency_matrix(weights=(qP1 + gd))
                A = -network.create_adjacency_matrix(weights=(-qN2 + gd))
            elif (s_dis == 'hybrid'):
                w = np.maximum(qP1, (qP1/2 + gd))
                am1 = network.create_adjacency_matrix(weights=w)
                w = np.maximum(-qN2, (-qN2/2 + gd))
                A = -network.create_adjacency_matrix(weights=w)
            elif (s_dis == 'powerlaw'):
                Pe1 = np.absolute(qP1/gd)
                Pe2 = np.absolute(qN2/gd)
                diag = gd*np.maximum(0, (1-0.1*Pe1)**5) + qP1
                offdiag = gd*np.maximum(0, (1-0.1*Pe2)**5) - qN2
                am1 = network.create_adjacency_matrix(weights=diag)
                A = -network.create_adjacency_matrix(weights=offdiag)
            A_diags = laplacian(am1)
            # Overwrite the diagonal
            A.setdiag(A_diags.diagonal())
            self._pure_A = A
        self._A = self._pure_A.copy()
        return A
