import numpy as np
from scipy.sparse.csgraph import laplacian
from openpnm.algorithms import ReactiveTransport
from openpnm.core import logging
logger = logging.getLogger(__name__)


class AdvectionDiffusion(ReactiveTransport):
    r"""
    A subclass of GenericTransport to simulate advection diffusion.

    """

    def __init__(self, settings={}, **kwargs):
        super().__init__(**kwargs)
        # Set some default settings
        self.settings.update({'quantity': 'pore.concentration',
                              'diffusive_conductance':
                              'throat.diffusive_conductance',
                              'hydraulic_conductance':
                              'throat.hydraulic_conductance',
                              'pressure': 'pore.pressure',
                              's_scheme': 'powerlaw'})
        # Apply any received settings to overwrite defaults
        self.settings.update(settings)

    def setup(self, phase=None, quantity='', diffusive_conductance='',
              hydraulic_conductance='', pressure='', s_scheme='', **kwargs):
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
        if s_scheme:
            self.settings['s_scheme'] = s_scheme
        super().setup(**kwargs)

    def _build_A(self, force=False):
        s_dis = self.settings['s_scheme']
        network = self.project.network
        phase = self.project.phases()[self.settings['phase']]
        conns = network['throat.conns']

        P = phase[self.settings['pressure']]
        gh = phase[self.settings['hydraulic_conductance']]
        gd = phase[self.settings['diffusive_conductance']]
        gd = np.tile(gd, 2)

        Qij = -gh*np.diff(P[conns], axis=1).squeeze()
        Qij = np.append(Qij, -Qij)

        if force:
            self._pure_A = None
        if self._pure_A is None:
            if (s_dis == 'upwind'):
                w = gd + np.maximum(0, -Qij)
                A = network.create_adjacency_matrix(weights=w)
            elif (s_dis == 'hybrid'):
                w = np.maximum(0, np.maximum(-Qij, gd-Qij/2))
                A = network.create_adjacency_matrix(weights=w)
            elif (s_dis == 'powerlaw'):
                Peij = np.absolute(Qij/gd)
                w = gd*np.maximum(0, (1-0.1*Peij)**5) + np.maximum(0, -Qij)
                A = network.create_adjacency_matrix(weights=w)
            A = laplacian(A)
            self._pure_A = A
        self.A = self._pure_A.copy()
