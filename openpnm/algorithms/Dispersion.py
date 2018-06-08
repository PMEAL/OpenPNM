"""
===============================================================================
module Dispersion: Class for solving advection diffusion
===============================================================================
"""

import numpy as np
from scipy.sparse.csgraph import laplacian
from openpnm.algorithms import ReactiveTransport
from openpnm.core import logging
logger = logging.getLogger(__name__)


class Dispersion(ReactiveTransport):
    r'''
    A subclass of GenericTransport to simulate advection diffusion.

    Examples
    --------
    >>>

    '''
    def __init__(self, settings={}, **kwargs):
        super().__init__(**kwargs)
        self.settings.update({'quantity': 'pore.concentration',
                              'hydraulic_conductance':
                              'throat.hydraulic_conductance',
                              'diffusive_conductance':
                              'throat.diffusive_conductance',
                              'pressure': 'pore.pressure'})
        self.settings.update(settings)

    def _build_A(self, force=False):
        r"""
        """
        network = self.project.network
        phase = self.project.phases()[self.settings['phase']]
        conns = network['throat.conns']

        P = phase[self.settings['pressure']]
        gh = phase[self.settings['hydraulic_conductance']]
        gd = phase[self.settings['diffusive_conductance']]
        gd = np.tile(gd, 2)

        Qij = -gh*np.diff(P[conns], axis=1).squeeze()
        Qij = np.append(Qij, -Qij)
        Peij = Qij/gd

        Peij[(Peij < 1e-10) & (Peij >= 0)] = 1e-10
        Peij[(Peij > -1e-10) & (Peij <= 0)] = -1e-10
        Qij = Peij*gd

        if force:
            self._pure_A = None
        if self._pure_A is None:
            w = -Qij / (1 - np.exp(Peij))
            A = network.create_adjacency_matrix(weights=w)
            A = laplacian(A)
            self._pure_A = A
        self.A = self._pure_A.copy()
