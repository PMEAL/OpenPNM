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
        logger.info('Create ' + self.__class__.__name__ + ' Object')
        self.settings.update({'quantity': 'pore.mole_fraction',
                              'hydraulic_conductance':
                              'throat.hydraulic_conductance',
                              'diffusivity': 'pore.diffusivity',
                              'molar_density': 'pore.molar_density',
                              'pressure': 'pore.pressure'})
        self.settings.update(settings)

    def _build_A(self, force=False):
        r"""
        """
        network = self.project.network
        phase = self.project.phases()[self.settings['phase']]
        P = phase[self.settings['pressure']]
        D = np.mean(phase[self.settings['diffusivity']])
        gh_0 = phase[self.settings['hydraulic_conductance']]
        conns1 = network['throat.conns']
        conns2 = np.flip(conns1, axis=1)

        # Calculating effective length and area
        dp = network['pore.diameter']
        Vp = network['pore.volume']
        Vt = network['throat.volume']
        Lt = network['throat.length']
        Le_0 = Lt + (dp[conns1[:, 0]] + dp[conns1[:, 1]])/2
        Le = np.tile(Le_0, 2)
        Ae = (Vt + Vp[conns1[:, 0]] + Vp[conns1[:, 1]]) / Le_0
        Ae = np.tile(Ae, 2)

        Qij1 = gh_0*np.diff(P[conns1], axis=1).squeeze()
        Qij1 = np.append(Qij1, -Qij1)
        Uij1 = Qij1 / Ae
        Peij1 = Uij1 * Le / D

        Qij2 = gh_0*np.diff(P[conns2], axis=1).squeeze()
        Qij2 = np.append(Qij2, -Qij2)
        Uij2 = Qij2 / Ae
        Peij2 = Uij2 * Le / D

        Peij1[(Peij1 < 1e-10) & (Peij1 >= 0)] = 1e-10
        Peij1[(Peij1 > -1e-10) & (Peij1 <= 0)] = -1e-10
        Peij1[Peij1 > 100] = 100
        Peij1[Peij1 < -100] = -100

        Peij2[(Peij2 < 1e-10) & (Peij2 >= 0)] = 1e-10
        Peij2[(Peij2 > -1e-10) & (Peij2 <= 0)] = -1e-10
        Peij2[Peij2 > 100] = 100
        Peij2[Peij2 < -100] = -100

        if force:
            self._pure_A = None
        if self._pure_A is None:
            w = -Qij1 + Qij1 / (1 - np.exp(Peij1))
            am1 = -network.create_adjacency_matrix(weights=w)
            w = -Qij2 / (1 - np.exp(Peij2))
            A = -network.create_adjacency_matrix(weights=w)
            A_diags = laplacian(am1)
            # Overwrite the diagonal
            A.setdiag(A_diags.diagonal())
            self._pure_A = A
        self.A = self._pure_A.copy()
