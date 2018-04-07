"""
===============================================================================
module Dispersion: Class for solving advection diffusion
===============================================================================
"""

import scipy as sp
from scipy.sparse.csgraph import laplacian
from openpnm.algorithms import GenericTransport
from openpnm.core import logging
logger = logging.getLogger(__name__)


class Dispersion(GenericTransport):
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
                              'conductance': 'throat.hydraulic_conductance',
                              'pressure': 'pore.pressure'})
        self.settings.update(settings)

    def setup(self, **kwargs):
        r"""

        """
        super().setup()

    def build_A(self):
        r"""
        """
        network = self.project.network
        phase = self.project.phases()[self.settings['phase']]

        # Get conns for triu and full matrix
        pore_ij = network['throat.conns']
        conns = sp.vstack((pore_ij, sp.flip(pore_ij, axis=1)))

        # Fetch phase properties, including pressure
        P = phase['pore.pressure']
        D = sp.mean(phase['pore.diffusivity'])

        # Fetch geometric properties for pores and throats
        Dp = network['pore.diameter']
        Vp = network['pore.volume']
        Vt = network['throat.volume']
        Lt = network['throat.length']

        # Calculating effective length and area specifically for dispersion
        Le = Lt + sp.mean(Dp[pore_ij], axis=1)
        Ae = (Vt + sp.sum(Vp[pore_ij], axis=1)) / Le

        # Calculate conditions in each throat, for both directions
        g = phase['throat.conductance']
        g = sp.tile(g, 2)
        Le_ik = sp.tile(Le, 2)
        Ae_ik = sp.tile(Ae, 2)
        Q_ik = g*sp.diff(P[conns], axis=1).squeeze()
        u_ik = Q_ik / Ae_ik
        Pe_ik = u_ik * Le_ik / D

        # Condition numerical extremes in Pe_ik array
        negs = Pe_ik < 0  # Note locations of negative numbers
        temp = sp.absolute(Pe_ik)
        Pe_ik = sp.clip(temp, 1e-5, 10)  # Clip large and near-zeros
        Pe_ik[negs] = -1.0*Pe_ik[negs]  # Replace negative numbers

        # Finally use Pe_ik to calculate effective conductance in each throat
        vals = -Q_ik + Q_ik / (1 - sp.exp(Pe_ik))

        # Build an adjacency matrix and pass to scipy's laplacian function
        am = network.create_adjacency_matrix(weights=vals, fmt='coo')
        A = laplacian(am)
        self.A = A
        return A
