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

    def build_A(self):
        r"""
        """
        network = self.project.network
        phase = self.project.phases()[self.settings['phase']]

        # Get conns for triu and full matrix
        pore_ij = network['throat.conns']
        pore_ij = sp.flip(pore_ij, axis=1)
        conns = sp.vstack((pore_ij, sp.flip(pore_ij, axis=1)))

        # Fetch phase properties, including pressure
        # TODO: Adding rands to prevent error at delta P = 0...could be better?
        P = phase['pore.pressure'] + sp.rand(self.Np)*1e-30
        # TODO: convert to a throat vector to account for spatial variation
        D = sp.mean(phase['pore.diffusivity'])

        # Fetch geometric properties for pores and throats
        Dp = network['pore.diameter']
        Vp = network['pore.volume']
        Vt = network['throat.volume']
        Lt = network['throat.length']

        # Calculate effective length and area specifically for dispersion
        # TODO: Move this to a pore scale model and make more rigorous
        Le = Lt + sp.mean(Dp[pore_ij], axis=1)
        Ae = (Vt + sp.sum(Vp[pore_ij], axis=1))/Le
        Le_ik = sp.tile(Le, 2)
        Ae_ik = sp.tile(Ae, 2)

        # Calculate flow conditions in each throat, for both directions
        g = phase[self.settings['conductance']]
        g = sp.tile(g, 2)
        Q_ik = g*sp.diff(P[conns], axis=1).squeeze()
        u_ik = Q_ik/Ae_ik
        Pe_ik = u_ik*Le_ik/D

        # Condition numerical extremes in Pe_ik array
        negs = Pe_ik < 0  # Note locations of negative numbers
        temp = sp.absolute(Pe_ik)
        Pe_ik = sp.clip(temp, 1e-10, 100)  # Clip large and near-zeros values
        Pe_ik[negs] = -1.0*Pe_ik[negs]  # Replace negative numbers

        # Calculate the-off diagonal terms
        off_diags = Q_ik/(1 - sp.exp(Pe_ik))
        # Build an adjacency matrix and pass to scipy's laplacian function
        am = network.create_adjacency_matrix(weights=off_diags, fmt='coo')
        A = laplacian(am)

        # Now add -Q_ik to each element of the diagonal
        diag = A.diagonal()
        sp.add.at(diag, conns[:, 0], -Q_ik)
        A.setdiag(diag)

        self.A = A
        return A
