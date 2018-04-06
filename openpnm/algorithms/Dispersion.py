"""
===============================================================================
module Dispersion: Class for solving advection diffusion
===============================================================================
"""

import scipy as sp
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
        
        # Get adjacancy and incidence matrices
        adj_mat = network.create_adjacency_matrix(fmt='csr')
        inc_mat = network.create_incidence_matrix(fmt='csr')
        
        K = network['throat.conns']
        # Phase must have 'pressure' attached to it.
        P = phase[self.settings['pressure']]
        dp = network['pore.diameter']
        Vp = network['pore.volume']
        Vt = network['throat.volume']
        Lt = network['throat.length']
        D = sp.mean(phase['pore.diffusivity'])
        
        # Calculating effective length and area
        Le = Lt + (dp[K[:,0]] + dp[K[:,1]])/2
        Ae = (Vt + Vp[K[:,0]] + Vp[K[:,1]]) / Le
        
        I, J, V = [], [], []
        for i in range(network.Np):   
            neighbor_pores = adj_mat[i,:].nonzero()[1]
            neighbor_throats = inc_mat[i,:].nonzero()[1]
            
            neighbor_conns = K[neighbor_throats]
            neighbor_pores_w_throat_ordering = neighbor_conns[neighbor_conns!=i]
            Pi, Pk = P[i], P[neighbor_pores_w_throat_ordering]
            Q_ik = (Pi - Pk) * phase[self.settings['conductance']][neighbor_throats]
            Le_ik = Le[neighbor_throats]
            Ae_ik = Ae[neighbor_throats]
            u_ik = Q_ik / Ae_ik
            Pe_ik = u_ik * Le_ik / D
     
            # Get rid of exp overflow when calculating condunctances
            Pe_ik[(Pe_ik < 1e-10) & (Pe_ik >= 0)] = 1e-10
            Pe_ik[(Pe_ik > -1e-10) & (Pe_ik <= 0)] = -1e-10
            Pe_ik[Pe_ik > 100] = 100
            Pe_ik[Pe_ik < -100] = -100
                
            # Coefficients of center pore
            val = sp.sum(-Q_ik + Q_ik / (1 - sp.exp(Pe_ik)))
            V.append(val)
            I.append(i)
            J.append(i)
            
            # Coefficients of neighbour pores
            val = -Q_ik / (1 - sp.exp(Pe_ik))
            V.extend(val)
            I.extend([i]*len(neighbor_pores))
            J.extend(neighbor_pores_w_throat_ordering)
        A = sp.sparse.coo_matrix((V,(I,J)),shape=(network.Np, network.Np))
        self.A = A
        return A
