"""
===============================================================================
module Dispersion: Class for solving dispersion
===============================================================================
"""

import scipy as sp
from openpnm.algorithms import GenericLinearTransport
from openpnm.core import logging
logger = logging.getLogger(__name__)


class Dispersion(GenericLinearTransport):
    r'''
    A subclass of GenericLinearTransport to simulate dispersion.

    Examples
    --------
    >>>     
    
    '''
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        logger.info('Create ' + self.__class__.__name__ + ' Object')
        
    def setup(self, conductance='diffusive_conductance',
              quantity='mole_fraction', **kwargs):
        r"""

        """
        super().setup(conductance=conductance, quantity=quantity)

    def build_A(self):
        r"""
        """
        network = self.simulation.network
        phase = self.simulation.phases[self.settings['phase']]
        
        # Get adjacancy and incidence matrices
        adj_mat = network.create_adjacency_matrix(fmt='csr')
        inc_mat = network.create_incidence_matrix(fmt='csr')
        
        K = network['throat.conns']
        P = phase['pore.pressure'] # Phase must have 'pressure' attached to it.
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
            Q_ik = (Pi - Pk) * phase['throat.conductance'][neighbor_throats]
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
        
        if 'pore.neumann' in self.keys():
            pass  # Do nothing to A, only b changes
        if 'pore.dirichlet' in self.keys():
            # Find all entries on rows associated with dirichlet pores
            P_bc = self.toindices(self['pore.dirichlet'])
            indrow = sp.in1d(A.row, P_bc)
            A.data[indrow] = 0  # Remove entries from A for all BC rows
            datadiag = A.diagonal()  # Add diagonal entries back into A
            datadiag[P_bc] = sp.ones_like(P_bc, dtype=float)
            A.setdiag(datadiag)
            A.eliminate_zeros()  # Remove 0 entries
        self.A = A
        return A
