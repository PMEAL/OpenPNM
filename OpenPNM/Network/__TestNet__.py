"""
===============================================================================
TestNet: Generate simple cubic network for testing purposes
===============================================================================

"""

import scipy as sp
from OpenPNM.Network import GenericNetwork
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)


class TestNet(GenericNetwork):
    r"""
    A small cubic network for quick testing purposes

    Parameters
    ----------
    This class accepts no arguments

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.generate()

    def generate(self):
        r"""
        Create test network, of cubic geometry [5, 5, 5]

        Parameters
        ----------
        This network type accepts no arguments
        """
        self._generate_setup()
        self._generate_pores()
        self._generate_throats()
        self._add_labels()
        return self

    def _generate_setup(self):
        r"""
        Perform applicable preliminary checks and calculations required for
        generation
        """
        self._Nx = 5
        self._Ny = 5
        self._Nz = 5
        self._Lc = 1
        self._Lx = sp.float16(self._Nx*self._Lc)
        self._Ly = sp.float16(self._Ny*self._Lc)
        self._Lz = sp.float16(self._Nz*self._Lc)

    def _generate_pores(self):
        r"""
        Generate the pores (coordinates, numbering and types)
        """
        Nx = self._Nx
        Ny = self._Ny
        Nz = self._Nz
        Lc = self._Lc
        Np = Nx*Ny*Nz
        ind = sp.arange(0, Np)
        self['pore.all'] = sp.ones_like(ind, dtype=bool)
        unraveled_index = sp.unravel_index(ind, dims=(Nx, Ny, Nz), order='F')
        pore_coords = Lc/2+Lc*sp.array(unraveled_index, dtype=sp.float64).T
        self['pore.coords'] = pore_coords

    def _generate_throats(self):
        r"""
        Generate the throats (connections, numbering and types)
        """
        Nx = self._Nx
        Ny = self._Ny
        Nz = self._Nz
        Np = Nx*Ny*Nz
        ind = sp.arange(0, Np)
        # Generate throats based on pattern of the adjacency matrix
        tpore1_1 = ind[(ind % Nx) < (Nx-1)]
        tpore2_1 = tpore1_1 + 1
        tpore1_2 = ind[(ind % (Nx*Ny)) < (Nx*(Ny-1))]
        tpore2_2 = tpore1_2 + Nx
        tpore1_3 = ind[(ind % Np) < (Nx*Ny*(Nz-1))]
        tpore2_3 = tpore1_3 + Nx*Ny
        tpore1 = sp.hstack((tpore1_1, tpore1_2, tpore1_3))
        tpore2 = sp.hstack((tpore2_1, tpore2_2, tpore2_3))
        connections = sp.vstack((tpore1, tpore2)).T
        connections = connections[sp.lexsort((connections[:, 1], connections[:, 0]))]
        self['throat.all'] = sp.ones_like(sp.arange(0, sp.shape(tpore1)[0]),
                                          dtype=bool)
        self['throat.conns'] = connections

    def _add_labels(self):
        coords = self['pore.coords']
        self['pore.front'] = self.tomask(coords[:, 0] <= self._Lc)
        self['pore.left'] = self.tomask(coords[:, 1] <= self._Lc)
        self['pore.bottom'] = self.tomask(coords[:, 2] <= self._Lc)
        self['pore.back'] = self.tomask(coords[:, 0] >= (self._Lc*(self._Nx-1)))
        self['pore.right'] = self.tomask(coords[:, 1] >= (self._Lc*(self._Ny-1)))
        self['pore.top'] = self.tomask(coords[:, 2] >= (self._Lc*(self._Nz-1)))
        for item in ['top', 'bottom', 'left', 'right', 'front', 'back']:
            ps = self.pores(item)
            ts = self.find_neighbor_throats(ps)
            ps = self.find_connected_pores(ts)
            ps0 = self['pore.'+item][ps[:, 0]]
            ps1 = self['pore.'+item][ps[:, 1]]
            ts = ts[ps1*ps0]
            self['throat.'+item] = self.tomask(throats=ts)
