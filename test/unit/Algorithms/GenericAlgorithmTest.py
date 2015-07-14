import OpenPNM
import numpy as np


class GenericAlgorithmTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[5, 5, 5])
        self.phase = OpenPNM.Phases.GenericPhase(network=self.net)
        self.alg = OpenPNM.Algorithms.GenericLinearTransport(network=self.net,
                                                             phase=self.phase)

    def test_set_BC_modes(self):
        BC1_pores = np.arange(25, 35)
        BC1_pores = np.arange(25, 35)
        self.alg.set_boundary_conditions(component=self.phase,
                                         bctype='Dirichlet',
                                         bcvalue=0.8,
                                         pores=BC1_pores)
        ptest = self.alg.pores('pore.' + self.phase.name + '_Dirichlet')
        assert np.all(ptest == BC1_pores)
        BC2_pores = np.arange(43, 50)
        self.alg.set_boundary_conditions(component=self.phase,
                                         bctype='Dirichlet',
                                         bcvalue=0.8,
                                         pores=BC2_pores,
                                         mode='merge')
        ptest = self.alg.pores('pore.' + self.phase.name + '_Dirichlet')
        assert np.all(ptest == np.concatenate((BC1_pores, BC2_pores)))
        BC3_pores = np.arange(4, 9)
        self.alg.set_boundary_conditions(bctype='Dirichlet',
                                         bcvalue=0.8,
                                         pores=BC3_pores,
                                         mode='overwrite')
        ptest = self.alg.pores('pore.' + self.phase.name + '_Dirichlet')
        assert np.all(ptest == BC3_pores)
        BC4_pores = [11, 90]
        self.alg.set_boundary_conditions(bctype='Neumann',
                                         bcvalue=0.5,
                                         pores=BC4_pores,
                                         mode='overwrite')
        ptest = self.alg.pores('pore.' + self.phase.name + '_Neumann')
        assert np.all(ptest == BC4_pores)
        self.alg.set_boundary_conditions(bctype='Dirichlet',
                                         pores=BC1_pores,
                                         bcvalue=0.3)
        ptest = self.alg.pores('pore.' + self.phase.name + '_Dirichlet')
        self.alg.set_boundary_conditions(bctype='Dirichlet',
                                         pores=self.alg.Ps,
                                         mode='remove')
        Dp = np.sum(self.alg['pore.' + self.phase.name + '_Dirichlet'])
        assert Dp == 0
        self.alg.set_boundary_conditions(bctype='Neumann',
                                         mode='remove')
        label = 'pore.' + self.phase.name + '_Neumann'
        assert (label not in self.alg.labels())
