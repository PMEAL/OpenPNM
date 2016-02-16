import OpenPNM
import numpy as np
import OpenPNM.Physics.models as pm


class GenericLinearTransportTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[5, 5, 5])
        self.phase = OpenPNM.Phases.GenericPhase(network=self.net)
        Ps = self.net.Ps
        Ts = self.net.Ts
        self.phys = OpenPNM.Physics.GenericPhysics(network=self.net,
                                                   phase=self.phase,
                                                   pores=Ps, throats=Ts)
        self.phys['throat.cond'] = 5e-8
        self.alg = OpenPNM.Algorithms.GenericLinearTransport(network=self.net,
                                                             phase=self.phase)

    def test_set_BC_modes(self):
        BC1_pores = np.arange(25, 35)
        BC1_pores = np.arange(25, 35)
        self.alg.set_boundary_conditions(bctype='Dirichlet',
                                         bcvalue=0.8,
                                         pores=BC1_pores)
        ptest = self.alg.pores('pore.Dirichlet')
        assert np.all(ptest == BC1_pores)
        BC2_pores = np.arange(43, 50)
        self.alg.set_boundary_conditions(bctype='Dirichlet',
                                         bcvalue=0.8,
                                         pores=BC2_pores,
                                         mode='merge')
        ptest = self.alg.pores('pore.Dirichlet')
        assert np.all(ptest == np.concatenate((BC1_pores, BC2_pores)))
        BC3_pores = np.arange(4, 9)
        self.alg.set_boundary_conditions(bctype='Dirichlet',
                                         bcvalue=0.8,
                                         pores=BC3_pores,
                                         mode='overwrite')
        ptest = self.alg.pores('pore.Dirichlet')
        assert np.all(ptest == BC3_pores)
        BC4_pores = [11, 90]
        self.alg.set_boundary_conditions(bctype='Neumann',
                                         bcvalue=0.5,
                                         pores=BC4_pores,
                                         mode='overwrite')
        ptest = self.alg.pores('pore.Neumann')
        assert np.all(ptest == BC4_pores)
        self.alg.set_boundary_conditions(bctype='Dirichlet',
                                         pores=BC1_pores,
                                         bcvalue=0.3)
        ptest = self.alg.pores('pore.Dirichlet')
        self.alg.set_boundary_conditions(bctype='Dirichlet',
                                         pores=self.alg.Ps,
                                         mode='remove')
        Dp = np.sum(self.alg['pore.Dirichlet'])
        assert Dp == 0
        self.alg.set_boundary_conditions(bctype='Neumann',
                                         mode='remove')
        label = 'pore.Neumann'
        assert (label not in self.alg.labels())

    def test_super_pore_conductance(self):
        g_super = []
        BC1_pores = np.arange(20, 30)
        self.alg.set_boundary_conditions(bctype='Dirichlet',
                                         bcvalue=0.4,
                                         pores=BC1_pores)
        BC2_pores = np.arange(45, 66)
        self.alg.set_boundary_conditions(bctype='Neumann_group',
                                         bcvalue=1.4e-10,
                                         pores=BC2_pores)
        g_super.append(2e-12)
        BC3_pores = np.arange(87, 94)
        self.alg.set_boundary_conditions(bctype='Neumann_group',
                                         bcvalue=-0.9e-10,
                                         pores=BC3_pores)
        g_super.append(np.ones(len(BC3_pores)) * 1.5e-12)
        BC4_pores = np.arange(3, 7)
        self.alg.set_boundary_conditions(bctype='Neumann_group',
                                         bcvalue=0.1e-10,
                                         pores=BC4_pores)
        g_super.append(np.array([6.42e-13]))
        self.alg.run(conductance='throat.cond',
                     quantity='pore.mole_fraction',
                     super_pore_conductance=g_super)
        self.alg.return_results()
        r1 = self.alg.rate(BC1_pores)[0]
        r2 = self.alg.rate(BC2_pores)[0]
        r3 = self.alg.rate(BC3_pores)[0]
        r4 = self.alg.rate(BC4_pores)[0]
        assert np.absolute(r1 + r2 + r3 + r4) < 1e-20
        assert np.size(self.alg.super_pore_conductance[0]) == 1
        assert np.size(self.alg.super_pore_conductance[1]) == 7
        assert np.size(self.alg.super_pore_conductance[2]) == 1

    def test_source_term_modes(self):
        self.phys['pore.item1'] = 0.5e-12
        self.phys['pore.item2'] = 2.5
        self.phys['pore.item3'] = -1.4e-11
        self.phys.models.add(propname='pore.A',
                             model=pm.generic_source_term.power_law,
                             A1='pore.item1',
                             A2='pore.item2',
                             A3='pore.item3',
                             x='mole_fraction',
                             return_rate=False,
                             regen_mode='on_demand')
        self.phys.models.add(propname='pore.B',
                             model=pm.generic_source_term.linear,
                             A1='pore.item1',
                             A2='pore.item3',
                             x='mole_fraction',
                             return_rate=False,
                             regen_mode='on_demand')
        S1_pores = np.arange(25, 35)
        self.alg.set_source_term(source_name=['pore.A', 'pore.B'],
                                 pores=S1_pores)
        mask1 = ~np.isnan(self.alg['pore.source_nonlinear_s1_A'])
        mask2 = ~np.isnan(self.alg['pore.source_nonlinear_s2_A'])
        assert np.all(self.alg.Ps[mask1] == S1_pores)
        assert np.all(self.alg.Ps[mask2] == S1_pores)
        self.alg.set_source_term(source_name='pore.A',
                                 pores=[26], x0=np.ones(self.phys.Np),
                                 mode='update')
        assert self.alg['pore.source_nonlinear_s1_A'][26] == 1.25e-12
        S2_pores = np.array([30, 31])
        self.alg.set_source_term(source_name='pore.A',
                                 pores=S2_pores,
                                 mode='overwrite')
        mask1 = ~np.isnan(self.alg['pore.source_nonlinear_s1_A'])
        assert np.all(self.alg.Ps[mask1] == S2_pores)
        self.alg.set_source_term(source_name='pore.B',
                                 pores=S1_pores,
                                 mode='remove')
        mask1 = np.isnan(self.alg['pore.source_nonlinear_s1_B'])
        assert np.all(self.alg.Ps[mask1] == self.alg.Ps)
        self.alg.set_source_term(source_name=['pore.A', 'pore.B'],
                                 pores=self.alg.Ps,
                                 mode='remove')
        assert ('pore.source_B' in self.alg.labels())
        assert ('pore.source_A' in self.alg.labels())
        self.alg.set_source_term(source_name=['pore.A', 'pore.B'],
                                 mode='remove')
        assert ('pore.source_B' not in self.alg.labels())
        assert ('pore.source_A' not in self.alg.labels())
