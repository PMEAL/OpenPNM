import OpenPNM
import numpy as np
import OpenPNM.Physics.models as pm


class GenericSourceTermTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[5, 5, 5])
        self.phase = OpenPNM.Phases.GenericPhase(network=self.net)
        Ps = self.net.Ps
        Ts = self.net.Ts
        self.phys = OpenPNM.Physics.GenericPhysics(network=self.net,
                                                   phase=self.phase,
                                                   pores=Ps, throats=Ts)
        self.phys['throat.diffusive_conductance'] = 5e-8
        self.phase['pore.mole_fraction'] = 0.
        self.alg = OpenPNM.Algorithms.GenericLinearTransport(network=self.net,
                                                             phase=self.phase)
        BC_pores = np.arange(20, 30)
        self.S_pores = np.arange(55, 85)
        self.alg.set_boundary_conditions(bctype='Dirichlet',
                                         bcvalue=0.4,
                                         pores=BC_pores)

    def test_linear(self):
        self.phys['pore.item1'] = 0.5e-11
        self.phys['pore.item2'] = 1.5e-12
        self.phys.models.add(propname='pore.source1',
                             model=pm.generic_source_term.linear,
                             A1='pore.item1',
                             A2='pore.item2',
                             x='mole_fraction',
                             return_rate=False,
                             regen_mode='on_demand')
        self.phys.models.add(propname='pore.source2',
                             model=pm.generic_source_term.linear,
                             A1='pore.item1',
                             A2='pore.item2',
                             x='mole_fraction',
                             return_rate=True,
                             regen_mode='on_demand')
        self.alg.set_source_term(source_name='pore.source1',
                                 pores=self.S_pores,
                                 mode='overwrite')
        self.alg.run(conductance='throat.diffusive_conductance',
                     quantity='pore.mole_fraction',
                     super_pore_conductance=None)
        self.alg.return_results()
        self.phys.regenerate(props='pore.source1')
        self.phys.regenerate(props='pore.source2')
        X = self.phase['pore.mole_fraction']
        r1 = np.round(np.sum(0.5e-11 * X[self.S_pores] + 1.5e-12), 20)
        r2 = np.round(np.sum(self.phys['pore.source2'][self.S_pores]), 20)
        r3 = np.round(self.alg.rate(pores=self.S_pores)[0], 20)
        assert r1 == r2
        assert r2 == -r3

    def test_power_law(self):
        self.phys['pore.item1'] = 0.5e-12
        self.phys['pore.item2'] = 2.5
        self.phys['pore.item3'] = -1.4e-11
        self.phys.models.add(propname='pore.source1',
                             model=pm.generic_source_term.power_law,
                             A1='pore.item1',
                             A2='pore.item2',
                             A3='pore.item3',
                             x='mole_fraction',
                             return_rate=False,
                             regen_mode='on_demand')
        self.phys.models.add(propname='pore.source2',
                             model=pm.generic_source_term.power_law,
                             A1='pore.item1',
                             A2='pore.item2',
                             A3='pore.item3',
                             x='mole_fraction',
                             return_rate=True,
                             regen_mode='on_demand')
        self.alg.set_source_term(source_name='pore.source1',
                                 pores=self.S_pores,
                                 mode='overwrite')
        self.alg.run(conductance='throat.diffusive_conductance',
                     quantity='pore.mole_fraction',
                     super_pore_conductance=None)
        self.alg.return_results()
        self.phys.regenerate(props='pore.source1')
        self.phys.regenerate(props='pore.source2')
        X = self.phase['pore.mole_fraction']
        r1 = np.round(np.sum(0.5e-12 * X[self.S_pores] ** 2.5 - 1.4e-11), 20)
        r2 = np.round(np.sum(self.phys['pore.source2'][self.S_pores]), 20)
        r3 = np.round(self.alg.rate(pores=self.S_pores)[0], 20)
        assert r1 == r2
        assert r2 == -r3

    def test_exponential(self):
        self.phys['pore.item1'] = 0.8e-11
        self.phys['pore.item2'] = 3
        self.phys['pore.item3'] = 0.5
        self.phys['pore.item4'] = 2
        self.phys['pore.item5'] = -0.34
        self.phys['pore.item6'] = 2e-14
        self.phys.models.add(propname='pore.source1',
                             model=pm.generic_source_term.exponential,
                             A1='pore.item1',
                             A2='pore.item2',
                             A3='pore.item3',
                             A4='pore.item4',
                             A5='pore.item5',
                             A6='pore.item6',
                             x='mole_fraction',
                             return_rate=False,
                             regen_mode='on_demand')
        self.phys.models.add(propname='pore.source2',
                             model=pm.generic_source_term.exponential,
                             A1='pore.item1',
                             A2='pore.item2',
                             A3='pore.item3',
                             A4='pore.item4',
                             A5='pore.item5',
                             A6='pore.item6',
                             x='mole_fraction',
                             return_rate=True,
                             regen_mode='on_demand')
        self.alg.set_source_term(source_name='pore.source1',
                                 pores=self.S_pores,
                                 mode='overwrite')
        self.alg.run(conductance='throat.diffusive_conductance',
                     quantity='pore.mole_fraction',
                     super_pore_conductance=None)
        self.alg.return_results()
        self.phys.regenerate(props='pore.source1')
        self.phys.regenerate(props='pore.source2')
        X = self.phase['pore.mole_fraction']
        r1 = np.round(np.sum(0.8e-11 * 3 ** (0.5 * X[self.S_pores] ** 2 -
                      0.34) + 2e-14), 20)
        r2 = np.round(np.sum(self.phys['pore.source2'][self.S_pores]), 20)
        r3 = np.round(self.alg.rate(pores=self.S_pores)[0], 20)
        assert r1 == r2
        assert r2 == -r3

    def test_natural_exponential(self):
        self.phys['pore.item1'] = 0.8e-11
        self.phys['pore.item2'] = 0.5
        self.phys['pore.item3'] = 2
        self.phys['pore.item4'] = -0.34
        self.phys['pore.item5'] = 2e-14
        self.phys.models.add(propname='pore.source1',
                             model=pm.generic_source_term.natural_exponential,
                             A1='pore.item1',
                             A2='pore.item2',
                             A3='pore.item3',
                             A4='pore.item4',
                             A5='pore.item5',
                             x='mole_fraction',
                             return_rate=False,
                             regen_mode='on_demand')
        self.phys.models.add(propname='pore.source2',
                             model=pm.generic_source_term.natural_exponential,
                             A1='pore.item1',
                             A2='pore.item2',
                             A3='pore.item3',
                             A4='pore.item4',
                             A5='pore.item5',
                             x='mole_fraction',
                             return_rate=True,
                             regen_mode='on_demand')
        self.alg.set_source_term(source_name='pore.source1',
                                 pores=self.S_pores,
                                 mode='overwrite')
        self.alg.run(conductance='throat.diffusive_conductance',
                     quantity='pore.mole_fraction',
                     super_pore_conductance=None)
        self.alg.return_results()
        self.phys.regenerate(props='pore.source1')
        self.phys.regenerate(props='pore.source2')
        X = self.phase['pore.mole_fraction']
        r1 = np.round(np.sum(0.8e-11 * np.exp(0.5 * X[self.S_pores] ** 2 -
                      0.34) + 2e-14), 20)
        r2 = np.round(np.sum(self.phys['pore.source2'][self.S_pores]), 20)
        r3 = np.round(self.alg.rate(pores=self.S_pores)[0], 20)
        assert r1 == r2
        assert r2 == -r3

    def test_logarithm(self):
        self.phys['pore.item1'] = 0.16e-13
        self.phys['pore.item2'] = 10
        self.phys['pore.item3'] = 4
        self.phys['pore.item4'] = 1.4
        self.phys['pore.item5'] = 0.133
        self.phys['pore.item6'] = -5.1e-13
        self.phys.models.add(propname='pore.source1',
                             model=pm.generic_source_term.logarithm,
                             A1='pore.item1',
                             A2='pore.item2',
                             A3='pore.item3',
                             A4='pore.item4',
                             A5='pore.item5',
                             A6='pore.item6',
                             x='mole_fraction',
                             return_rate=False,
                             regen_mode='on_demand')
        self.phys.models.add(propname='pore.source2',
                             model=pm.generic_source_term.logarithm,
                             A1='pore.item1',
                             A2='pore.item2',
                             A3='pore.item3',
                             A4='pore.item4',
                             A5='pore.item5',
                             A6='pore.item6',
                             x='mole_fraction',
                             return_rate=True,
                             regen_mode='on_demand')
        self.alg.set_source_term(source_name='pore.source1',
                                 pores=self.S_pores,
                                 mode='overwrite')
        self.alg.run(conductance='throat.diffusive_conductance',
                     quantity='pore.mole_fraction',
                     super_pore_conductance=None)
        self.alg.return_results()
        self.phys.regenerate(props='pore.source1')
        self.phys.regenerate(props='pore.source2')
        X = self.phase['pore.mole_fraction']
        r1 = np.round(np.sum(0.16e-13 * np.log(4 * X[self.S_pores] ** (1.4) +
                             0.133) / np.log(10) - 5.1e-13), 20)
        r2 = np.round(np.sum(self.phys['pore.source2'][self.S_pores]), 20)
        r3 = np.round(self.alg.rate(pores=self.S_pores)[0], 20)
        assert r1 == r2
        assert r2 == -r3

    def test_natural_logarithm(self):
        self.phys['pore.item1'] = 0.16e-14
        self.phys['pore.item2'] = 4
        self.phys['pore.item3'] = 1.4
        self.phys['pore.item4'] = 0.133
        self.phys['pore.item5'] = -5.1e-14
        self.phys.models.add(propname='pore.source1',
                             model=pm.generic_source_term.natural_logarithm,
                             A1='pore.item1',
                             A2='pore.item2',
                             A3='pore.item3',
                             A4='pore.item4',
                             A5='pore.item5',
                             x='mole_fraction',
                             return_rate=False,
                             regen_mode='on_demand')
        self.phys.models.add(propname='pore.source2',
                             model=pm.generic_source_term.natural_logarithm,
                             A1='pore.item1',
                             A2='pore.item2',
                             A3='pore.item3',
                             A4='pore.item4',
                             A5='pore.item5',
                             x='mole_fraction',
                             return_rate=True,
                             regen_mode='on_demand')
        self.alg.set_source_term(source_name='pore.source1',
                                 pores=self.S_pores,
                                 mode='overwrite')
        self.alg.run(conductance='throat.diffusive_conductance',
                     quantity='pore.mole_fraction',
                     super_pore_conductance=None)
        self.alg.return_results()
        self.phys.regenerate(props='pore.source1')
        self.phys.regenerate(props='pore.source2')
        X = self.phase['pore.mole_fraction']
        r1 = np.round(np.sum(0.16e-14 * np.log(4 * X[self.S_pores] ** (1.4) +
                             0.133) - 5.1e-14), 20)
        r2 = np.round(np.sum(self.phys['pore.source2'][self.S_pores]), 20)
        r3 = np.round(self.alg.rate(pores=self.S_pores)[0], 20)
        assert r1 == r2
        assert r2 == -r3
