import openpnm as op
import numpy as np
import openpnm.models.physics as pm


class GenericSourceTermTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5])
        Ps = self.net.Ps
        Ts = self.net.Ts
        self.geo = op.geometry.GenericGeometry(network=self.net, pores=Ps,
                                               throats=Ts)
        self.phase = op.phases.GenericPhase(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)
        self.phys['throat.diffusive_conductance'] = 5e-8
        self.phase['pore.mole_fraction'] = 0.0
        self.BC_pores = np.arange(20, 30)
        self.source_pores = np.arange(55, 85)

    def test_linear(self):
        self.phys['pore.item1'] = 0.5e-11
        self.phys['pore.item2'] = 1.5e-12
        self.phys.add_model(propname='pore.source1',
                            model=pm.generic_source_term.linear,
                            A1='pore.item1',
                            A2='pore.item2',
                            quantity='pore.mole_fraction',
                            regen_mode='normal')
        self.alg = op.algorithms.ReactiveTransport(network=self.net,
                                                   phase=self.phase)
        self.alg.settings.update({'conductance': 'throat.diffusive_conductance',
                                  'quantity': 'pore.mole_fraction'})
        self.alg.set_BC(bctype='dirichlet', bcvalues=0.4,
                        pores=self.BC_pores)
        self.alg.set_source(propname='pore.source1',
                            pores=self.source_pores)
        self.alg.run()
        self.phase.update(self.alg.results())
        self.phys.regenerate_models(propnames='pore.source1')
        X = self.phase['pore.mole_fraction']
        r1 = np.round(np.sum(0.5e-11 * X[self.source_pores] + 1.5e-12), 20)
        r2 = np.round(np.sum(self.phys['pore.source1.rate'][self.source_pores]), 20)
#        r3 = np.round(self.alg.rate(pores=self.S_pores)[0], 20)
        assert r1 == r2
#        assert r2 == -r3

    def test_power_law(self):
        self.phys['pore.item1'] = 0.5e-12
        self.phys['pore.item2'] = 2.5
        self.phys['pore.item3'] = -1.4e-11
        self.phys.add_model(propname='pore.source1',
                            model=pm.generic_source_term.power_law,
                            A1='pore.item1',
                            A2='pore.item2',
                            A3='pore.item3',
                            quantity='pore.mole_fraction',
                            regen_mode='normal')
        self.alg = op.algorithms.ReactiveTransport(network=self.net,
                                                   phase=self.phase)
        self.alg.set_BC(bctype='dirichlet', bcvalues=0.4,
                        pores=self.BC_pores)
        self.alg.set_source(propname='pore.source1',
                            pores=self.source_pores)
        self.alg.settings.update({'conductance': 'throat.diffusive_conductance',
                                  'quantity': 'pore.mole_fraction'})
        self.alg.run()
        self.phase.update(self.alg.results())
        self.phys.regenerate_models(propnames='pore.source1')
        X = self.phase['pore.mole_fraction']
        r1 = np.sum(0.5e-12 * X[self.source_pores]**2.5 - 1.4e-11)
        r2 = np.sum(self.phys['pore.source1.rate'][self.source_pores])
#        r3 = np.round(self.alg.rate(pores=self.S_pores)[0], 20)
        assert r1 == r2
#        assert r2 == -r3

    def test_exponential(self):
        self.phys['pore.item1'] = 0.8e-11
        self.phys['pore.item2'] = 3
        self.phys['pore.item3'] = 0.5
        self.phys['pore.item4'] = 2
        self.phys['pore.item5'] = -0.34
        self.phys['pore.item6'] = 2e-14
        self.phys.add_model(propname='pore.source1',
                            model=pm.generic_source_term.exponential,
                            A1='pore.item1',
                            A2='pore.item2',
                            A3='pore.item3',
                            A4='pore.item4',
                            A5='pore.item5',
                            A6='pore.item6',
                            quantity='pore.mole_fraction',
                            regen_mode='normal')
        self.alg = op.algorithms.ReactiveTransport(network=self.net,
                                                   phase=self.phase)
        self.alg.set_BC(bctype='dirichlet', bcvalues=0.4,
                        pores=self.BC_pores)
        self.alg.set_source(propname='pore.source1',
                            pores=self.source_pores)
        self.alg.settings.update({'conductance': 'throat.diffusive_conductance',
                                  'quantity': 'pore.mole_fraction'})
        self.alg.run()
        self.phase.update(self.alg.results())
        self.phys.regenerate_models(propnames='pore.source1')
        X = self.phase['pore.mole_fraction']
        r1 = np.round(np.sum(0.8e-11 * 3 ** (0.5 * X[self.source_pores]**2 -
                      0.34) + 2e-14), 20)
        r2 = np.round(np.sum(self.phys['pore.source1.rate'][self.source_pores]), 20)
#        r3 = np.round(self.alg.rate(pores=self.S_pores)[0], 20)
        assert r1 == r2
#        assert r2 == -r3

    def test_natural_exponential(self):
        self.phys['pore.item1'] = 0.8e-11
        self.phys['pore.item2'] = 0.5
        self.phys['pore.item3'] = 2
        self.phys['pore.item4'] = -0.34
        self.phys['pore.item5'] = 2e-14
        self.phys.add_model(propname='pore.source1',
                            model=pm.generic_source_term.natural_exponential,
                            A1='pore.item1',
                            A2='pore.item2',
                            A3='pore.item3',
                            A4='pore.item4',
                            A5='pore.item5',
                            quantity='pore.mole_fraction',
                            regen_mode='normal')
        self.alg = op.algorithms.ReactiveTransport(network=self.net,
                                                   phase=self.phase)
        self.alg.set_BC(bctype='dirichlet', bcvalues=0.4,
                        pores=self.BC_pores)
        self.alg.set_source(propname='pore.source1',
                            pores=self.source_pores)
        self.alg.settings.update({'conductance': 'throat.diffusive_conductance',
                                  'quantity': 'pore.mole_fraction'})
        self.alg.run()
        self.phase.update(self.alg.results())
        self.phys.regenerate_models(propnames='pore.source1')
        X = self.phase['pore.mole_fraction']
        r1 = np.round(np.sum(0.8e-11 * np.exp(0.5 * X[self.source_pores]**2 -
                      0.34) + 2e-14), 20)
        r2 = np.round(np.sum(self.phys['pore.source1.rate'][self.source_pores]), 20)
#        r3 = np.round(self.alg.rate(pores=self.S_pores)[0], 20)
        assert r1 == r2
#        assert r2 == -r3

    def test_logarithm(self):
        self.phys['pore.item1'] = 0.16e-13
        self.phys['pore.item2'] = 10
        self.phys['pore.item3'] = 4
        self.phys['pore.item4'] = 1.4
        self.phys['pore.item5'] = 0.133
        self.phys['pore.item6'] = -5.1e-13
        self.phys.add_model(propname='pore.source1',
                            model=pm.generic_source_term.logarithm,
                            A1='pore.item1',
                            A2='pore.item2',
                            A3='pore.item3',
                            A4='pore.item4',
                            A5='pore.item5',
                            A6='pore.item6',
                            quantity='pore.mole_fraction',
                            regen_mode='normal')
        self.alg = op.algorithms.ReactiveTransport(network=self.net,
                                                   phase=self.phase)
        self.alg.set_BC(bctype='dirichlet', bcvalues=0.4,
                        pores=self.BC_pores)
        self.alg.set_source(propname='pore.source1',
                            pores=self.source_pores)
        self.alg.settings.update({'conductance': 'throat.diffusive_conductance',
                                  'quantity': 'pore.mole_fraction'})
        self.alg.run()
        self.phase.update(self.alg.results())
        self.phys.regenerate_models(propnames='pore.source1')
        X = self.phase['pore.mole_fraction']
        r1 = np.round(np.sum(0.16e-13*np.log(4*X[self.source_pores]**(1.4) +
                             0.133) / np.log(10) - 5.1e-13), 20)
        r2 = np.round(np.sum(self.phys['pore.source1.rate'][self.source_pores]), 20)
#        r3 = np.round(self.alg.rate(pores=self.S_pores)[0], 20)
        assert r1 == r2
#        assert r2 == -r3

    def test_natural_logarithm(self):
        self.phys['pore.item1'] = 0.16e-14
        self.phys['pore.item2'] = 4
        self.phys['pore.item3'] = 1.4
        self.phys['pore.item4'] = 0.133
        self.phys['pore.item5'] = -5.1e-14
        self.phys.add_model(propname='pore.source1',
                            model=pm.generic_source_term.natural_logarithm,
                            A1='pore.item1',
                            A2='pore.item2',
                            A3='pore.item3',
                            A4='pore.item4',
                            A5='pore.item5',
                            quantity='pore.mole_fraction',
                            regen_mode='on_demand')
        self.alg = op.algorithms.ReactiveTransport(network=self.net,
                                                   phase=self.phase)
        self.alg.set_BC(bctype='dirichlet', bcvalues=0.4,
                        pores=self.BC_pores)
        self.alg.set_source(propname='pore.source1',
                            pores=self.source_pores)
        self.alg.settings.update({'conductance': 'throat.diffusive_conductance',
                                  'quantity': 'pore.mole_fraction'})
        self.alg.run()
        self.phase.update(self.alg.results())
        self.phys.regenerate_models(propnames='pore.source1')
        X = self.phase['pore.mole_fraction']
        r1 = np.round(np.sum(0.16e-14*np.log(4*X[self.source_pores]**(1.4) +
                             0.133) - 5.1e-14), 20)
        r2 = np.round(np.sum(self.phys['pore.source1.rate'][self.source_pores]), 20)
#        r3 = np.round(self.alg.rate(pores=self.S_pores)[0], 20)
        assert r1 == r2
#        assert r2 == -r3


if __name__ == '__main__':

    t = GenericSourceTermTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
