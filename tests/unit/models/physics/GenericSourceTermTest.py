import collections
import numpy as np
import openpnm as op
import openpnm.models.physics as pm
from numpy.testing import assert_allclose
from sympy import ln as sym_ln
from sympy import symbols


class GenericSourceTermTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5])
        Ps = self.net.Ps
        Ts = self.net.Ts
        self.geo = op.geometry.GenericGeometry(network=self.net, pores=Ps,
                                               throats=Ts)
        self.phase = op.phase.GenericPhase(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)
        self.phys['throat.diffusive_conductance'] = 5e-8
        self.phase['pore.mole_fraction'] = 0.0
        self.BC_pores = np.arange(20, 30)
        self.source_pores = np.arange(55, 85)

    def test_default_values_should_give_zero_rate(self):
        sources = ["linear", "power_law", "exponential", "natural_exponential",
                   "logarithm", "natural_logarithm"]
        self.alg = op.algorithms.ReactiveTransport(network=self.net,
                                                   phase=self.phase)
        self.alg.settings._update({'conductance': 'throat.diffusive_conductance',
                                   'quantity': 'pore.mole_fraction'})
        self.alg.set_value_BC(values=0.4, pores=self.BC_pores)
        self.alg.set_source(propname='pore.source_term',
                            pores=self.source_pores)
        # To avoid nans in logarithm-based source term models
        self.phase['pore.mole_fraction'] = 0.1
        for source in sources:
            self.phys.add_model(propname="pore.source_term",
                                model=getattr(pm.source_terms, source),
                                X="pore.mole_fraction")
            assert self.phys["pore.source_term.rate"].mean() == 0
            assert self.phys["pore.source_term.S1"].mean() == 0
            assert self.phys["pore.source_term.S2"].mean() == 0

    def test_linear(self):
        self.phys['pore.item1'] = 0.5e-11
        self.phys['pore.item2'] = 1.5e-12
        self.phys.add_model(propname='pore.source1',
                            model=pm.source_terms.linear,
                            A1='pore.item1',
                            A2='pore.item2',
                            X='pore.mole_fraction',
                            regen_mode='normal')
        self.alg = op.algorithms.ReactiveTransport(network=self.net,
                                                   phase=self.phase)
        self.alg.settings._update({'conductance': 'throat.diffusive_conductance',
                                   'quantity': 'pore.mole_fraction'})
        self.alg.set_value_BC(values=0.4, pores=self.BC_pores)
        self.alg.set_source(propname='pore.source1',
                            pores=self.source_pores)
        self.alg.run()
        self.phase.update(self.alg.results())
        self.phys.regenerate_models(propnames='pore.source1')
        X = self.phase['pore.mole_fraction']
        r1 = np.round(np.sum(0.5e-11 * X[self.source_pores] + 1.5e-12), 20)
        r2 = np.round(np.sum(self.phys['pore.source1.rate'][self.source_pores]), 20)
        assert r1 == r2

    def test_power_law(self):
        self.phys['pore.item1'] = 0.5e-12
        self.phys['pore.item2'] = 2.5
        self.phys['pore.item3'] = -1.4e-11
        self.phys.add_model(propname='pore.source1',
                            model=pm.source_terms.power_law,
                            A1='pore.item1',
                            A2='pore.item2',
                            A3='pore.item3',
                            X='pore.mole_fraction',
                            regen_mode='normal')
        self.alg = op.algorithms.ReactiveTransport(network=self.net,
                                                   phase=self.phase)
        self.alg.set_value_BC(values=0.4, pores=self.BC_pores)
        self.alg.set_source(propname='pore.source1',
                            pores=self.source_pores)
        self.alg.settings._update({'conductance': 'throat.diffusive_conductance',
                                   'quantity': 'pore.mole_fraction'})
        self.alg.run()
        self.phase.update(self.alg.results())
        self.phys.regenerate_models(propnames='pore.source1')
        X = self.phase['pore.mole_fraction']
        r1 = np.sum(0.5e-12 * X[self.source_pores]**2.5 - 1.4e-11)
        r2 = np.sum(self.phys['pore.source1.rate'][self.source_pores])
        assert r1 == r2

    def test_exponential(self):
        self.phys['pore.item1'] = 0.8e-11
        self.phys['pore.item2'] = 3
        self.phys['pore.item3'] = 0.5
        self.phys['pore.item4'] = 2
        self.phys['pore.item5'] = -0.34
        self.phys['pore.item6'] = 2e-14
        self.phys.add_model(propname='pore.source1',
                            model=pm.source_terms.exponential,
                            A1='pore.item1',
                            A2='pore.item2',
                            A3='pore.item3',
                            A4='pore.item4',
                            A5='pore.item5',
                            A6='pore.item6',
                            X='pore.mole_fraction',
                            regen_mode='normal')
        self.alg = op.algorithms.ReactiveTransport(network=self.net,
                                                   phase=self.phase)
        self.alg.set_value_BC(values=0.4, pores=self.BC_pores)
        self.alg.set_source(propname='pore.source1',
                            pores=self.source_pores)
        self.alg.settings._update({'conductance': 'throat.diffusive_conductance',
                                   'quantity': 'pore.mole_fraction'})
        self.alg.run()
        self.phase.update(self.alg.results())
        self.phys.regenerate_models(propnames='pore.source1')
        X = self.phase['pore.mole_fraction']
        r1 = np.sum(0.8e-11 * 3 ** (0.5 * X[self.source_pores]**2 - 0.34) + 2e-14)
        r2 = np.sum(self.phys['pore.source1.rate'][self.source_pores])
        assert_allclose(actual=r2, desired=r1, rtol=1e-7)

    def test_natural_exponential(self):
        self.phys['pore.item1'] = 0.8e-11
        self.phys['pore.item2'] = 0.5
        self.phys['pore.item3'] = 2
        self.phys['pore.item4'] = -0.34
        self.phys['pore.item5'] = 2e-14
        self.phys.add_model(propname='pore.source1',
                            model=pm.source_terms.natural_exponential,
                            A1='pore.item1',
                            A2='pore.item2',
                            A3='pore.item3',
                            A4='pore.item4',
                            A5='pore.item5',
                            X='pore.mole_fraction',
                            regen_mode='normal')
        self.alg = op.algorithms.ReactiveTransport(network=self.net,
                                                   phase=self.phase)
        self.alg.set_value_BC(values=0.4, pores=self.BC_pores)
        self.alg.set_source(propname='pore.source1',
                            pores=self.source_pores)
        self.alg.settings._update({'conductance': 'throat.diffusive_conductance',
                                   'quantity': 'pore.mole_fraction'})
        self.alg.run()
        self.phase.update(self.alg.results())
        self.phys.regenerate_models(propnames='pore.source1')
        X = self.phase['pore.mole_fraction']
        r1 = np.sum(0.8e-11 * np.exp(0.5 * X[self.source_pores]**2 - 0.34) + 2e-14)
        r2 = np.sum(self.phys['pore.source1.rate'][self.source_pores])
        assert_allclose(actual=r2, desired=r1)

    def test_logarithm(self):
        self.phys['pore.item1'] = 0.16e-13
        self.phys['pore.item2'] = 10
        self.phys['pore.item3'] = 4
        self.phys['pore.item4'] = 1.4
        self.phys['pore.item5'] = 0.133
        self.phys['pore.item6'] = -5.1e-13
        self.phys.add_model(propname='pore.source1',
                            model=pm.source_terms.logarithm,
                            A1='pore.item1',
                            A2='pore.item2',
                            A3='pore.item3',
                            A4='pore.item4',
                            A5='pore.item5',
                            A6='pore.item6',
                            X='pore.mole_fraction',
                            regen_mode='normal')
        self.alg = op.algorithms.ReactiveTransport(network=self.net,
                                                   phase=self.phase)
        self.alg.set_value_BC(values=0.4, pores=self.BC_pores)
        self.alg.set_source(propname='pore.source1',
                            pores=self.source_pores)
        self.alg.settings._update({'conductance': 'throat.diffusive_conductance',
                                   'quantity': 'pore.mole_fraction'})
        self.alg.run()
        self.phase.update(self.alg.results())
        self.phys.regenerate_models(propnames='pore.source1')
        X = self.phase['pore.mole_fraction']
        r1 = np.sum(0.16e-13 * np.log(4*X[self.source_pores]**(1.4) + 0.133)
                    / np.log(10) - 5.1e-13)
        r2 = np.sum(self.phys['pore.source1.rate'][self.source_pores])
        assert_allclose(actual=r2, desired=r1)

    def test_natural_logarithm(self):
        self.phys['pore.item1'] = 0.16e-14
        self.phys['pore.item2'] = 4
        self.phys['pore.item3'] = 1.4
        self.phys['pore.item4'] = 0.133
        self.phys['pore.item5'] = -5.1e-14
        self.phys.add_model(propname='pore.source1',
                            model=pm.source_terms.natural_logarithm,
                            A1='pore.item1',
                            A2='pore.item2',
                            A3='pore.item3',
                            A4='pore.item4',
                            A5='pore.item5',
                            X='pore.mole_fraction',
                            regen_mode='on_demand')
        self.alg = op.algorithms.ReactiveTransport(network=self.net,
                                                   phase=self.phase)
        self.alg.set_value_BC(values=0.4, pores=self.BC_pores)
        self.alg.set_source(propname='pore.source1',
                            pores=self.source_pores)
        self.alg.settings._update({'conductance': 'throat.diffusive_conductance',
                                   'quantity': 'pore.mole_fraction'})
        self.alg.run()
        self.phase.update(self.alg.results())
        self.phys.regenerate_models(propnames='pore.source1')
        X = self.phase['pore.mole_fraction']
        r1 = np.sum(0.16e-14*np.log(4*X[self.source_pores]**1.4 + 0.133) - 5.1e-14)
        r2 = np.sum(self.phys['pore.source1.rate'][self.source_pores])
        assert r1 == r2

    def test_general_symbolic(self):
        # natural log function
        y = "a*ln(b*x**c + d)+e"
        phys = self.phys
        phys['pore.item1'] = 0.16e-14
        phys['pore.item2'] = 4
        phys['pore.item3'] = 1.4
        phys['pore.item4'] = 0.133
        phys['pore.item5'] = -5.1e-14
        phys.add_model(propname='pore.source1',
                       model=pm.source_terms.natural_logarithm,
                       A1='pore.item1',
                       A2='pore.item2',
                       A3='pore.item3',
                       A4='pore.item4',
                       A5='pore.item5',
                       X='pore.mole_fraction',
                       regen_mode='normal')
        arg_map = {'a': 'pore.item1',
                   'b': 'pore.item2',
                   'c': 'pore.item3',
                   'd': 'pore.item4',
                   'e': 'pore.item5'}
        phys.add_model(propname='pore.general',
                       model=op.models.physics.generic_source_term.general_symbolic,
                       eqn=y, x='pore.mole_fraction', **arg_map)
        assert np.allclose(phys['pore.source1.rate'], phys['pore.general.rate'])
        assert np.allclose(phys['pore.source1.S1'], phys['pore.general.S1'])
        assert np.allclose(phys['pore.source1.S2'], phys['pore.general.S2'])

    def test_butler_volmer_kinetics(self):
        np.random.seed(10)
        self.net["pore.reaction_area"] = np.random.rand(self.net.Np)
        self.phys['pore.electrolyte_voltage'] = np.random.rand(self.net.Np) * 0.2
        self.phys['pore.solid_voltage'] = 1.1
        self.phys['pore.open_circuit_voltage'] = 1.2
        self.phase['pore.electrolyte_concentration'] = np.random.rand(self.net.Np)
        BV_params = {
            "z": 4,
            "j0": 1e-3,
            "c_ref": 1000,
            "alpha_anode": 0.4,
            "alpha_cathode": 0.6
        }
        self.phys.add_model(propname='pore.rxn_BV_c',
                            model=pm.source_terms.butler_volmer_conc,
                            X="pore.electrolyte_concentration", **BV_params)
        self.phys.add_model(propname='pore.rxn_BV_v',
                            model=pm.source_terms.butler_volmer_voltage,
                            X="pore.electrolyte_voltage", **BV_params)
        # Check Butler-Volmer model (concentration)
        S1_BV_c = self.phys["pore.rxn_BV_c.S1"]
        S2_BV_c = self.phys["pore.rxn_BV_c.S2"]
        rate_BV_c = self.phys["pore.rxn_BV_c.rate"]
        eps = np.finfo(float).eps * np.abs(rate_BV_c).max()
        assert_allclose(S2_BV_c.mean(), 0, atol=eps)
        assert_allclose(S1_BV_c.mean(), -0.06977055)
        assert_allclose(rate_BV_c.mean(), -0.03541646)
        # Check Butler-Volmer model (voltage)
        S1_BV_v = self.phys["pore.rxn_BV_v.S1"]
        S2_BV_v = self.phys["pore.rxn_BV_v.S2"]
        rate_BV_v = self.phys["pore.rxn_BV_v.rate"]
        assert_allclose(S2_BV_v.mean(), 226884.32)
        assert_allclose(S1_BV_v.mean(), -1277463.52)
        assert_allclose(rate_BV_v.mean(), -13668.675)
        # The two Butler-Volmer models must only differ by z*F (unit conversion)
        assert_allclose(rate_BV_v, rate_BV_c * BV_params["z"] * 96485.33212)


if __name__ == '__main__':

    t = GenericSourceTermTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
