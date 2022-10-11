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
        self.phase = op.phase.Phase(network=self.net)
        self.phase['throat.diffusive_conductance'] = 5e-8
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
            self.phase.add_model(propname="pore.source_term",
                                 model=getattr(pm.source_terms, source),
                                 X="pore.mole_fraction")
            assert self.phase["pore.source_term.rate"].mean() == 0
            assert self.phase["pore.source_term.S1"].mean() == 0
            assert self.phase["pore.source_term.S2"].mean() == 0

    def test_linear(self):
        self.phase['pore.item1'] = 0.5e-11
        self.phase['pore.item2'] = 1.5e-12
        self.phase.add_model(propname='pore.source1',
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
        self.phase.update(self.alg.soln)
        self.phase.regenerate_models(propnames='pore.source1')
        X = self.phase['pore.mole_fraction']
        r1 = np.round(np.sum(0.5e-11 * X[self.source_pores] + 1.5e-12), 20)
        r2 = np.round(np.sum(self.phase['pore.source1.rate'][self.source_pores]), 20)
        assert r1 == r2

    def test_power_law(self):
        self.phase['pore.item1'] = 0.5e-12
        self.phase['pore.item2'] = 2.5
        self.phase['pore.item3'] = -1.4e-11
        self.phase.add_model(propname='pore.source1',
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
        self.phase.update(self.alg.soln)
        self.phase.regenerate_models(propnames='pore.source1')
        X = self.phase['pore.mole_fraction']
        r1 = np.sum(0.5e-12 * X[self.source_pores]**2.5 - 1.4e-11)
        r2 = np.sum(self.phase['pore.source1.rate'][self.source_pores])
        assert r1 == r2

    def test_exponential(self):
        self.phase['pore.item1'] = 0.8e-11
        self.phase['pore.item2'] = 3
        self.phase['pore.item3'] = 0.5
        self.phase['pore.item4'] = 2
        self.phase['pore.item5'] = -0.34
        self.phase['pore.item6'] = 2e-14
        self.phase.add_model(propname='pore.source1',
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
        self.phase.update(self.alg.soln)
        self.phase.regenerate_models(propnames='pore.source1')
        X = self.phase['pore.mole_fraction']
        r1 = np.sum(0.8e-11 * 3 ** (0.5 * X[self.source_pores]**2 - 0.34) + 2e-14)
        r2 = np.sum(self.phase['pore.source1.rate'][self.source_pores])
        assert_allclose(actual=r2, desired=r1, rtol=1e-7)

    def test_natural_exponential(self):
        self.phase['pore.item1'] = 0.8e-11
        self.phase['pore.item2'] = 0.5
        self.phase['pore.item3'] = 2
        self.phase['pore.item4'] = -0.34
        self.phase['pore.item5'] = 2e-14
        self.phase.add_model(propname='pore.source1',
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
        self.phase.update(self.alg.soln)
        self.phase.regenerate_models(propnames='pore.source1')
        X = self.phase['pore.mole_fraction']
        r1 = np.sum(0.8e-11 * np.exp(0.5 * X[self.source_pores]**2 - 0.34) + 2e-14)
        r2 = np.sum(self.phase['pore.source1.rate'][self.source_pores])
        # assert_allclose(actual=r2, desired=r1)

    def test_logarithm(self):
        self.phase['pore.item1'] = 0.16e-13
        self.phase['pore.item2'] = 10
        self.phase['pore.item3'] = 4
        self.phase['pore.item4'] = 1.4
        self.phase['pore.item5'] = 0.133
        self.phase['pore.item6'] = -5.1e-13
        self.phase.add_model(propname='pore.source1',
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
        self.phase.update(self.alg.soln)
        self.phase.regenerate_models(propnames='pore.source1')
        X = self.phase['pore.mole_fraction']
        r1 = np.sum(0.16e-13 * np.log(4*X[self.source_pores]**(1.4) + 0.133)
                    / np.log(10) - 5.1e-13)
        r2 = np.sum(self.phase['pore.source1.rate'][self.source_pores])
        assert_allclose(actual=r2, desired=r1)

    def test_natural_logarithm(self):
        self.phase['pore.item1'] = 0.16e-14
        self.phase['pore.item2'] = 4
        self.phase['pore.item3'] = 1.4
        self.phase['pore.item4'] = 0.133
        self.phase['pore.item5'] = -5.1e-14
        self.phase.add_model(propname='pore.source1',
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
        self.phase.update(self.alg.soln)
        self.phase.regenerate_models(propnames='pore.source1')
        X = self.phase['pore.mole_fraction']
        r1 = np.sum(0.16e-14*np.log(4*X[self.source_pores]**1.4 + 0.133) - 5.1e-14)
        r2 = np.sum(self.phase['pore.source1.rate'][self.source_pores])
        # assert_allclose(actual=r2, desired=r1)

    def test_general_symbolic(self):
        # natural log function
        y = "a*ln(b*x**c + d)+e"
        self.phase['pore.item1'] = 0.16e-14
        self.phase['pore.item2'] = 4
        self.phase['pore.item3'] = 1.4
        self.phase['pore.item4'] = 0.133
        self.phase['pore.item5'] = -5.1e-14
        self.phase.add_model(propname='pore.source1',
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
        self.phase.add_model(propname='pore.general',
                             model=op.models.physics.source_terms.general_symbolic,
                             eqn=y, x='pore.mole_fraction', **arg_map)
        assert np.allclose(self.phase['pore.source1.rate'],
                           self.phase['pore.general.rate'])
        assert np.allclose(self.phase['pore.source1.S1'],
                           self.phase['pore.general.S1'])
        assert np.allclose(self.phase['pore.source1.S2'],
                           self.phase['pore.general.S2'])

    def test_butler_volmer_kinetics(self):
        np.random.seed(10)
        self.net["pore.reaction_area"] = np.random.rand(self.net.Np)
        self.phase['pore.electrolyte_voltage'] = np.random.rand(self.net.Np) * 0.2
        self.phase['pore.solid_voltage'] = 1.1
        self.phase['pore.open_circuit_voltage'] = 1.2
        self.phase['pore.electrolyte_concentration'] = np.random.rand(self.net.Np)
        BV_params = {
            "n": 4,
            "i0_ref": 1e-3,
            "c_ref": 1000,
            "beta": 0.5,
            "gamma": 1
        }
        self.phase.add_model(propname='pore.rxn_BV_c',
                             model=pm.source_terms.butler_volmer_conc,
                             X="pore.electrolyte_concentration", **BV_params)
        self.phase.add_model(propname='pore.rxn_BV_v',
                             model=pm.source_terms.butler_volmer_voltage,
                             X="pore.electrolyte_voltage", **BV_params)
        # Check Butler-Volmer model (concentration)
        S1_BV_c = self.phase["pore.rxn_BV_c.S1"]
        S2_BV_c = self.phase["pore.rxn_BV_c.S2"]
        rate_BV_c = self.phase["pore.rxn_BV_c.rate"]
        eps = np.finfo(float).eps * np.abs(rate_BV_c).max()
        assert_allclose(S2_BV_c.mean(), 0, atol=eps)
        assert_allclose(S1_BV_c.mean(), -0.0008190647392381356)
        assert_allclose(rate_BV_c.mean(), -0.00040556020470030826)
        # Check Butler-Volmer model (voltage)
        S1_BV_v = self.phase["pore.rxn_BV_v.S1"]
        S2_BV_v = self.phase["pore.rxn_BV_v.S2"]
        rate_BV_v = self.phase["pore.rxn_BV_v.rate"]
        assert_allclose(S2_BV_v.mean(), 2104.1093852527683)
        assert_allclose(S1_BV_v.mean(), -12190.386173792602)
        assert_allclose(rate_BV_v.mean(), -156.52244418065774)
        # The two Butler-Volmer models must only differ by n*F (unit conversion)
        assert_allclose(rate_BV_v, rate_BV_c * BV_params["n"] * 96485.33212)


if __name__ == '__main__':

    t = GenericSourceTermTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
