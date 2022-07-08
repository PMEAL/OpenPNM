import openpnm as op
import numpy as np


class MixtureTest:
    def setup_class(self):
        ws = op.Workspace()
        ws.clear()

    def test_init_binary_gas(self):
        net = op.network.Demo()
        N2 = op.phase.GasByName(network=net, species='n2', name='pure_N2')
        O2 = op.phase.GasByName(network=net, species='o2', name='pure_O2')
        air = op.phase.BinaryGas(network=net, components=[N2, O2], name='air')

    def test_set_component_by_property(self):
        net = op.network.Demo()

    def test_liquid_mixture_2comps(self):
        net = op.network.Demo()
        A = op.phase.LiquidByName(network=net, species='ethanol')
        B = op.phase.LiquidByName(network=net, species='h2o')
        vodka = op.phase.LiquidMixture(network=net, components=[A, B])
        vodka.x(A.name, 0.4)
        vodka.x(B.name, 0.6)
        vodka.regenerate_models()
        # Make sure the density is between the two components
        assert vodka['pore.density'].mean() > A['pore.density'].mean()
        assert vodka['pore.density'].mean() < B['pore.density'].mean()


    def test_liquid_mixture_3comps(self):
        pass

    def test_add_and_remove_component_method(self):
        net = op.network.Demo()
        o2 = op.phase.GasByName(network=net, species='o2', name='pure_O2')
        n2 = op.phase.GasByName(network=net, species='n2', name='pure_N2')
        air = op.phase.GasMixture(network=net, components=[n2, o2])
        air.remove_comp(n2)
        air.remove_comp(o2)

    def test_check_health(self):
        self.air['pore.mole_fraction.pure_N2'] = 0.79
        self.air['pore.mole_fraction.pure_O2'] = 0.21
        h = self.air.check_mixture_health()
        assert h.health is True
        self.air['pore.mole_fraction.pure_O2'] = 0.25
        h = self.air.check_mixture_health()
        assert h.health is False

    def test_getitem(self):
        d = self.air['pore.mole_fraction']
        set_a = set(['pure_N2', 'pure_O2'])
        assert set_a.difference(set(d.keys())) == set()


if __name__ == '__main__':

    t = MixtureTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
