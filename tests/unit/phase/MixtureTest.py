import openpnm as op
import numpy as np


class MixtureTest:
    def setup_class(self):
        ws = op.Workspace()
        ws.clear()

    def test_init_binary_gas(self):
        net = op.network.Demo()
        N2 = op.phase.Species(network=net, species='n2', name='pure_N2')
        O2 = op.phase.Species(network=net, species='o2', name='pure_O2')
        air = op.phase.BinaryGas(network=net, components=[N2, O2], name='air')

    def test_liquid_mixture_two_comps(self):
        net = op.network.Demo()
        A = op.phase.StandardLiquid(network=net, species='ethanol')
        B = op.phase.StandardLiquid(network=net, species='h2o')
        vodka = op.phase.LiquidMixture(network=net, components=[A, B])
        vodka.x(A.name, 0.4)
        vodka.x(B.name, 0.6)
        assert len(vodka.components) == 2

    def test_standard_liquid_mixture(self):
        net = op.network.Demo()
        A = op.phase.StandardLiquid(network=net, species='ethanol')
        A.regenerate_models()
        B = op.phase.StandardLiquid(network=net, species='h2o')
        B.regenerate_models()
        vodka = op.phase.StandardLiquidMixture(network=net, components=[A, B])
        vodka.x(A.name, 0.4)
        vodka.x(B.name, 0.6)
        vodka.regenerate_models()


    def test_add_and_remove_component_method(self):
        net = op.network.Demo()
        o2 = op.phase.Species(network=net, species='o2', name='pure_O2')
        n2 = op.phase.Species(network=net, species='n2', name='pure_N2')
        air = op.phase.GasMixture(network=net, components=[n2, o2])
        air.remove_comp(n2)
        assert len(air.components) == 1
        air.remove_comp(o2)
        assert len(air.components) == 0

    def test_check_health(self):
        net = op.network.Demo()
        o2 = op.phase.Species(network=net, species='o2', name='pure_O2')
        n2 = op.phase.Species(network=net, species='n2', name='pure_N2')
        air = op.phase.GasMixture(network=net, components=[n2, o2])
        air['pore.mole_fraction.pure_N2'] = 0.79
        air['pore.mole_fraction.pure_O2'] = 0.21
        h = air.check_mixture_health()
        assert h.health is True
        air['pore.mole_fraction.pure_O2'] = 0.25
        h = air.check_mixture_health()
        assert h.health is False

    def test_getitem(self):
        net = op.network.Demo()
        o2 = op.phase.Species(network=net, species='o2', name='pure_O2')
        n2 = op.phase.Species(network=net, species='n2', name='pure_N2')
        air = op.phase.GasMixture(network=net, components=[n2, o2])
        d = air['pore.mole_fraction']
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
