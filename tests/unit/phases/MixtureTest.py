import openpnm as op
import scipy as sp
import pytest


class MixtureTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[10, 10, 10])

    def test_mixture_init(self):
        ws = op.Workspace()
        ws.clear()
        self.net = op.network.Cubic(shape=[10, 10, 10])
        N2 = op.phases.components.gases.N2(network=self.net, name='pure_N2')
        O2 = op.phases.components.gases.O2(network=self.net, name='pure_O2')
        CO2 = op.phases.components.gases.CO2(network=self.net, name='pure_CO2')
        H2 = op.phases.components.gases.H2(network=self.net, name='pure_H2')
        air = op.phases.GenericMixture(network=self.net,
                                       components=[N2, O2, H2, CO2],
                                       name='air_mixture')
        assert sp.all(air['pore.mole_fraction.all'] == 0.0)
        air.set_mole_fraction(N2, 0.790)
        air.set_mole_fraction(O2, 0.209)
        air.set_mole_fraction(CO2, 0.001)
        assert sp.all(air['pore.mole_fraction.all'] == 1.0)

    def test_check_health(self):
        ws = op.Workspace()
        ws.clear()
        self.net = op.network.Cubic(shape=[10, 10, 10])
        N2 = op.phases.components.gases.N2(network=self.net, name='pure_N2')
        O2 = op.phases.components.gases.O2(network=self.net, name='pure_O2')
        CO2 = op.phases.components.gases.CO2(network=self.net, name='pure_CO2')
        H2 = op.phases.components.gases.H2(network=self.net, name='pure_H2')
        air = op.phases.GenericMixture(network=self.net,
                                       components=[N2, O2, H2, CO2],
                                       name='air_mixture')
        air.set_mole_fraction(N2, 0.790)
        air.set_mole_fraction(O2, 0.209)
        air.set_mole_fraction(CO2, 0.001)
        h = air.check_mixture_health()
        assert h.health is True
        air.set_mole_fraction(N2, 0.790)
        air.set_mole_fraction(O2, 0.209)
        air.set_mole_fraction(CO2, 0.001)
        h = air.check_mixture_health()
        assert h.health is True
        air.set_mole_fraction(CO2, 0.002)
        h = air.check_mixture_health()
        assert h.health is False
        air.set_mole_fraction(CO2, 0.000)
        h = air.check_mixture_health()
        assert h.health is False
        assert len(h['mole_fraction_too_low']) == air.Np
        air.set_mole_fraction(CO2, 0.001)
        air['pore.mole_fraction.'+CO2.name][0] = 0.0
        h = air.check_mixture_health()
        assert h.health is False
        assert len(h['mole_fraction_too_low']) == 1


if __name__ == '__main__':

    t = MixtureTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
