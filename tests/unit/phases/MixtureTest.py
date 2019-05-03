import openpnm as op
import scipy as sp
import pytest


class MixtureTest:
    def setup_class(self):
        ws = op.Workspace()
        ws.clear()
        self.net = op.network.Cubic(shape=[10, 10, 10])
        self.net = op.network.Cubic(shape=[10, 10, 10])
        self.N2 = op.phases.species.gases.N2(network=self.net, name='pure_N2')
        self.O2 = op.phases.species.gases.O2(network=self.net, name='pure_O2')
        self.CO2 = op.phases.species.gases.CO2(network=self.net, name='pure_CO2')
        self.H2 = op.phases.species.gases.H2(network=self.net, name='pure_H2')
        self.air = op.phases.mixtures.GenericMixture(network=self.net,
                                                     components=[self.N2,
                                                                 self.O2,
                                                                 self.H2,
                                                                 self.CO2],
                                                     name='air_mixture')

    def test_set_mole_fraction(self):
        self.air.set_mole_fraction(self.N2, 0.790)
        self.air.set_mole_fraction(self.O2, 0.209)
        self.air.set_mole_fraction(self.CO2, 0.001)
        assert sp.all(self.air['pore.mole_fraction.all'] == 1.0)

    def test_check_health(self):
        self.air.set_mole_fraction(self.N2, 0.790)
        self.air.set_mole_fraction(self.O2, 0.209)
        self.air.set_mole_fraction(self.CO2, 0.001)
        h = self.air.check_mixture_health()
        assert h.health is True
        self.air.set_mole_fraction(self.CO2, 0.002)
        h = self.air.check_mixture_health()
        assert h.health is False
        self.air.set_mole_fraction(self.CO2, 0.000)
        h = self.air.check_mixture_health()
        assert h.health is False
        assert len(h['mole_fraction_too_low']) == self.air.Np
        self.air.set_mole_fraction(self.CO2, 0.001)
        self.air['pore.mole_fraction.'+self.CO2.name][0] = 0.0
        h = self.air.check_mixture_health()
        assert h.health is False
        assert len(h['mole_fraction_too_low']) == 1

    def test_getitem(self):
        d = self.air['pore.mole_fraction']
        set_a = set(['pore.mole_fraction.pure_N2', 'pore.mole_fraction.pure_O2',
                     'pore.mole_fraction.pure_H2', 'pore.mole_fraction.pure_CO2',
                     'pore.mole_fraction.all'])
        assert set_a.difference(set(d.keys())) == set()


if __name__ == '__main__':

    t = MixtureTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
