import openpnm as op
import numpy as np


class MixtureTest:
    def setup_class(self):
        ws = op.Workspace()
        ws.clear()
        self.net = op.network.Cubic(shape=[10, 10, 10])
        self.net = op.network.Cubic(shape=[10, 10, 10])
        self.N2 = op.phase.GasByName(network=self.net,
                                     species='n2', name='pure_N2')
        self.O2 = op.phase.GasByName(network=self.net,
                                     species='o2', name='pure_O2')
        self.air = op.phase.GasMixture(network=self.net,
                                       components=[self.N2, self.O2],
                                       name='air_mixture')

    def test_set_component(self):
        self.CO2 = op.phase.GasByName(network=self.net, species='co2',
                                      name='pure_CO2')
        self.air.components = self.CO2
        self.air['pore.mole_fraction.pure_N2'] = 0.79
        self.air['pore.mole_fraction.pure_O2'] = 0.20
        self.air['pore.mole_fraction.pure_CO2'] = 0.01
        conc = np.zeros(self.air.Np)
        for i in self.air['pore.mole_fraction'].values():
            conc += i
        assert np.all(conc == 1)
        assert len(self.air.components) == 3
        del self.air['pore.mole_fraction.pure_CO2']
        assert len(self.air.components) == 2

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
