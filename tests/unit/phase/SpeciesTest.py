import openpnm as op
from openpnm.phase import


class SpeciesTest:
    def setup_class(self):
        ws = op.Workspace()
        ws.clear()
        self.net = op.network.Cubic(shape=[10, 10, 10])
        self.net = op.network.Cubic(shape=[10, 10, 10])
        self.N2 = mixtures.species.gases.N2(network=self.net, name='pure_N2')
        self.O2 = mixtures.species.gases.O2(network=self.net, name='pure_O2')
        self.CO2 = mixtures.species.gases.CO2(network=self.net, name='pure_CO2')
        self.H2 = mixtures.species.gases.H2(network=self.net, name='pure_H2')
        self.air = mixtures.GenericMixture(network=self.net,
                                           components=[self.N2, self.O2,
                                                       self.H2, self.CO2],
                                           name='air_mixture')

    def test_lookup_mixture(self):
        assert self.N2.mixture is self.air


if __name__ == '__main__':

    t = SpeciesTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
